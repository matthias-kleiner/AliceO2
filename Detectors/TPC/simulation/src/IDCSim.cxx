// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "TPCSimulation/IDCSim.h"
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TFile.h"
#include "TPCBase/Mapper.h"
#include <fmt/format.h>
#include "Framework/Logger.h"

void o2::tpc::IDCSim::integrateDigitsForOneTF(const gsl::span<const o2::tpc::Digit>& digits)
{
  resetIDCs();

  // loop over digits from one sector for ALL Time Frames
  const unsigned int switchAfterTB = getLastTimeBinForSwitch();

  if (mAddInterval) {
    // decrease the size of the vector if the last integration intervall is empty
    if (switchAfterTB == (mIntegrationIntervalsPerTF - 1) * mTimeStampsPerIntegrationInterval) {
      for (unsigned int ireg = 0; ireg < Mapper::NREGIONS; ++ireg) {
        mIDCs[mBufferIndex][ireg].resize(mMaxIDCs[ireg] - Mapper::PADSPERREGION[ireg]);
      }
    } else {
      for (auto& idcs : mIDCs[mBufferIndex]) {
        idcs.resize(idcs.capacity());
      }
    }
  }

  for (const auto& digit : digits) {
    const o2::tpc::CRU cru(digit.getCRU());
    const unsigned int region = cru.region();
    const int timeStamp = digit.getTimeStamp();
    if (timeStamp < switchAfterTB) {
      const unsigned int indexPad = getIndex(timeStamp, region, digit.getRow(), digit.getPad());
      mIDCs[mBufferIndex][region][indexPad] += digit.getChargeFloat();
    } else {
      const unsigned int indexPad = getIndex(timeStamp - switchAfterTB, region, digit.getRow(), digit.getPad());
      mIDCs[!mBufferIndex][region][indexPad] += digit.getChargeFloat();
    }
  }

  mBufferIndex = !mBufferIndex; // switch buffer index
  setNewOffset();               // set offset
}

unsigned int o2::tpc::IDCSim::getLastTimeBinForSwitch() const
{
  const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsRemainder);
  return (totaloffs >= mTimeStampsPerIntegrationInterval) ? mIntegrationIntervalsPerTF * mTimeStampsPerIntegrationInterval - mTimeBinsOff : (mIntegrationIntervalsPerTF - mAddInterval) * mTimeStampsPerIntegrationInterval - mTimeBinsOff;
}

void o2::tpc::IDCSim::setNewOffset()
{
  const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsRemainder);
  mTimeBinsOff = (totaloffs >= mTimeStampsPerIntegrationInterval) ? (totaloffs - static_cast<int>(mTimeStampsPerIntegrationInterval)) : totaloffs;
}

/// set all IDC values to 0
void o2::tpc::IDCSim::resetIDCs()
{
  for (auto& idcs : mIDCs[!mBufferIndex]) {
    std::fill(idcs.begin(), idcs.end(), 0);
  }
}

void o2::tpc::IDCSim::dumpIDCs(const char* outFileName, const char* outName) const
{
  TFile fOut(outFileName, "RECREATE");
  fOut.WriteObject(this, outName);
  fOut.Close();
}

void o2::tpc::IDCSim::createDebugTree(const char* nameTree) const
{
  const Mapper& mapper = Mapper::instance();
  o2::utils::TreeStreamRedirector pcstream(nameTree, "RECREATE");
  pcstream.GetFile()->cd();

  int cru = mSector * Mapper::NREGIONS;

  // loop over data from regions
  for (const auto& idcs : mIDCs[!mBufferIndex]) {
    unsigned int sectorTmp = mSector;
    const o2::tpc::CRU cruTmp(cru);
    unsigned int region = cruTmp.region();
    const unsigned long padsPerCRU = Mapper::PADSPERREGION[region];
    std::vector<int> vRow(padsPerCRU);
    std::vector<int> vPad(padsPerCRU);
    std::vector<float> vXPos(padsPerCRU);
    std::vector<float> vYPos(padsPerCRU);
    std::vector<float> vGlobalXPos(padsPerCRU);
    std::vector<float> vGlobalYPos(padsPerCRU);
    std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

    for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
      const GlobalPadNumber globalNum = Mapper::GLOBALPADOFFSET[region] + iPad;
      const auto& padPosLocal = mapper.padPos(globalNum);
      vRow[iPad] = padPosLocal.getRow();
      vPad[iPad] = padPosLocal.getPad();
      vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
      vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();

      const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), cruTmp.sector());
      vGlobalXPos[iPad] = globalPos.X();
      vGlobalYPos[iPad] = globalPos.Y();
    }

    for (unsigned int integrationInterval = 0; integrationInterval < mIntegrationIntervalsPerTF; ++integrationInterval) {
      for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
        idcsPerTimeBin[iPad] = (idcs)[iPad + integrationInterval * Mapper::PADSPERREGION[region]];
      }

      pcstream << "tree"
               << "cru=" << cru
               << "sector=" << sectorTmp
               << "region=" << region
               << "integrationInterval=" << integrationInterval
               << "IDC.=" << idcsPerTimeBin
               << "pad.=" << vPad
               << "row.=" << vRow
               << "lx.=" << vXPos
               << "ly.=" << vYPos
               << "gx.=" << vGlobalXPos
               << "gy.=" << vGlobalYPos
               << "\n";
    }
    ++cru;
  }
  pcstream.Close();
}

void o2::tpc::IDCSim::createDebugTreeForAllCRUs(const char* nameTree, const char* filename)
{
  const Mapper& mapper = Mapper::instance();
  o2::utils::TreeStreamRedirector pcstream(nameTree, "RECREATE");
  pcstream.GetFile()->cd();

  TFile fInp(filename, "READ");
  for (TObject* keyAsObj : *fInp.GetListOfKeys()) {
    const auto key = dynamic_cast<TKey*>(keyAsObj);
    LOGP(info, "Key name: {} Type: {}", key->GetName(), key->GetClassName());

    if (std::strcmp(o2::tpc::IDCSim::Class()->GetName(), key->GetClassName()) != 0) {
      LOGP(info, "skipping object. wrong class.");
      continue;
    }

    IDCSim* idcsim = (IDCSim*)fInp.Get(key->GetName());
    const unsigned int sector = idcsim->getSector();
    unsigned int cru = sector * Mapper::NREGIONS;

    // loop over data from regions
    for (const auto& idcs : idcsim->get()) {
      int sectorTmp = sector;
      const o2::tpc::CRU cruTmp(cru);
      unsigned int region = cruTmp.region();
      const unsigned long padsPerCRU = Mapper::PADSPERREGION[region];
      std::vector<int> vRow(padsPerCRU);
      std::vector<int> vPad(padsPerCRU);
      std::vector<float> vXPos(padsPerCRU);
      std::vector<float> vYPos(padsPerCRU);
      std::vector<float> vGlobalXPos(padsPerCRU);
      std::vector<float> vGlobalYPos(padsPerCRU);
      std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

      for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
        const GlobalPadNumber globalNum = Mapper::GLOBALPADOFFSET[region] + iPad;
        const auto& padPosLocal = mapper.padPos(globalNum);
        vRow[iPad] = padPosLocal.getRow();
        vPad[iPad] = padPosLocal.getPad();
        vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
        vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();

        const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), cruTmp.sector());
        vGlobalXPos[iPad] = globalPos.X();
        vGlobalYPos[iPad] = globalPos.Y();
      }

      for (unsigned int integrationInterval = 0; integrationInterval < idcsim->getNIntegrationIntervalsPerTF(); ++integrationInterval) {
        for (unsigned int iPad = 0; iPad < padsPerCRU; ++iPad) {
          idcsPerTimeBin[iPad] = (idcs)[iPad + integrationInterval * Mapper::PADSPERREGION[region]];
        }

        pcstream << "tree"
                 << "cru=" << cru
                 << "sector=" << sectorTmp
                 << "region=" << region
                 << "integrationInterval=" << integrationInterval
                 << "IDC.=" << idcsPerTimeBin
                 << "pad.=" << vPad
                 << "row.=" << vRow
                 << "lx.=" << vXPos
                 << "ly.=" << vYPos
                 << "gx.=" << vGlobalXPos
                 << "gy.=" << vGlobalYPos
                 << "\n";
      }
      ++cru;
    }
    delete idcsim;
  }
  pcstream.Close();
}
