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

  mBufferIndex = !mBufferIndex;  // switch buffer index
  mTimeBinsOff = getNewOffset(); // set offset
}

unsigned int o2::tpc::IDCSim::getLastTimeBinForSwitch() const
{
  const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsReminder);
  return totaloffs >= mTimeStampsPerIntegrationInterval ? mIntegrationIntervalsPerTF * mTimeStampsPerIntegrationInterval : (mIntegrationIntervalsPerTF - mAddInterval) * mTimeStampsPerIntegrationInterval - mTimeBinsOff;
}

int o2::tpc::IDCSim::getNewOffset() const
{
  const int totaloffs = mTimeBinsOff + static_cast<int>(mTimeStampsReminder);
  return totaloffs >= mTimeStampsPerIntegrationInterval ? (totaloffs - static_cast<int>(mTimeStampsPerIntegrationInterval)) : totaloffs;
}

/// set all IDC values to 0
void o2::tpc::IDCSim::resetIDCs()
{
  for (auto& idcs : mIDCs[!mBufferIndex]) {
    std::fill(idcs.begin(), idcs.end(), 0);
  }
}

void o2::tpc::IDCSim::dumpIDCs(const int timeframe)
{
  const std::string name = fmt::format("idcs_obj_{:02}_{:02}.root", mSector, timeframe);
  TFile fOut(name.data(), "RECREATE");
  int cru = mSector * Mapper::NREGIONS;
  for (const auto& idcs : mIDCs[!mBufferIndex]) {
    fOut.WriteObject(&idcs, Form("cru_%i", cru));
    ++cru;
  }
  fOut.Close();
}

void o2::tpc::IDCSim::createDebugTree(const int timeframe)
{
  const static Mapper& mapper = Mapper::instance();

  const std::string nameTree = fmt::format("idcs_tree_{:02}_{:02}.root", mSector, timeframe);
  o2::utils::TreeStreamRedirector pcstream(nameTree.data(), "RECREATE");
  pcstream.GetFile()->cd();

  int cru = mSector * Mapper::NREGIONS;
  for (const auto& idcs : mIDCs[!mBufferIndex]) {
    int sectorTmp = mSector;
    const o2::tpc::CRU cruTmp(cru);
    unsigned int region = cruTmp.region();
    const int padsPerCRU = Mapper::PADSPERREGION[region];
    std::vector<int> vRow(padsPerCRU);
    std::vector<int> vPad(padsPerCRU);
    std::vector<float> vXPos(padsPerCRU);
    std::vector<float> vYPos(padsPerCRU);
    std::vector<float> vGlobalXPos(padsPerCRU);
    std::vector<float> vGlobalYPos(padsPerCRU);
    std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

    for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
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

    for (int iTimeBin = 0; iTimeBin < mIntegrationIntervalsPerTF; ++iTimeBin) {
      for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
        idcsPerTimeBin[iPad] = (idcs)[iPad + iTimeBin * Mapper::PADSPERREGION[region]];
      }
      int cruiTmp = cru;
      pcstream << "tree"
               << "cru=" << cruiTmp
               << "sector=" << sectorTmp
               << "region=" << region
               << "timeBin=" << iTimeBin
               << "IDCs.=" << idcsPerTimeBin
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
