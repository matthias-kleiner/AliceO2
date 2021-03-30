// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TPCSimulation/IDCSim.h"
#include "TFile.h"
#include "TPCBase/Mapper.h"

void o2::tpc::IDCSim::integrateDigitsForOneTF(const gsl::span<const o2::tpc::Digit>& digits)
{
  resetIDCs();

  // loop over digits from one sector for ALL Time Frames
  const unsigned int switchAfterTB = getLastTimeBinForSwitch();
  for (const auto& digit : digits) {
    const o2::tpc::CRU cru(digit.getCRU());
    const unsigned int region = cru.region();
    const int timeStamp = digit.getTimeStamp();
    if (timeStamp < switchAfterTB) {
      const unsigned int indexPad = getIndex(timeStamp, region, digit.getRow(), digit.getPad());
      mIDCsTmp[mBufferIndex][region][indexPad] += digit.getChargeFloat();
    } else {
      const unsigned int indexPad = getIndex(timeStamp - switchAfterTB, region, digit.getRow(), digit.getPad());
      mIDCsTmp[!mBufferIndex][region][indexPad] += digit.getChargeFloat();
    }
  }

  mBufferIndex = !mBufferIndex;  // switch buffer index
  mTimeBinsOff = getNewOffset(); // set offset
}

void o2::tpc::IDCSim::dumpIDCs(const int timeframe)
{
  const std::string name = fmt::format("idcs_obj_{:02}_{:02}.root", mSector, timeframe);
  TFile fOut(name.data(), "RECREATE");
  int cru = mSector * mRegions;
  for (const auto& idcs : mIDCsTmp[!mBufferIndex]) {
    fOut.WriteObject(&idcs, Form("cru_%i", cru));
    ++cru;
  }
  fOut.Close();
}

void o2::tpc::IDCSim::createDebugTree(const int timeframe)
{
  const static Mapper& mapper = Mapper::instance();

  LOGP(info, "mSector debug {}", mSector);
  const std::string nameTree = fmt::format("idcs_tree_{:02}_{:02}.root", mSector, timeframe);
  o2::utils::TreeStreamRedirector pcstream(nameTree.data(), "RECREATE");
  pcstream.GetFile()->cd();

  int cru = mSector * mRegions;
  for (const auto& idcs : mIDCsTmp[!mBufferIndex]) {
    int sectorTmp = mSector;
    const o2::tpc::CRU cruTmp(cru);
    unsigned int region = cruTmp.region();
    const int padsPerCRU = mPadsPerRegion[region];
    std::vector<int> vRow(padsPerCRU);
    std::vector<int> vPad(padsPerCRU);
    std::vector<float> vXPos(padsPerCRU);
    std::vector<float> vYPos(padsPerCRU);
    std::vector<float> vGlobalXPos(padsPerCRU);
    std::vector<float> vGlobalYPos(padsPerCRU);
    std::vector<float> idcsPerTimeBin(padsPerCRU); // idcs for one time bin

    for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
      const GlobalPadNumber globalNum = mGlobalPadOffs[region] + iPad;
      const auto& padPosLocal = mapper.padPos(globalNum);
      vRow[iPad] = padPosLocal.getRow();
      vPad[iPad] = padPosLocal.getPad();
      vXPos[iPad] = mapper.getPadCentre(padPosLocal).X();
      vYPos[iPad] = mapper.getPadCentre(padPosLocal).Y();

      const GlobalPosition2D globalPos = mapper.LocalToGlobal(LocalPosition2D(vXPos[iPad], vYPos[iPad]), mSector);
      vGlobalXPos[iPad] = globalPos.X();
      vGlobalYPos[iPad] = globalPos.Y();
    }
    int cruiTmp = cru;

    for (int iTimeBin = 0; iTimeBin < mIntegrationIntervalsPerTF; ++iTimeBin) {
      for (int iPad = 0; iPad < padsPerCRU; ++iPad) {
        idcsPerTimeBin[iPad] = (idcs)[iPad + iTimeBin * mPadsPerRegion[region]];
      }

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
