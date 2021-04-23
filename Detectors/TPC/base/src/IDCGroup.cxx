// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "TPCBase/IDCGroup.h"
#include "CommonUtils/TreeStreamRedirector.h" // for debugging
#include "TFile.h"

void o2::tpc::IDCGroup::initStorage()
{
  mNIDCsPerCRU = std::accumulate(mPadsPerRow.begin(), mPadsPerRow.end(), decltype(mPadsPerRow)::value_type(0));
  mIDCsGrouped.resize(mNIDCsPerCRU);

  for (unsigned int i = 1; i < mRows; ++i) {
    const unsigned int lastInd = i - 1;
    mOffsRow[i] = mOffsRow[lastInd] + mPadsPerRow[lastInd];
  }
}

void o2::tpc::IDCGroup::dumpToTree(const char* outname) const
{
  o2::utils::TreeStreamRedirector pcstream(outname, "RECREATE");
  pcstream.GetFile()->cd();
  for (int integrationInterval = 0; integrationInterval < getNIntegrationIntervals(); ++integrationInterval) {
    for (int irow = 0; irow < mRows; ++irow) {
      for (int ipad = 0; ipad < mPadsPerRow[irow]; ++ipad) {
        float idc = (*this)(irow, ipad, integrationInterval);
        pcstream << "idcs"
                 << "row=" << irow
                 << "pad=" << ipad
                 << "IDC=" << idc
                 << "\n";
      }
    }
  }
  pcstream.Close();
}
