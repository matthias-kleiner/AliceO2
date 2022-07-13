// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file DebugStreamer.h
/// \brief Definition of class for writing debug informations
/// \author Matthias Kleiner <mkleiner@ikf.uni-frankfurt.de>

#ifndef ALICEO2_TPC_DEBUGSTREAMER_H_
#define ALICEO2_TPC_DEBUGSTREAMER_H_

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
#include "CommonUtils/ConfigurableParamHelper.h"
#include "CommonUtils/TreeStreamRedirector.h"
#endif

namespace o2::gpu
{

#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
/// struct for setting and storing the streamer level
struct ParameterDebugStreamer : public o2::conf::ConfigurableParamHelper<ParameterDebugStreamer> {

  /// struct defining the flags which can be used to check if a certain debug streamer is used
  enum StreamFlags {
    streamdEdx = 1 << 0, ///< stream corrections and cluster properties used for the dE/dx
  };

  StreamFlags StreamLevel{}; /// flag to store what will be streamed
  O2ParamDef(ParameterDebugStreamer, "DebugStreamerParam");
};

inline ParameterDebugStreamer::StreamFlags operator|(ParameterDebugStreamer::StreamFlags a, ParameterDebugStreamer::StreamFlags b) { return static_cast<ParameterDebugStreamer::StreamFlags>(static_cast<int>(a) | static_cast<int>(b)); }
inline ParameterDebugStreamer::StreamFlags operator&(ParameterDebugStreamer::StreamFlags a, ParameterDebugStreamer::StreamFlags b) { return static_cast<ParameterDebugStreamer::StreamFlags>(static_cast<int>(a) & static_cast<int>(b)); }
inline ParameterDebugStreamer::StreamFlags operator~(ParameterDebugStreamer::StreamFlags a) { return static_cast<ParameterDebugStreamer::StreamFlags>(~static_cast<int>(a)); }
#endif

/// class to enable streaming debug information to root files
class DebugStreamer
{
#if !defined(GPUCA_GPUCODE) && !defined(GPUCA_STANDALONE)
 public:
  /// constructor
  /// \param outFile output file name without .root suffix
  /// \param option RECREATE or UPDATE
  DebugStreamer(const char* outFile, const char* option);

  /// \return returns streamer object
  o2::utils::TreeStreamRedirector& stream() { return *mTreeStreamer; }

  /// \return returns streamer level i.e. what will be written to file
  static ParameterDebugStreamer::StreamFlags getStreamFlags() { return ParameterDebugStreamer::Instance().StreamLevel; }

  ///< return returns unique ID for each CPU thread to give each thread an own output file
  static size_t getCPUID();

  /// \return returns number of trees in the streamer
  int getNTrees() const;

  /// \return returns an unique branch name which is not already written in the file
  /// \param tree name of the tree for which to get a unique branch name
  std::string getUniqueTreeName(const char* tree) const;

  /// set directly the debug level
  static void setStreamFlags(const ParameterDebugStreamer::StreamFlags streamFlags) { o2::conf::ConfigurableParam::setValue("DebugStreamerParam", "StreamLevel", static_cast<int>(streamFlags)); }

  /// enable specific streamer flag
  static void enableStream(const ParameterDebugStreamer::StreamFlags streamFlag);

  /// disable a specific streamer flag
  static void disableStream(const ParameterDebugStreamer::StreamFlags streamFlag);

  /// check if streamer for specific flag is enabled
  static bool checkStream(const ParameterDebugStreamer::StreamFlags streamFlag) { return ((getStreamFlags() & streamFlag) == streamFlag); }

  /// merge trees with the same content structure, but different naming
  /// \param inpFile input file containing several trees with the same content
  /// \param outFile contains the merged tree from the input file in one branch
  /// \param option setting which is used for the merging
  static void mergeTrees(const char* inpFile, const char* outFile, const char* option = "fast");

 private:
  std::unique_ptr<o2::utils::TreeStreamRedirector> mTreeStreamer; ///< streamer which is used for the debugging
#else
  /* Empty for GPU */
#endif
};

} // namespace o2::gpu

#endif
