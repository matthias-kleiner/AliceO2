set(MODULE_NAME "TPCSimulation")

O2_SETUP(NAME ${MODULE_NAME})

link_directories( ${LINK_DIRECTORIES})

set(SRCS
   src/Detector.cxx
   src/DigitMCMetaData.cxx
   src/DigitContainer.cxx
   src/DigitGlobalPad.cxx
   src/Digitizer.cxx
   src/DigitTime.cxx
   src/ElectronTransport.cxx
   src/GEMAmplification.cxx
   src/PadResponse.cxx
   src/Point.cxx
   src/SAMPAProcessing.cxx
   src/SpaceCharge.cxx
)

set(HEADERS
   include/${MODULE_NAME}/Detector.h
   include/${MODULE_NAME}/DigitMCMetaData.h
   include/${MODULE_NAME}/DigitContainer.h
   include/${MODULE_NAME}/DigitGlobalPad.h
   include/${MODULE_NAME}/Digitizer.h
   include/${MODULE_NAME}/DigitTime.h
   include/${MODULE_NAME}/ElectronTransport.h
   include/${MODULE_NAME}/GEMAmplification.h
   include/${MODULE_NAME}/PadResponse.h
   include/${MODULE_NAME}/Point.h
   include/${MODULE_NAME}/SAMPAProcessing.h
   include/${MODULE_NAME}/SpaceCharge.h
)
Set(LINKDEF src/TPCSimulationLinkDef.h)
Set(LIBRARY_NAME ${MODULE_NAME})
set(BUCKET_NAME tpc_simulation_bucket)

O2_GENERATE_LIBRARY()

install(
  DIRECTORY files
  DESTINATION share/Detectors/TPC/
)

# todo we repeat ourselves because the macro O2_GENERATE_LIBRARY dares deleting the variables we pass to it.
set(BUCKET_NAME tpc_simulation_bucket)
set(LIBRARY_NAME ${MODULE_NAME})

#list(LENGTH Exe_Names _length)
#math(EXPR _length ${_length}-1)

#foreach (_file RANGE 0 ${_length}) # loop over a range because we traverse 2 lists and not 1
  #list(GET Exe_Names ${_file} _name)
  #list(GET Exe_Source ${_file} _src)
  ## FIXME: propably a copy paste remnant, remove?
##  Set(DEPENDENCIES CCDB)
  #O2_GENERATE_EXECUTABLE(
      #EXE_NAME ${_name}
      #SOURCES ${_src}
      #MODULE_LIBRARY_NAME ${LIBRARY_NAME}
      #BUCKET_NAME ${BUCKET_NAME}
  #)
#endforeach (_file RANGE 0 ${_length})

set(TEST_SRCS
   test/testTPCDigitContainer.cxx
   test/testTPCElectronTransport.cxx
   test/testTPCGEMAmplification.cxx
   test/testTPCSAMPAProcessing.cxx
   test/testTPCSimulation.cxx
)

O2_GENERATE_TESTS(
  BUCKET_NAME ${BUCKET_NAME}
  MODULE_LIBRARY_NAME ${MODULE_NAME}
  TEST_SRCS ${TEST_SRCS}
)

# add the TPC run sim as a unit test (if simulation was enabled)
if (HAVESIMULATION)
  add_test_wrap(NAME tpcsim_G4 COMMAND ${CMAKE_BINARY_DIR}/bin/o2-sim-tpc -n 2 -e TGeant4)
  set_tests_properties(tpcsim_G4 PROPERTIES PASS_REGULAR_EXPRESSION "Macro finished succesfully")
  add_test_wrap(NAME tpcsim_G3 COMMAND ${CMAKE_BINARY_DIR}/bin/o2-sim-tpc -n 2 -e TGeant3)
  set_tests_properties(tpcsim_G3 PROPERTIES PASS_REGULAR_EXPRESSION "Macro finished succesfully")
  # sets the necessary environment
  set_tests_properties(tpcsim_G3 tpcsim_G4  PROPERTIES ENVIRONMENT VMCWORKDIR=${CMAKE_SOURCE_DIR})
endif()
