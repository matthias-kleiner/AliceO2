O2_SETUP(NAME MCHMappingSegContour)
set(BUCKET_NAME mch_mapping_segcontour_bucket)

set(SRCS
        src/CathodeSegmentationContours.cxx
        src/CathodeSegmentationSVGWriter.cxx
        src/SegmentationContours.cxx)

set(LIBRARY_NAME MCHMappingSegContour3)

O2_GENERATE_LIBRARY()

O2_GENERATE_EXECUTABLE(
        EXE_NAME o2-mch-mapping-svg-segmentation3
        SOURCES src/SVGSegmentation.cxx
        BUCKET_NAME mch_mapping_segcontour_bucket
        MODULE_LIBRARY_NAME ${LIBRARY_NAME})

