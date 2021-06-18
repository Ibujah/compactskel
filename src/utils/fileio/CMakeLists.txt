set(LIBRARY_NAME ${FILEIO_LIB})

set(SOURCE_FILES SkeletonFile.cpp
                 BoundaryFile2D.cpp)

# make the library
add_library(
    ${LIBRARY_NAME}
    ${INCLUDE_FILES}
    ${SOURCE_FILES}
	)

include_directories(${CMAKE_SOURCE_DIR}/src/lib/
					${Boost_INCLUDE_DIR})

target_link_libraries(${LIBRARY_NAME} ${MATHTOOLS_LIB}
									  ${SHAPE_LIB}
									  ${BOUNDARY_LIB}
									  ${SKELETON_LIBS})


SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE CXX)
