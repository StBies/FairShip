# Create a library called "libDtAlignmentTools" which includes the source files given in

# Set the build type. Possibilities are None, Debug, Release,
# RelWithDebInfo and MinSizeRel
SET(CMAKE_BUILD_TYPE Release)

set(INCLUDE_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/charmdet/DtAlignment/cpp_utils/include
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
) 

link_directories( ${LINK_DIRECTORIES})

# add the list of source code files -- careful, when a new file is added!
AUX_SOURCE_DIRECTORY( ${CMAKE_SOURCE_DIR}/charmdet/DtAlignment/cpp_utils/src	library_sources )
set(OWN_DIR ${CMAKE_SOURCE_DIR}/charmdet/DtAlignment/cpp_utils) 
set(SRCS ${OWN_DIR}/src/track_fit.cpp)

# set(HEADERS )
Set(LINKDEF ${CMAKE_SOURCE_DIR}/charmdet/DtAlignment/cpp_utils/src/linfit_linkDef.h)
Set(LIBRARY_NAME DtAlignmentTools)
Set(DEPENDENCIES Core Matrix MathCore)

GENERATE_LIBRARY()