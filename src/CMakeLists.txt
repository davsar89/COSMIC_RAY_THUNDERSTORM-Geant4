#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MOS_test C CXX Fortran)
enable_language (Fortran)

execute_process(COMMAND "source ~/.bashrc")

set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")




#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED gdml ui_all vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()



# CADMesh
#find_package(cadmesh)
#include_directories(/Home/siv29/dsa030/Desktop/GEANT4/CADMesh/install/include/)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB sourcesf ${PROJECT_SOURCE_DIR}/src/*.f)
file(GLOB sourcesf2 ${PROJECT_SOURCE_DIR}/src/msis/*.for)
file(GLOB sourcesf90 ${PROJECT_SOURCE_DIR}/src/parma/*.f90)
file(GLOB sourcesc ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB sourcesc2 ${PROJECT_SOURCE_DIR}/src/C_functions/*.c)
file(GLOB sourcesh ${PROJECT_SOURCE_DIR}/include/*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(mos_test mos_test.cc ${sourcesf} ${sourcesf2} ${sourcesf90} ${sourcesc} ${sourcesc2} ${sourcesh} ${sources} ${headers})
target_link_libraries(mos_test ${Geant4_LIBRARIES} )
#target_link_libraries(fegs_bgo ${cadmesh_LIBRARIES} )

target_compile_options(mos_test PRIVATE -std=c++11 -mfpmath=sse -msse2 -ffp-contract=off)
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS mos_test DESTINATION bin)

