project(RotPriorH-Examples CXX)


# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})

add_library(exp_help ../utils/generatePointCloudPlanar.h ../utils/generatePointCloudPlanar.cpp)
add_library(exp_help_gen ../utils/generatePointCloud.h ../utils/generatePointCloud.cpp)


# Example with zero translation
add_executable(example_basic_zero ${CMAKE_CURRENT_SOURCE_DIR}/example_basic_zero.cpp)
target_link_libraries(example_basic_zero 
        RotPriorH 
        exp_help_gen)

# Example with plane
add_executable(example_basic_plane ${CMAKE_CURRENT_SOURCE_DIR}/example_basic_plane.cpp)
target_link_libraries(example_basic_plane 
        RotPriorH 
        exp_help )
        
  
# Example with plane and certifier
add_executable(example_certifier_plane ${CMAKE_CURRENT_SOURCE_DIR}/example_certifier_plane.cpp)
target_link_libraries(example_certifier_plane 
        RotPriorH 
        exp_help 
        IterCertAlg)
        
      
