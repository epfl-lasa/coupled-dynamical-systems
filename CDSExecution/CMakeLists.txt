cmake_minimum_required(VERSION 2.4.6)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries

project(cds_execution)

find_package(catkin REQUIRED COMPONENTS roscpp mathlib gmr cds)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS mathlib gmr cds
)

FILE(GLOB SRC_FILES src/*.cpp)

include_directories(include ${catkin_INCLUDE_DIRS})

#common commands for building c++ executables and libraries
add_library(${PROJECT_NAME} ${SRC_FILES})
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

#SET(V_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)
#SET(V_LIB_DIRS ${PROJECT_SOURCE_DIR}/lib)
#SET(V_LIBS ${PROJECT_NAME})
#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.in ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake @ONLY)
