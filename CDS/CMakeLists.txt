cmake_minimum_required(VERSION 2.4.6)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries

project(cds)

find_package(catkin REQUIRED COMPONENTS roscpp mathlib gmr)

catkin_package(
 INCLUDE_DIRS include
 LIBRARIES ${PROJECT_NAME}
 CATKIN_DEPENDS mathlib gmr
)

include_directories(include ${catkin_INCLUDE_DIRS})

FILE(GLOB SRC_FILES src/*.cpp)
#common commands for building c++ executables and libraries
add_library(${PROJECT_NAME} ${SRC_FILES})

