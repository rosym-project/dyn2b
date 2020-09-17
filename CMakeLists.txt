project(dyn2b)
cmake_minimum_required(VERSION 2.8)

option(BUILD_TEST "Build unit tests" Off)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(PackageRegistry)

include_directories(
  include
)

add_subdirectory(src)


if(BUILD_TEST)
  enable_testing()

  find_package(Check)
  add_subdirectory(test)
endif()


install(
  DIRECTORY include/
  DESTINATION include
  FILES_MATCHING PATTERN "*.h"
)

register_package(
  INCLUDE_DIRS
    include
  LIBRARIES
    dyn2b
    dyn2b_example
)
