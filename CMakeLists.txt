project(dyn2b)
cmake_minimum_required(VERSION 2.8)

option(BUILD_TEST "Build unit tests" Off)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include_directories(
  include
)

add_subdirectory(src)


if(BUILD_TEST)
  enable_testing()

  find_package(Check)
  add_subdirectory(test)
endif()