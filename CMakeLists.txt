cmake_minimum_required(VERSION 3.6)
project(ctest_again)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES bcf2.cpp)
add_executable(ctest_again ${SOURCE_FILES})