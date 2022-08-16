#
# Demonstrates and tests how other projects can use sotb
#
cmake_minimum_required(VERSION 3.16)

project(TestBuild LANGUAGES C Fortran)

find_package(sotb REQUIRED)

add_executable(test_build main.c)
target_link_libraries(test_build sotb::sotb)
