cmake_minimum_required(VERSION 3.17)
project(4islaki)

set(CMAKE_CXX_STANDARD 17)
file(GLOB run "Lab_4/main.cpp")
file(GLOB test "functions/*")

add_executable(4islaki ${run} ${test})