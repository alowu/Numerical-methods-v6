cmake_minimum_required(VERSION 3.17)
project(4islaki)

set(CMAKE_CXX_STANDARD 17)
file(GLOB run "Lab_2/main.cpp")

add_executable(4islaki ${run})