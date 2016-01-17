cmake_minimum_required(VERSION 3.3)
project(assignment2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(SOURCE_FILES Chromosome.h DataReader.h Node.h GA2.cpp Chromosome.cpp DataReader.cpp Node.cpp )
add_executable(assignment2 ${SOURCE_FILES})