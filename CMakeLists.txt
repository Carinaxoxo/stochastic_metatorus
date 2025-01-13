cmake_minimum_required(VERSION 3.25)
project(StochasticMetaTorus)

set(CMAKE_CXX_STANDARD 17)

add_executable(StochasticMetaTorus main.cpp metatorus.cpp metatorus.h bbst.cpp bbst.h random.h)
