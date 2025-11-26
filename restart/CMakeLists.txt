cmake_minimum_required(VERSION 3.31)
project(restart)

set(CMAKE_CXX_STANDARD 20)

add_executable(restart main.cpp
        metatorus.cpp
        metatorus.h
        routing.cpp
        experiments.cpp
        experiments.h)

#cmake_minimum_required(VERSION 3.20)
#project(restart)
#
#set(CMAKE_CXX_STANDARD 17)
#
## 1. Define the library (The Logic)
## We combine metatorus.cpp and routing.cpp into one library
#add_library(core_lib
#        metatorus.cpp
#        routing.cpp
#)
#
## 2. Define the executable (The Runner)
#add_executable(MetaTorus_Run main.cpp)
#
## 3. Link them together
#target_link_libraries(MetaTorus_Run core_lib)