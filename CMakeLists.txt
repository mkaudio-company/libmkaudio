cmake_minimum_required(VERSION 3.22)
project(mkaudio VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(mkaudio
    iirfilter.cpp
    saturation.cpp
    simulation.cpp
)

target_include_directories(mkaudio PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)
