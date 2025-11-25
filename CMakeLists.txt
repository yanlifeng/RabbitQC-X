cmake_minimum_required(VERSION 3.12)
project(RabbitQC-X VERSION 1.0.0 LANGUAGES C CXX)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Platform selection: x86 or sunway
if(NOT DEFINED PLATFORM)
    set(PLATFORM "x86" CACHE STRING "Target platform: x86 or sunway")
endif()

message(STATUS "Building for platform: ${PLATFORM}")

# Validate platform selection
if(NOT ${PLATFORM} STREQUAL "x86" AND NOT ${PLATFORM} STREQUAL "sunway")
    message(FATAL_ERROR "Invalid PLATFORM: ${PLATFORM}. Must be 'x86' or 'sunway'")
endif()

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Common include directory
include_directories(${CMAKE_SOURCE_DIR}/common)

# Build the selected platform
if(${PLATFORM} STREQUAL "x86")
    message(STATUS "Configuring RabbitQCPlus for x86 platform")
    add_subdirectory(x86)
elseif(${PLATFORM} STREQUAL "sunway")
    message(STATUS "Configuring SWQC for Sunway platform")
    add_subdirectory(sunway)
endif()

