cmake_minimum_required(VERSION 3.12)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Use Sunway MPI compilers
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

# Source files from src directory (excluding main.cpp)
file(GLOB SRC_CXX_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER SRC_CXX_SOURCES EXCLUDE REGEX ".*main\\.cpp$")
file(GLOB SRC_C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

# Source files from slave directory
file(GLOB SLAVE_CXX_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/slave/*.cpp")
file(GLOB SLAVE_C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/slave/lib/*.c")

# Create object libraries for different compilation modes
# Host objects (main code)
add_library(host_objects OBJECT ${SRC_CXX_SOURCES})
target_compile_options(host_objects PRIVATE
    -mhost
    -DVerbose
    -std=c++11
    -g
    -O3
    -w
)
target_include_directories(host_objects PRIVATE
    ${CMAKE_SOURCE_DIR}  # Add project root for platform-specific includes
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/common/include
    ${CMAKE_SOURCE_DIR}/common/lib
)
target_compile_definitions(host_objects PRIVATE PLATFORM_SUNWAY)

# Host C objects
add_library(host_c_objects OBJECT ${SRC_C_SOURCES})
target_compile_options(host_c_objects PRIVATE
    -mhost
    -g
    -O3
    -w
)

# Slave objects (accelerator code)
add_library(slave_objects OBJECT ${SLAVE_CXX_SOURCES})
target_compile_options(slave_objects PRIVATE
    -mslave
    -msimd
    -DVerbose
    -std=c++11
    -g
    -O3
    -w
)
target_include_directories(slave_objects PRIVATE
    ${CMAKE_SOURCE_DIR}  # Add project root for platform-specific includes
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/common/include
    ${CMAKE_SOURCE_DIR}/common/lib
)
target_compile_definitions(slave_objects PRIVATE PLATFORM_SUNWAY)

# Slave C objects
add_library(slave_c_objects OBJECT ${SLAVE_C_SOURCES})
target_compile_options(slave_c_objects PRIVATE
    -mslave
    -msimd
    -g
    -O3
    -w
)

# Create executable with hybrid mode (using unified main.cpp from root src)
add_executable(rabbitqc-x
    ${CMAKE_SOURCE_DIR}/src/main.cpp
    $<TARGET_OBJECTS:host_objects>
    $<TARGET_OBJECTS:host_c_objects>
    $<TARGET_OBJECTS:slave_objects>
    $<TARGET_OBJECTS:slave_c_objects>
)

# Set target name and platform macro
set_target_properties(rabbitqc-x PROPERTIES OUTPUT_NAME "rabbitqc-x")
target_compile_definitions(rabbitqc-x PRIVATE PLATFORM_SUNWAY)
target_include_directories(rabbitqc-x PRIVATE
    ${CMAKE_SOURCE_DIR}  # Add project root for src/Globals.h, src/cmdinfo.h, etc.
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/common/include
    ${CMAKE_SOURCE_DIR}/common/lib
)

# Link flags for hybrid mode
target_link_options(rabbitqc-x PRIVATE
    -mhybrid
    -static
)

# Find static zlib library for Sunway platform
# Sunway requires static linking due to -static flag
set(ZLIB_USE_STATIC_LIBS ON)
find_library(ZLIB_LIBRARY_STATIC 
    NAMES libz.a z
    PATHS /usr/sw/lib /usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/lib
    NO_DEFAULT_PATH
)

if(ZLIB_LIBRARY_STATIC)
    message(STATUS "Found static zlib: ${ZLIB_LIBRARY_STATIC}")
else()
    # Fallback to find_package
    find_package(ZLIB REQUIRED)
    set(ZLIB_LIBRARY_STATIC ${ZLIB_LIBRARIES})
endif()

# Link libraries
target_link_libraries(rabbitqc-x PRIVATE
    ${ZLIB_LIBRARY_STATIC}
    pthread
    rt
)

# Installation
install(TARGETS rabbitqc-x DESTINATION bin)

