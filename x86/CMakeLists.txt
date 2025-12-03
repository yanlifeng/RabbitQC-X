cmake_minimum_required(VERSION 3.12)

# Option to enable MPI support
option(USE_MPI_IO "Enable MPI I/O support" OFF)

# Option to enable verbose debug output
option(VERBOSE "Enable verbose debug output" OFF)

# If MPI is enabled, use MPI compiler wrappers
if(USE_MPI_IO)
    message(STATUS "MPI I/O support enabled - using MPI compiler wrappers")
    find_package(MPI REQUIRED)
    
    # Use MPI compiler wrappers
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    
    message(STATUS "Using MPI C compiler: ${MPI_C_COMPILER}")
    message(STATUS "Using MPI C++ compiler: ${MPI_CXX_COMPILER}")
else()
    message(STATUS "MPI I/O support disabled - using standard compilers")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set C standard
set(CMAKE_C_STANDARD 11)

# Check GCC version and CPU flags for instruction set selection
if(NOT DEFINED BIOCONDA)
    # Detect CPU features
    execute_process(COMMAND lscpu OUTPUT_VARIABLE LSCPU_OUTPUT)
    
    string(FIND "${LSCPU_OUTPUT}" "sse2" SSE2_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "sse4_1" SSE4_1_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "avx" AVX_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "avx2" AVX2_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "avx512f" AVX512F_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "avx512bw" AVX512BW_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "avx512vl" AVX512VL_FOUND)
    string(FIND "${LSCPU_OUTPUT}" "popcnt" POPCNT_FOUND)
    
    # Determine instruction set
    if(NOT SSE2_FOUND EQUAL -1 AND NOT SSE4_1_FOUND EQUAL -1 AND 
       NOT AVX_FOUND EQUAL -1 AND NOT AVX2_FOUND EQUAL -1 AND 
       NOT AVX512F_FOUND EQUAL -1 AND NOT AVX512BW_FOUND EQUAL -1 AND 
       NOT AVX512VL_FOUND EQUAL -1 AND NOT POPCNT_FOUND EQUAL -1)
        message(STATUS "Using AVX512 instruction set")
        set(INSTRUCT_SET "-DVec512")
        set(ARCH_FLAGS "-march=native")
    elseif(NOT SSE2_FOUND EQUAL -1 AND NOT SSE4_1_FOUND EQUAL -1 AND 
           NOT AVX_FOUND EQUAL -1 AND NOT AVX2_FOUND EQUAL -1 AND 
           NOT POPCNT_FOUND EQUAL -1)
        message(STATUS "Using AVX2 instruction set")
        set(INSTRUCT_SET "-DVec256")
        set(ARCH_FLAGS "-march=native")
    else()
        message(STATUS "Using compiler automatic vectorization")
        set(INSTRUCT_SET "")
        set(ARCH_FLAGS "-march=native")
    endif()
else()
    message(STATUS "BIOCONDA build: Using SSE4.2 instruction set")
    set(INSTRUCT_SET "")
    set(ARCH_FLAGS "-msse4.2")
endif()

# Thrust include directory
set(THRUST_INCDIR ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/thrust-1.17.0)

# Check GCC version for constexpr options
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Detected GCC version: ${GCC_VERSION}")
    
    # -fconstexpr-ops-limit is only available in GCC 9+
    if(GCC_VERSION VERSION_GREATER_EQUAL "9.0")
        set(CONSTEXPR_FLAG "-fconstexpr-ops-limit=99000100")
        message(STATUS "Using -fconstexpr-ops-limit (GCC 9+)")
    else()
        # For older GCC, use -fconstexpr-depth instead (available since GCC 4.6)
        set(CONSTEXPR_FLAG "-fconstexpr-depth=512")
        message(STATUS "Using -fconstexpr-depth for GCC < 9")
    endif()
else()
    set(CONSTEXPR_FLAG "")
endif()

# Source files from platform-specific directory
# Collect all source files (excluding main.cpp from platform src)
file(GLOB CXX_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER CXX_SOURCES EXCLUDE REGEX ".*main\\.cpp$")
file(GLOB C_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

# Create executable with unified main.cpp from root src directory
add_executable(rabbitqc-x ${CMAKE_SOURCE_DIR}/src/main.cpp ${CXX_SOURCES} ${C_SOURCES})

# Set target name
set_target_properties(rabbitqc-x PROPERTIES OUTPUT_NAME "rabbitqc-x")

# Include directories
target_include_directories(rabbitqc-x PRIVATE
    ${CMAKE_SOURCE_DIR}  # Add project root for platform-specific includes
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/huffman
    ${CMAKE_CURRENT_SOURCE_DIR}/include/chunkdecoding
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gzip
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${THRUST_INCDIR}
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/common/include
    ${CMAKE_SOURCE_DIR}/common/lib
)

# Platform macro for conditional compilation
target_compile_definitions(rabbitqc-x PRIVATE PLATFORM_X86)

# Add MPI definitions if enabled
if(USE_MPI_IO)
    target_compile_definitions(rabbitqc-x PRIVATE USE_MPI_IO)
    message(STATUS "Added USE_MPI_IO compile definition")
endif()

# Add Verbose definition if enabled
if(VERBOSE)
    target_compile_definitions(rabbitqc-x PRIVATE Verbose)
    message(STATUS "Verbose debug output enabled")
endif()

# Common compiler flags for both C and C++
target_compile_options(rabbitqc-x PRIVATE
    ${INSTRUCT_SET}
    ${ARCH_FLAGS}
    -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP
    -DNDUBUG
    -O3
    -w
    -fopenmp
    -Wall
    -Wextra
)

# C++ specific flags (only applied to C++ files)
target_compile_options(rabbitqc-x PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${CONSTEXPR_FLAG}>
    $<$<COMPILE_LANGUAGE:CXX>:-Wno-terminate>
    $<$<COMPILE_LANGUAGE:CXX>:-Wno-class-memaccess>
)

# C-specific flags (only applied to C files)
target_compile_options(rabbitqc-x PRIVATE
    $<$<COMPILE_LANGUAGE:C>:-Wno-unknown-pragmas>
    $<$<COMPILE_LANGUAGE:C>:-Wcast-qual>
)

# Find required libraries
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

# Link libraries
target_link_libraries(rabbitqc-x PRIVATE
    Threads::Threads
    ZLIB::ZLIB
    OpenMP::OpenMP_CXX
    rt
    stdc++fs
    dl
)

# Add MPI libraries if enabled
if(USE_MPI_IO)
    target_link_libraries(rabbitqc-x PRIVATE MPI::MPI_CXX)
    message(STATUS "Linked MPI libraries")
endif()

# Installation
install(TARGETS rabbitqc-x DESTINATION bin)

