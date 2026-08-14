cmake_minimum_required(VERSION 3.16)
project(magic LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find omp
find_package(OpenMP REQUIRED)

# NetCDF (Cross-platform using CMake config instead of nc-config)
find_package(netCDF REQUIRED)

# Source code
file(GLOB_RECURSE SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
    ${CMAKE_SOURCE_DIR}/headers/*.hpp
)

add_executable(magic ${SRC_FILES})

# Should work for both windows and linux, we hope...
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(magic PRIVATE 
        -Wall -g -O3 
        -ffast-math
        # Avoid -march=native in shared repos unless explicitly requested
    )
elseif(MSVC)
    target_compile_options(magic PRIVATE 
        /W4 /O2 /fp:fast
    )
endif()

# Options
option(MAGIC_HI_PRECISION "Use high precision" ON)
if(MAGIC_HI_PRECISION)
    target_compile_definitions(magic PRIVATE MAGIC_PRECISION_HI)
else()
    target_compile_definitions(magic PRIVATE MAGIC_PRECISION_LO)
endif()

# Now do the linking
target_link_libraries(magic PRIVATE
    OpenMP::OpenMP_CXX
    netCDF::netcdf
)

# We only add the math library on unix
# For windows it should be there automatically
if(UNIX)
    target_link_libraries(magic PRIVATE m)
endif()
