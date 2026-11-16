cmake_minimum_required(VERSION 3.16)
project(magic LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(MAGIC_HI_PRECISION "Use high precision" ON)
option(MAGIC_PYTHON "Build the Python module (set by pip/uv when installing specmagic)" OFF)

# Find omp
find_package(OpenMP REQUIRED)

# NetCDF (Cross-platform using CMake config instead of nc-config)
find_package(netCDF REQUIRED)

# Source code: everything but the CLI's main goes into a library that both the CLI and
# the Python module use, so that they run the same code
file(GLOB_RECURSE SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/*.cpp
    ${CMAKE_SOURCE_DIR}/headers/*.hpp
)
list(REMOVE_ITEM SRC_FILES ${CMAKE_SOURCE_DIR}/src/main.cpp)

add_library(magic_core STATIC ${SRC_FILES})
set_target_properties(magic_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Compiler flags, precision and libraries of everything built here
function(magic_settings target visibility)

    # Should work for both windows and linux, we hope...
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${target} ${visibility}
            -Wall -g -O3
            -ffast-math
            # Avoid -march=native in shared repos unless explicitly requested
        )
    elseif(MSVC)
        target_compile_options(${target} ${visibility}
            /W4 /O2 /fp:fast
        )
    endif()

    if(MAGIC_HI_PRECISION)
        target_compile_definitions(${target} ${visibility} MAGIC_PRECISION_HI)
    else()
        target_compile_definitions(${target} ${visibility} MAGIC_PRECISION_LO)
    endif()

    # Now do the linking
    target_link_libraries(${target} ${visibility}
        OpenMP::OpenMP_CXX
        netCDF::netcdf
    )

    # We only add the math library on unix
    # For windows it should be there automatically
    if(UNIX)
        target_link_libraries(${target} ${visibility} m)
    endif()
endfunction()

magic_settings(magic_core PUBLIC)

# The CLI
add_executable(magic ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(magic PRIVATE magic_core)

# The Python module, specmagic._magic
if(MAGIC_PYTHON)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_magic ${CMAKE_SOURCE_DIR}/python/bindings.cpp)
    target_link_libraries(_magic PRIVATE magic_core)
    install(TARGETS _magic DESTINATION specmagic)
endif()
