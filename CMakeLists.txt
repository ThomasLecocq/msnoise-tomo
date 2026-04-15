cmake_minimum_required(VERSION 3.18)
project(msnoise_tomo LANGUAGES C CXX)

# Python + NumPy
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# Get NumPy include dir via Python
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "NumPy include dir: ${NUMPY_INCLUDE_DIR}")

include_directories(${NUMPY_INCLUDE_DIR})

# Output dir = package/lib
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/msnoise_tomo/lib)

# -------------------------
# mk_MatPaths
# -------------------------
add_library(mk_MatPaths MODULE src/mk_MatPaths.c)
target_include_directories(mk_MatPaths PRIVATE ${NumPy_INCLUDE_DIRS})
set_target_properties(mk_MatPaths PROPERTIES PREFIX "")

# -------------------------
# mkMatSmoothing
# -------------------------
add_library(mkMatSmoothing MODULE src/mkMatSmoothing.c)
target_include_directories(mkMatSmoothing PRIVATE ${NumPy_INCLUDE_DIRS})
set_target_properties(mkMatSmoothing PROPERTIES PREFIX "")

# -------------------------
# vg_fta (C++)
# -------------------------
add_library(vg_fta MODULE
    src/configparser.cpp
    src/fft_NR.cpp
    src/fta_param.cpp
    src/libfta.cpp
    src/readsac.cpp
    src/vg_fta.cpp
)

target_include_directories(vg_fta PRIVATE ${NumPy_INCLUDE_DIRS})
set_target_properties(vg_fta PROPERTIES PREFIX "")

# -------------------------
# Windows (.def exports)
# -------------------------
if (WIN32)
    set_target_properties(mk_MatPaths PROPERTIES LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/mk_MatPaths.def")
    set_target_properties(mkMatSmoothing PROPERTIES LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/mkMatSmoothing.def")
    set_target_properties(vg_fta PROPERTIES LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/vg_fta.def")
endif()