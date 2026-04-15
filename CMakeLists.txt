cmake_minimum_required(VERSION 3.18)
project(msnoise_tomo LANGUAGES C CXX)

# -------------------------
# Python
# -------------------------
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# -------------------------
# NumPy include dir (via Python)
# -------------------------
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "NumPy include dir: ${NUMPY_INCLUDE_DIR}")

# -------------------------
# mk_MatPaths
# -------------------------
Python_add_library(mk_MatPaths MODULE src/mk_MatPaths.c)
target_include_directories(mk_MatPaths PRIVATE ${NUMPY_INCLUDE_DIR})
set_target_properties(mk_MatPaths PROPERTIES PREFIX "")
target_link_libraries(mk_MatPaths PRIVATE Python::Module)

# -------------------------
# mkMatSmoothing
# -------------------------
Python_add_library(mkMatSmoothing MODULE src/mkMatSmoothing.c)
target_include_directories(mkMatSmoothing PRIVATE ${NUMPY_INCLUDE_DIR})
set_target_properties(mkMatSmoothing PROPERTIES PREFIX "")
target_link_libraries(mkMatSmoothing PRIVATE Python::Module)

# -------------------------
# vg_fta (C++)
# -------------------------
Python_add_library(vg_fta MODULE
    src/configparser.cpp
    src/fft_NR.cpp
    src/fta_param.cpp
    src/libfta.cpp
    src/readsac.cpp
    src/vg_fta.cpp
)
target_include_directories(vg_fta PRIVATE ${NUMPY_INCLUDE_DIR})
set_target_properties(vg_fta PROPERTIES PREFIX "")
target_link_libraries(vg_fta PRIVATE Python::Module)

# -------------------------
# Windows: export symbols via .def files
# -------------------------
if (WIN32)
    set_target_properties(mk_MatPaths PROPERTIES
        LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/mk_MatPaths.def"
    )
    set_target_properties(mkMatSmoothing PROPERTIES
        LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/mkMatSmoothing.def"
    )
    set_target_properties(vg_fta PROPERTIES
        LINK_FLAGS "/DEF:${CMAKE_SOURCE_DIR}/src/vg_fta.def"
    )
endif()

# -------------------------
# Install into Python package
# -------------------------
install(TARGETS
    mk_MatPaths
    mkMatSmoothing
    vg_fta
    LIBRARY DESTINATION msnoise_tomo/lib
    RUNTIME DESTINATION msnoise_tomo/lib
)