cmake_minimum_required(VERSION 3.24)
project(PETL VERSION 1.0.0 LANGUAGES CXX)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Use an out-of-source build: cmake -S . -B build")
endif()

option(PETL_FETCH_DEPENDENCIES "Download missing dependencies into the build tree" ON)
option(PETL_ENABLE_OPENMP "Enable OpenMP parallel loops when available" ON)
include(cmake/Dependencies.cmake)

add_library(petl_core STATIC src/Truss.cpp src/config.cpp src/homo3D.cpp src/plsm.cpp)
target_include_directories(petl_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_compile_features(petl_core PUBLIC cxx_std_17)
target_link_libraries(petl_core PUBLIC petl_eigen petl_libigl PRIVATE petl_mma)
set_target_properties(petl_core PROPERTIES CXX_EXTENSIONS OFF)
if(MSVC)
    target_compile_options(petl_core PUBLIC /utf-8 /bigobj)
    target_compile_definitions(petl_core PUBLIC NOMINMAX _CRT_SECURE_NO_WARNINGS)
endif()

if(PETL_ENABLE_OPENMP)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(petl_core PUBLIC OpenMP::OpenMP_CXX)
        target_link_libraries(petl_mma PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(STATUS "OpenMP unavailable; using serial execution")
    endif()
endif()

add_executable(petl app/main.cpp)
target_link_libraries(petl PRIVATE petl_core)
set_target_properties(petl PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(CTest)
if(BUILD_TESTING)
    add_test(NAME cli_help COMMAND petl --help)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(NAME integration COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration.py"
            "$<TARGET_FILE:petl>" "${CMAKE_CURRENT_SOURCE_DIR}"
            "${CMAKE_CURRENT_BINARY_DIR}/test-output")
        set_tests_properties(integration PROPERTIES TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2")
    else()
        message(STATUS "Python 3 not found; integration tests will not be registered")
    endif()
endif()
