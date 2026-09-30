include(FetchContent)

# Upstream source lives only in the build tree (or a user-supplied directory).
# PETL uses headers from Eigen/libigl and only the MMA solver translation unit;
# do not configure upstream demos, viewers, tests or unrelated dependencies.
function(petl_source name url sha256 marker)
    string(TOUPPER "${name}" upper_name)
    set(PETL_${upper_name}_SOURCE_DIR "" CACHE PATH "Existing ${name} source tree (optional)")
    if(PETL_${upper_name}_SOURCE_DIR)
        set(source_dir "${PETL_${upper_name}_SOURCE_DIR}")
    elseif(PETL_FETCH_DEPENDENCIES)
        FetchContent_Declare(${name}
            URL "${url}"
            URL_HASH "SHA256=${sha256}"
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
            SOURCE_SUBDIR petl-no-upstream-build)
        FetchContent_MakeAvailable(${name})
        set(source_dir "${${name}_SOURCE_DIR}")
    else()
        message(FATAL_ERROR "Missing ${name}: set PETL_${upper_name}_SOURCE_DIR or enable PETL_FETCH_DEPENDENCIES")
    endif()
    if(NOT EXISTS "${source_dir}/${marker}")
        message(FATAL_ERROR "Invalid ${name} source directory: ${source_dir} (missing ${marker})")
    endif()
    set(${name}_SOURCE_DIR "${source_dir}" PARENT_SCOPE)
endfunction()

# Eigen 3.4.0 does not ship unsupported/Eigen/NNLS. Use an unmodified release
# containing that module instead of the original locally patched Eigen copy.
petl_source(eigen
    "https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz"
    315c881e19e17542a7d428c5aa37d113c89b9500d350c433797b730cd449c056
    unsupported/Eigen/NNLS)
add_library(petl_eigen INTERFACE)
target_include_directories(petl_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

# libigl v2.5.0; marching cubes retains the original copyleft implementation.
petl_source(libigl
    "https://codeload.github.com/libigl/libigl/tar.gz/fdaac01bcc52888994f7afd029dcc045dd408484"
    0c727f745a7096eebbef17d221d5cb3f3148d4c57ed0df3df698a8cf6b1dffca
    include/igl/copyleft/marching_cubes.h)
add_library(petl_libigl INTERFACE)
target_include_directories(petl_libigl SYSTEM INTERFACE "${libigl_SOURCE_DIR}/include")
target_link_libraries(petl_libigl INTERFACE petl_eigen)

# This revision is identical to the previously bundled MMASolver.h/.cpp.
petl_source(mma
    "https://codeload.github.com/jdumas/mma/tar.gz/aa51333f942220ac98e5957accb1b7e60590ec6f"
    953202fd53b16225555ea2196bb8928e44942fe50aa0cb0270554b425f471258
    src/mma/MMASolver.cpp)
add_library(petl_mma STATIC "${mma_SOURCE_DIR}/src/mma/MMASolver.cpp")
target_include_directories(petl_mma SYSTEM PUBLIC "${mma_SOURCE_DIR}/src/mma")
target_compile_features(petl_mma PUBLIC cxx_std_17)
if(MSVC)
    target_compile_options(petl_mma PRIVATE /utf-8)
endif()
