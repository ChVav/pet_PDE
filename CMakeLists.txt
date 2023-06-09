cmake_minimum_required(VERSION 3.22)

# Specify the name of the project itself
project(PDE)

# Create target executables
add_executable(PDEsolver-dense)
add_executable(PDE-dense-bench) #testing benchmarking dense implementation, for now benchmarking outputted to the console
add_executable(PDEsolver-sparse-man src/common.cpp src/sparse-man.cpp)
add_executable(PDE-sparse-man-bench bench/PDE-sparse-man-bench.cpp)

# Make an additional directory for storing output PDEsolver
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

# Copy assets to output directory
add_custom_target(copy_assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets
)
add_dependencies(PDEsolver-dense copy_assets)
add_dependencies(PDE-dense-bench copy_assets)
add_dependencies(PDEsolver-sparse-man copy_assets)
add_dependencies(PDE-sparse-man-bench copy_assets)

# Add src and bench as a subdirectories, including contained (private) CMakeLists.txt files
add_subdirectory(src)
add_subdirectory(bench)

target_include_directories(PDEsolver-dense
    PRIVATE
        "src"
)

target_include_directories(PDE-dense-bench
	PRIVATE
		"bench"
)

target_include_directories(PDEsolver-sparse-man
	PRIVATE
		"src"
)

target_include_directories(PDE-sparse-man-bench
	PRIVATE
		"bench"
)

# build google benchmark (target: benchmark)
# do not build tests of benchmarking lib
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's test" FORCE)
add_subdirectory(vendor/google/benchmark)

# link targets to the benchmark library
target_link_libraries("${PROJECT_NAME}-dense-bench" benchmark)
target_link_libraries("${PROJECT_NAME}-sparse-man-bench" benchmark)

