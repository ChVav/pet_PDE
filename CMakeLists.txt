cmake_minimum_required(VERSION 3.22)

# Specify the name of the project itself
project(PDEsolver)

# Create a target, being an executable
add_executable(PDEsolver)

# Copy assets to output directory
add_custom_target(copy_assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/assets ${CMAKE_CURRENT_BINARY_DIR}/assets
)
add_dependencies(PDEsolver copy_assets)

# Add src as a subdirectory, meaning, include the file src/CMakeLists.txt
add_subdirectory(src)

target_include_directories(PDEsolver
    PRIVATE
        "src"
)
