cmake_minimum_required(VERSION 3.20)

project(WoScMCF_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(raylib REQUIRED)

add_executable(WoScMCF
    src/main.cpp
    mantis/mantis.cpp
    mantis/Delaunay_psm.cpp
)

target_compile_options(WoScMCF PRIVATE -march=native)

set_source_files_properties(
    mantis/mantis.cpp mantis/Delaunay_psm.cpp
    PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")

target_include_directories(WoScMCF PRIVATE mantis)
target_link_libraries(WoScMCF PRIVATE raylib)