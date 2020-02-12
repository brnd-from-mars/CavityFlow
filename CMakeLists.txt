cmake_minimum_required(VERSION 3.12)
project(NavierStokesEquation)

set(CMAKE_CXX_STANDARD 17)

find_package(SFML 2.5 COMPONENTS graphics REQUIRED)

add_executable(NavierStokesEquation main.cpp HSL.cpp HSL.hpp)
target_link_libraries(NavierStokesEquation sfml-graphics)
