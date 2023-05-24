cmake_minimum_required(VERSION 3.13)
project(shared LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(-g -Wall -Wextra -O3)
    message(STATUS "Added compiler options: -g -Wall -Wextra -O3 -pthread")
endif()

find_package(Boost REQUIRED)
include_directories(${BOOST_INCLUDE_DIR} include/potential 
include/vectors include/consts include/df
include/utility include/multithreading)

file(GLOB src_files src/*/*.cpp include/*/*.hpp)
add_library(mylib SHARED ${src_files})
set_target_properties(mylib PROPERTIES SUFFIX ".so")

target_link_libraries(mylib ${BOOST_LIBRARIES})

target_include_directories(mylib INTERFACE include/potential 
include/vectors include/df include/consts 
include/utility include/multithreading)



include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    unit_tests
    ${src_files})

target_link_libraries(
    unit_tests
    mylib
    GTest::gtest_main
    fftw3)

include(GoogleTest)
gtest_discover_tests(unit_tests)