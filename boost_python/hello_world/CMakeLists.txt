cmake_minimum_required(VERSION 3.8)

project(hello_pylib)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
# Find necessary packages
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
 
find_package(Boost COMPONENTS python3 REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
 
# Define the wrapper library that wraps our library
add_library(hello_pylib SHARED hello.cpp)
target_link_libraries(hello_pylib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
set_target_properties(hello_pylib PROPERTIES PREFIX "")
