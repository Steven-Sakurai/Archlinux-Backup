 #!/usr/local/bin/python2 
 # -*- coding: utf-8 -*-

import sys
filename = "CMakeLists.txt"
with open(filename, 'w') as f:
    f.write(r"""#cmake_minimum_required(VERSION 3.8)
project (PROJECT_NAME)
# new C++ standard
set (CMAKE_CXX_STANDARD 14) 
# Set the output folder where your program will be created
#set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
# set src files
file (GLOB SOURCES "./*.cpp")
# all headers put in 'include'
include_directories(include)
add_executable(output ${SOURCES})
""")
