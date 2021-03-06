cmake is designed to be "Write once, Run everywhere" 

## USAGE
For small projects,  

```bash
cmake .
make
./output
```

For bigger projects,   

```bash
mkdir build && cd build
cmake ..
make
./output
```


```
cmake_minimum_required(VERSION 3.5) # not necessary really
project(projectName)
set(CMAKE_CXX_STANDARD 14) # use new standard
add_executable(output main.cpp fun1.cpp fun2.cpp ...)
```

```
# add all src under (current) directory and save in a var 'DIR_SRC'
aux_source_directory(. DIR_SRC)
add_executable(output ${DIR_SRC})
```

```
# headers in 'include', sources in 'src'
include_directories(include)
# only use some files in 'src'
set(SOURCES src/main.cpp src/a.cpp ...)
# or use wildcard additions
file(GLOB SOURCES "src/*.cpp")
add_executable(output ${SOURCES})
```
  
Armadillo
```
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(output ${ARMADILLO_LIBRARIES})
```  

Boost
```
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
target_link_libraries(output ${Boost_LIBRARIES})
```


Example

```
project(TestCmake)
set(CMAKE_CXX_STANDARD 11)
aux_source_directory(./src DIR_SRC)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
include_directories(include)
add_executable(output ${DIR_SRC} main.cc)
target_link_libraries(output ${ARMADILLO_LIBRARIES})
target_link_libraries(output ${Boost_LIBRARIES})
```
