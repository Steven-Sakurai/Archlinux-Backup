cmake is designed to be "Write once, Run everywhere" 

## USEAGE
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
cmake_minimum_required (VERSION 3.5) # not necessary really
project (projectName)
set (CMAKE_CXX_STANDARD 14) # use new standard
add_executable (output main.cpp fun1.cpp fun2.cpp ...)
```

```
# add all src under (current) directory and save in a var 'DIR_SRC'
aux_source_directory (. DIR_SRC)
add_executable (output ${DIR_SRC})
```

```
# headers in 'include', sources in 'src'
include_directories(include)
# only use some files in 'src'
set (SOURCES src/main.cpp src/a.cpp ...)
# or use wildcard additions
file (GLOB SOURCES "src/*.cpp")
add_executable (output ${SOURCES})
```
