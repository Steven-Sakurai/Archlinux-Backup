project(heston_svj_mle)
set(CMAKE_CXX_STANDARD 11)
find_package(Boost COMPONENTS system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
include_directories(include)
add_executable(a.out test_initialMoments.cpp)
target_link_libraries(a.out ${Boost_LIBRARIES})

