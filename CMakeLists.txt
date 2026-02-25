cmake_minimum_required(VERSION 3.16)
project(HelloWorld CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output binaries to a predictable location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(include)

add_executable(hello
    src/ball.cpp
    src/main.cpp
    src/landscape.cpp
    src/outputcode.cpp
    src/set.cpp
    src/simple.hpp
    src/type.cpp
)

target_link_libraries(hello Eigen3::Eigen)
