cmake_minimum_required(VERSION 3.12)

# Nome do projeto
project(HippopotamusOptimization VERSION 1.0)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCES
    main.cpp
    HippopotamusOptimization.h
    TestFunctions.h
)

add_executable(HippoOpt ${SOURCES})

if(MSVC)
    target_compile_options(HippoOpt PRIVATE /W3)
else()
    target_compile_options(HippoOpt PRIVATE -Wall -Wextra)
endif()