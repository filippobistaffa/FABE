cmake_minimum_required(VERSION 3.9)
option(PROFILE "Enable Gperftools CPU profiler" OFF)
option(HEAP "Enable Gperftools heap profiler" OFF)

set(NAME fabe)
project(${NAME})
set(CMAKE_CXX_STANDARD 20)

set(PEDANTIC_COMPILE_FLAGS
    #-pedantic-errors
    -Wall
    -Wextra
    -pedantic
)

set(EXCLUDED_WARNINGS
    -Wno-sign-compare
    -Wno-unused-result
)

set(OPTIMIZATION_COMPILE_FLAGS
    -Ofast
    -march=native
    -funroll-loops
    -funsafe-loop-optimizations
    -falign-functions=16
    -falign-loops=16
    -fopenmp
)

set(CXX_SOURCES
    be.cpp
    conversion.cpp
    io.cpp
    main.cpp
    order.cpp
)

set(C_SOURCES
    libfa/fa.c
    libfa/hash.c
    libfa/ref.c
)

include_directories(libfa)

include(FetchContent)
FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

find_package(OpenMP REQUIRED)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(Gperftools)

add_executable(${NAME} ${C_SOURCES} ${CXX_SOURCES})
#set_source_files_properties(${CXX_SOURCES} PROPERTIES COMPILE_FLAGS -Wnrvo)
target_compile_options(${NAME} PRIVATE ${PEDANTIC_COMPILE_FLAGS} ${EXCLUDED_WARNINGS} ${OPTIMIZATION_COMPILE_FLAGS})
set_target_properties(${NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
target_link_libraries(${NAME} PRIVATE OpenMP::OpenMP_CXX fmt::fmt)

if(GPERFTOOLS_FOUND)
    LIST(APPEND LINK_LIBRARIES ${GPERFTOOLS_LIBRARIES})
    if(PROFILE)
        message(STATUS "Enabling Gperftools CPU profiler")
        LIST(APPEND LINK_LIBRARIES profiler)
        target_compile_definitions(${NAME} PRIVATE CPU_PROFILER)
        target_compile_options(${NAME} PRIVATE -g)
    endif()
    if(HEAP)
        message(STATUS "Enabling Gperftools heap profiler")
        target_compile_definitions(${NAME} PRIVATE HEAP_PROFILER)
    endif()
    target_link_libraries(${NAME} PRIVATE ${LINK_LIBRARIES})
endif()
