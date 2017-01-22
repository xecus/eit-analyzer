cmake_minimum_required(VERSION 2.6)
PROJECT(EIT)

IF( CMAKE_SIZEOF_VOID_P EQUAL 4 )
    MESSAGE( STATUS "32 bit architecture detected" )
ENDIF()

IF( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    MESSAGE( STATUS "64 bit architecture detected" )
ENDIF()

include_directories(~/eigen-eigen-da9b4e14c255)
include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(myapp_sources main.cpp)
ADD_EXECUTABLE(EIT ${myapp_sources})

target_link_libraries(EIT ncurses)
