cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(lattgillespie)

include_directories(/usr/local/include)

AUX_SOURCE_DIRECTORY(src SRC_FILES) 

file (GLOB HEADER_FILES "include/lattgillespie_bits/*.hpp")

add_library(lattgillespie SHARED 
	src/lattice_gillespie.cpp
	src/lattice.cpp
	src/lattice.hpp
	src/reactions.cpp
	src/reactions.hpp
	src/species.cpp
	src/species.hpp
)

target_include_directories(lattgillespie PRIVATE include/lattgillespie_bits)

install(TARGETS lattgillespie DESTINATION lib)

install(FILES 
	include/lattgillespie 
	DESTINATION include)

install(DIRECTORY 
	include/lattgillespie_bits 
	DESTINATION include)