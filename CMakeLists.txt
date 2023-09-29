cmake_minimum_required(VERSION 3.12) #specify minimum version required 

project(CFD3D)

# specify the standard you want to use
set(CMAKE_CXX__STANDARD 14)

# set the directories for source code
set(SOURCE_DIR src)
# set the diretcories for headear files
set(INCLUDE_DIR include)

# Include directories for header files
include_directories(${INCLUDE_DIR})

# list source files
set(SOURCES
    src/header.cpp
    src/main.cpp
    src/solvePrediction.cpp
    src/solvePressure.cpp
    src/solveCorrection.cpp
    src/calculateError.cpp
    src/exportToVTK.cpp
)

# list header files
set(HEADERS
    include/header.h
)

# create executable from source files
add_executable(CFD3D ${SOURCES} ${HEADERS})

# link additional dependencies
# target_link_libraries(CFD3D library_name)