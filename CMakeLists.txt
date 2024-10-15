cmake_minimum_required(VERSION 3.26)
project(parametric_maize_leaf)

set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED COMPONENTS Core)
find_package(GLM REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMesh CONFIG REQUIRED)
include(${CGAL_USE_FILE})
include_directories(${Eigen3_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})
add_definitions     ( ${EIGEN3_DEFINITIONS} )


add_executable(parametric_maize_leaf main.cpp
        leaf_generation.h
)

target_compile_definitions(parametric_maize_leaf PRIVATE -D_USE_MATH_DEFINES)

target_link_libraries(parametric_maize_leaf
        CGAL::CGAL
        GSL::gsl
        GSL::gslcblas
        glm::glm
        ${Eigen3_LIBRARIES}
        ${OPENMESH_LIBRARIES}
)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
