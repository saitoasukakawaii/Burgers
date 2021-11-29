set(INC_DIR ~/eigen-3.4.0 ~/boost_1_77_0)
include_directories(${INC_DIR})

add_library(Polynomial STATIC GradJacobiP.cpp GradVandermonde1D.cpp JacobiP.cpp JacobiGQ.cpp JacobiGL.cpp Jacobi1D.cpp Vandermonde1D.cpp)

target_include_directories(Polynomial
        INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
        )
