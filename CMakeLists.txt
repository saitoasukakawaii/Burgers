cmake_minimum_required(VERSION 3.16)
project(Burgers)


set(INC_DIR ./ ~/eigen-3.4.0 ~/boost_1_77_0)
set(CMAKE_CXX_STANDARD 17)


include_directories(${INC_DIR})

add_subdirectory(Polynomial)
list(APPEND EXTRA_LIBS Polynomial)



add_executable(Burgers main.cpp MeshGen1D.cpp GeometricFactors1D.cpp BurgersRHS1D.cpp)

target_link_libraries(Burgers PUBLIC ${EXTRA_LIBS})

target_include_directories(Burgers PUBLIC
        "${PROJECT_BINARY_DIR}"
        )


