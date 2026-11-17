cmake_minimum_required(VERSION 3.12.0)
project(lab02 LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_library(Poisson2D STATIC src/Poisson2D.cpp)
deal_ii_setup_target(Poisson2D)

add_executable(exercise-01 src/exercise-01.cpp)
target_link_libraries(exercise-01 Poisson2D)
deal_ii_setup_target(exercise-01)

add_library(DiffusionReaction STATIC src/DiffusionReaction.cpp)
deal_ii_setup_target(DiffusionReaction)

add_executable(exercise-02 src/exercise-02.cpp)
target_link_libraries(exercise-02 DiffusionReaction)
deal_ii_setup_target(exercise-02)