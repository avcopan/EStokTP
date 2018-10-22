project(EStokTP C Fortran)
cmake_minimum_required(VERSION 2.8)
find_package(BLAS REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-leading-underscore")

if(BLAS_FOUND)
    message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
endif()

add_library(estoktp
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/g09fopt.f
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/LineRead.f
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/croutines.c
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/mt.c
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/dihedral.c
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/molprofopt.f
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/zmat_routines.f)

add_executable(estoktp.x
    ${PROJECT_SOURCE_DIR}/sources/estoktp_1.8_distrib/estoktp.f)

target_link_libraries(estoktp.x estoktp blas)

install(TARGETS estoktp.x DESTINATION bin)
