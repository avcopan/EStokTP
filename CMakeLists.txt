project(EStokTP C Fortran)
cmake_minimum_required(VERSION 2.8)
find_package(Threads)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -fno-leading-underscore")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

add_library(estoktp
    ${PROJECT_SOURCE_DIR}/src/estoktp/g09fopt.f
    ${PROJECT_SOURCE_DIR}/src/estoktp/LineRead.f
    ${PROJECT_SOURCE_DIR}/src/estoktp/croutines.c
    ${PROJECT_SOURCE_DIR}/src/estoktp/mt.c
    ${PROJECT_SOURCE_DIR}/src/estoktp/molprofopt.f
    ${PROJECT_SOURCE_DIR}/src/estoktp/zmat_routines.f)

add_library(projRot
    ${PROJECT_SOURCE_DIR}/src/projRot/eigen2.c
    ${PROJECT_SOURCE_DIR}/src/projRot/eigen.c
    ${PROJECT_SOURCE_DIR}/src/projRot/nrutil.c
    ${PROJECT_SOURCE_DIR}/src/projRot/spline.c
    ${PROJECT_SOURCE_DIR}/src/projRot/invert_matrix3.c)

add_executable(estoktp.x
    ${PROJECT_SOURCE_DIR}/src/estoktp/estoktp.f)

add_executable(RPHt.exe
    ${PROJECT_SOURCE_DIR}/src/projRot/main.c)

target_link_libraries(estoktp.x estoktp)

target_link_libraries(RPHt.exe projRot m)

install(TARGETS estoktp.x DESTINATION bin)
install(TARGETS RPHt.exe DESTINATION bin)
