#
#  This is a sample CMakeLists.txt suitable for direct use with a new PETSc application or
#  you can add fragments of material below to an existing applications CMakeLists.txt
#
#  Usage:  To build ex1 from ex1.c
#    rm -fr build
#    mkdir build
#    cd build
#    cmake ..
#    make
#    ./ex1
#
#  By default it gets both the compiler information and the library information from PETSc
#
cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})

set(ADOLC $ENV{ADOLC_DIR})

set(LAPACK $ENV{LAPACK_DIR})

set(HDF5 $ENV{HDF5_DIR}/CMakeFiles)

# PKG_CONFIG_PATH is set mannually in bashrc

# Remove the lines below if you do not wish to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)

#  tells CMake to build the application ex1 from the source file ex1.c
#  this must appear AFTER the compilers are set

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
#find_package(ADOLC REQUIRED)

project(test)
#add_executable(original_weno testweno.cpp func.cpp ../src/integral.cpp ../src/weno/flux.cpp ../src/mesh.c ../src/weno/weno_basis.cpp ../src/output.c ../src/input.cpp ../src/hdf5_io.c)
add_executable(multilevel_weno multiweno.cpp func.cpp ../src/integral.cpp ../src/mesh.c ../src/weno/weno_multilevel.cpp ../src/output.c ../src/input.cpp ../src/hdf5_io.c ../src/weno/flux_multilevel.cpp)

include_directories(../include ../src)

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_search_module(ADOLC REQUIRED IMPORTED_TARGET adolc)
pkg_search_module(LAPACKE REQUIRED IMPORTED_TARGET lapacke)
pkg_search_module(HDF5 REQUIRED IMPORTED_TARGET hdf5)
#target_link_libraries(original_weno PkgConfig::PETSC PkgConfig::ADOLC PkgConfig::LAPACKE PkgConfig::HDF5)
target_link_libraries(multilevel_weno PkgConfig::PETSC PkgConfig::ADOLC PkgConfig::LAPACKE PkgConfig::HDF5)
