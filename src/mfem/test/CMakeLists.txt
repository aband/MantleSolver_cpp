cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})

#set(ADOLC $ENV{ADOLC_DIR})

set(LAPACK $ENV{LAPACK_DIR})

set(HDF5 $ENV{HDF5_CMAKE})

set(SRC $ENV{MantleSolver}/src)

set(CGNS $ENV{CGNS_DIR})

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

# Set Optimization flags
#set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")
set(CMAKE_CXX_FLAGS "-g")

project(test)

# Find packages using find_package 
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_search_module(LAPACKE REQUIRED IMPORTED_TARGET lapacke)
pkg_search_module(HDF5 REQUIRED IMPORTED_TARGET hdf5)
#pkg_search_module(CGNS REQUIRED IMPORTED_TARGET cgns)

# Build library
# Create a library from the files in src file
add_library(general_lib
        ${SRC}/general/integral.cpp 
        ${SRC}/general/mesh.c 
        ${SRC}/general/util.cpp
)

target_include_directories(general_lib
    PUBLIC ${SRC}/general
	 PUBLIC ${PETSC_DIR}/include
)

target_link_libraries(general_lib
    PkgConfig::PETSC
)

#Create a library from files in mfem file
add_library(mfem_lib
        ${SRC}/mfem/basis.cpp
        ${SRC}/mfem/Hdivmixed.cpp
        ${SRC}/mfem/brmixed.cpp
        ${SRC}/mfem/assemble.cpp
        ${SRC}/mfem/bndry.cpp
        ${SRC}/mfem/solve.cpp
        ${SRC}/mfem/error.cpp
		  ${SRC}/mfem/parallel/pbndry.cpp
		  ${SRC}/mfem/parallel/passemble.cpp
)

target_include_directories(mfem_lib
		  PUBLIC ${SRC}/mfem
		  PUBLIC ${SRC}/mfem/parallel
		  #		  PUBLIC /workspace/CGNS/src
		  #		  PUBLIC /workspace/CGNS/build/src
)

target_link_libraries(mfem_lib
    general_lib
)

add_library(io_lib
        ${SRC}/io/input.cpp
        ${SRC}/io/output.c
)

target_include_directories(io_lib
		  PUBLIC ${SRC}/io
		  PUBLIC ${PETSC_DIR}/include
)

target_link_libraries(io_lib
    general_lib
    PkgConfig::PETSC
    PkgConfig::HDF5
	 #    ${CGNS}/libcgns.a
#	 PkgConfig::CGNS
)

# Add to excutable
add_executable(testCoupled
        testCoupled.cpp 
		  myfunc.cpp
		  #myfunc_mixedBndry.cpp
)

add_executable(testSeparate
        testSeparate.cpp
        myfunc4.cpp
)

add_executable(testBasis
        testBasis.cpp
)

add_executable(testParallel
        testParallel.cpp
		  myfunc4.cpp
)

# Link static libraries generated previously
target_link_libraries(testCoupled 
    mfem_lib
    io_lib
    PkgConfig::LAPACKE 
)

target_link_libraries(testSeparate
    mfem_lib
    io_lib
    PkgConfig::LAPACKE 
)

target_link_libraries(testBasis
    mfem_lib
    io_lib
    PkgConfig::LAPACKE 
)

target_link_libraries(testParallel
    mfem_lib
    io_lib
    PkgConfig::LAPACKE 
)
