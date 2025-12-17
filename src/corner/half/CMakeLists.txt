cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(PETSC $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})

#set(ADOLC $ENV{ADOLC_DIR})

#set(LAPACK $ENV{LAPACK_DIR})

set(HDF5 $ENV{HDF5_DIR})

set(SRC $ENV{MantleSolver}/src)

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
set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra")
#set(CMAKE_CXX_FLAGS "-g")

project(test)

# Find packages using find_package 
find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
pkg_search_module(LAPACKE REQUIRED IMPORTED_TARGET lapacke)
pkg_search_module(HDF5 REQUIRED IMPORTED_TARGET hdf5)
#pkg_search_module(CGNS REQUIRED IMPORTED_TARGET cgns)

# Build library
# Create a library from the files in general folder
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

# Create library from io folder 
add_library(io_lib
        ${SRC}/io/input.cpp
        ${SRC}/io/output.c
#		  ${SRC}/io/cgns_io.c
)

target_include_directories(io_lib
		  PUBLIC ${SRC}/io
		  PUBLIC ${PETSC_DIR}/include
)

target_link_libraries(io_lib
    general_lib
    PkgConfig::PETSC
    PkgConfig::HDF5
#	 PkgConfig::CGNS
)

#Create a library from phase folder
add_library(phase_lib
		  ${SRC}/phase/eutectic_rescaled.cpp
)

target_include_directories(phase_lib
		  PUBLIC ${SRC}/phase
		  PUBLIC ${SRC}/include
		  )

#Create a library from files in mfem file
add_library(mfem_lib
        ${SRC}/mfem/basis.cpp
        ${SRC}/mfem/Hdivmixed.cpp
        ${SRC}/mfem/brmixed.cpp
        #${SRC}/mfem/assemble.cpp
        #${SRC}/mfem/bndry.cpp
        #${SRC}/mfem/solve.cpp
        #${SRC}/mfem/error.cpp
		  #${SRC}/mfem/shape.tpp
		  #${SRC}/mfem/parallel/passemble.cpp
		  #${SRC}/mfem/parallel/locmat.cpp
		  #${SRC}/mfem/parallel/psolve.cpp
		  #${SRC}/mfem/parallel/preconst.cpp
		  ${SRC}/mfem/clean/serial_util.cpp
		  ${SRC}/mfem/clean/serial_bndry.cpp
		  ${SRC}/mfem/clean/serial_locmat.cpp
		  ${SRC}/mfem/clean/serial_assemble.cpp
)

target_include_directories(mfem_lib
		  PUBLIC ${SRC}/mfem
		  PUBLIC ${SRC}/mfem/clean
)

target_link_libraries(mfem_lib
    phase_lib
    general_lib
)

# ===================================================================
add_library(couple_lib
		  ${SRC}/corner/couple.cpp
		  ${SRC}/corner/print.cpp
		  porosity.cpp
		  flow.cpp
		  )

target_include_directories(couple_lib
		  PUBLIC ${SRC}/phase
		  PUBLIC ${SRC}/general
		  PUBLIC ${SRC}/mfem
		  PUBLIC ${SRC}/mfem/clean
		  PUBLIC ${SRC}/corner
)

target_link_libraries(couple_lib
		  io_lib
		  general_lib
		  phase_lib
		  mfem_lib
		  PkgConfig::PETSC
)

add_executable(driver
		  ${SRC}/corner/driver.cpp
)

target_include_directories(driver
		  PUBLIC ${SRC}/corner
		  )

target_link_libraries(driver
		  io_lib
		  phase_lib
		  mfem_lib
		  couple_lib
		  PkgConfig::LAPACKE
		  PkgConfig::PETSC
		  phase_lib
)
