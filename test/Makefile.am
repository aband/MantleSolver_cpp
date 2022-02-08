bin_PROGRAMS = testweno
testweno_SOURCES = testweno.cpp ../src/weno/weno_basis.cpp ../src/mesh.c ../src/integral.c ../include/weno.h ../include/mesh.h

include $(PETSC_DIR)/lib/petsc/conf/variables

AM_LDFLAGS = -L${ADOLC_DIR}/lib64 ${PETSC_LIB_BASIC} -L${PETSC_LIB_DIR}  -ladolc -lm

AM_CPPFLAGS  = -DINCLUDE_PETSC $(PETSC_CC_INCLUDES)  ${PETSC_CCPPFLAGS} -I${ADOLC_DIR}/include

AM_CFLAGS = -O0
