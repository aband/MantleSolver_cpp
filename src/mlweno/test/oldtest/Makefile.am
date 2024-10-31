bin_PROGRAMS = testweno
testweno_SOURCES = testweno.cpp ../src/weno/weno_basis.cpp ../src/integral.cpp ../src/mesh.c ../include/weno.h ../include/integral.h ../include/mesh.h

AM_LDFLAGS = -L${ADOLC_DIR}/lib64 -lm -ladolc -lz -fstack-protector

AM_CPPFLAGS = -I${ADOLC_DIR}/include -O0

AM_CC = mpicxx

include ${PETSC_DIR}/lib/petsc/conf/variables

# Add petsc 

AM_LDFLAGS +=  ${PETSC_LIB_BASIC} ${PETSC_SNES_LIB} -L${PETSC_LIB_DIR} -lpetsc -lpthread

AM_CPPFLAGS += -DINCLUDE_PETSC $(PETSC_CC_INCLUDES) ${PETSC_CCPPFLAGS}
