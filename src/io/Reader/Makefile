#TARG = common
CODE = reader_$(TARG).a

# Set-up information ==========================================================

HOME = ..
include $(HOME)/Config/Make.cf.$(TARG)

# Files =======================================================================

SRC = expr.cpp reader.cpp
INC = dimen.h expr.h reader.h units.h units.inc

OBJ = $(patsubst %.cpp, %.o, $(SRC))

# Main Rules ==================================================================

.SUFFIXES:

.SUFFIXES: .o .cpp

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $<

all: $(CODE)

$(CODE): $(INC) $(SRC)
	@echo 'LIBRARY: Reader'
	@make $(OBJ)
	$(AR) $(ARFLAGS) $(CODE) $(OBJ)
	@rm -f $(OBJ)

# Miscellaneous Rules =========================================================

clean: 
	rm -f *.o *~ $(CODE)

allclean: clean
	rm -f *.a

