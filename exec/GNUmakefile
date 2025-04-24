#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..

# Data type selection (float, double, complex_float, complex_double)
DATATYPE ?= double

# We've made a common set of makefile definitions and rules in the home directory.
# By including that file, you have all the definitions you need.

include $(HOME)/GNUmakefile
CFLAGS += -O3
#CFLAGS += -g


# Add data type specific flags
ifeq ($(DATATYPE),float)
    CFLAGS += -DDATATYPE=float
else ifeq ($(DATATYPE),double)
    CFLAGS += -DDATATYPE=double
else ifeq ($(DATATYPE),complex_float)
    CFLAGS += -DDATATYPE=complex_float
else ifeq ($(DATATYPE),complex_double)
    CFLAGS += -DDATATYPE=complex_double
else
    $(error Invalid DATATYPE. Use one of: float, double, complex_float, complex_double)
endif

CXX = g++
# CXX = clang++
# Set DIM dynamically on command line, e.g. make extrude DIM=3, default is 2D
DIM ?= 2
CPPFLAGS = -DDIM=$(DIM)


fesolver.exe: GNUmakefile $(OBJS) FEMain.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o fesolver.exe FEMain.cpp  $(OBJS)

run: fesolver.exe 
	./fesolver.exe ../FEData/A
	./fesolver.exe ../FEData/B.2

# Set DIM as 3 on command line, i.e. make extrude DIM=3
extrude: fesolver.exe 
	./fesolver.exe ../FEData/B.2 --extrude

poly: fesolver.exe 
	./fesolver.exe ../FEData/test 2

clean:
	rm -fr $(OBJS) $(OBJS:.o=.d) *.exe *.dSYM *.vtk *.o *.d

# Add targets for different data types
run_float:
	$(MAKE) DATATYPE=float run

run_double:
	$(MAKE) DATATYPE=double run

run_complex_float:
	$(MAKE) DATATYPE=complex_float run

run_complex_double:
	$(MAKE) DATATYPE=complex_double run

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d
