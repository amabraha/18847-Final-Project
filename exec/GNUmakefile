# exec/Makefile

HOME = ./..
DOXYGEN_OUTPUT_DIR = $(HOME)/docs

# include the common rules from project root
include $(HOME)/GNUmakefile

# CFLAGS   += -O3
CFLAGS   += -g

CXX       = g++
#CXX     = clang++

# Set DIM dynamically on command line, default is 2D
DIM      ?= 2
CPPFLAGS  = -DDIM=$(DIM)

# Data type selection (float, double, complex_float, complex_double)
DATATYPE ?= double

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

#---------------------------------------------------------
# Original 2D solver (unchanged)
fesolver.exe: $(HOME)/GNUmakefile $(OBJS) FEMain.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ FEMain.cpp $(OBJS)

run: fesolver.exe
	./fesolver.exe ../FEData/A
	./fesolver.exe ../FEData/B.2

#---------------------------------------------------------
# Time-dependent solver
fesolver_TD.exe: GNUmakefile $(OBJS) FEMain_TD.cpp
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o fesolver_TD.exe FEMain_TD.cpp  $(OBJS)

run_TD: fesolver_TD.exe 
	rm -fr */*.vtk
	./fesolver_TD.exe ../FEData/td 0.05

#---------------------------------------------------------
# 3D extrusion
# Usage: make extrude DIM=3
extrude: fesolver.exe
	./fesolver.exe ../FEData/square 0.5

#---------------------------------------------------------
# Triangle integration
poly: fesolver.exe
	./fesolver.exe ../FEData/test 2

#--------------------------------------------------------------------
# 3D Convergence‐study targets (Part 10)
# Build the convergence driver object from FE3dConvergence.cpp
FE3dConvergence.o: FE3dConvergence.cpp
	$(CXX) $(CPPFLAGS) $(CFLAGS) -I$(HOME)/src/fe -I$(HOME)/src/spmat -c $< -o $@

# Link into the final 3D‐convergence executable
fe3d_convergence.exe: $(HOME)/GNUmakefile $(OBJS) FE3dConvergence.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ FE3dConvergence.o $(OBJS)

# Run the convergence study on square.poly
# Usage: make run3d DIM=3
run3d: fe3d_convergence.exe
	./fe3d_convergence.exe ../FEData/square 1 5

# ---------------------------------------------------------
# 2-D convergence study
fe2d_convergence.exe: GNUmakefile $(OBJS) FE2dConvergence.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o $@ FE2dConvergence.o $(OBJS)

# Run the convergence study on square.poly
run2d: fe2d_convergence.exe
	./fe2d_convergence.exe ../FEData/square 1.0 5

#---------------------------------------------------------
# convenience targets for data types
run_float:
	$(MAKE) DATATYPE=float run

run_double:
	$(MAKE) DATATYPE=double run

run_complex_float:
	$(MAKE) DATATYPE=complex_float run

run_complex_double:
	$(MAKE) DATATYPE=complex_double run

#---------------------------------------------------------
# generic compile rule for .o from .cpp
%.o: %.cpp $(HOME)/GNUmakefile
	$(CXX) $(CPPFLAGS) $(CFLAGS) -c $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

#---------------------------------------------------------
# doxygen, generate docs in ../docs (DOXYGEN_OUTPUT_DIR)
docs:
	doxygen $(HOME)/doxygen.config

docs_clean:
	rm -rf $(DOXYGEN_OUTPUT_DIR)
#---------------------------------------------------------
clean:
	rm -fr $(OBJS) $(OBJS:.o=.d) *.exe FE3dConvergence.o fe3d_convergence.exe FE2dConvergence.o FE2dConvergence.d  *.vtk
