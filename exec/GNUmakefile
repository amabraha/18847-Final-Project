#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..

# We've made a common set of makefile definitions and rules in the home directory.
# By including that file, you have all the definitions you need.

include $(HOME)/GNUmakefile
CFLAGS += -O3
#CFLAGS += -g

fesolver.exe: GNUmakefile $(OBJS) FEMain.cpp
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o fesolver.exe FEMain.cpp  $(OBJS)

run: fesolver.exe 
	./fesolver.exe ../FEData/A
	./fesolver.exe ../FEData/B.2

fesolver_TD.exe: GNUmakefile $(OBJS) FEMain_TD.cpp
	  $(CXX) $(CPPFLAGS) $(CFLAGS) -o fesolver_TD.exe FEMain_TD.cpp  $(OBJS)

run_TD: fesolver_TD.exe 
	./fesolver_TD.exe ../FEData/A

clean:
	rm -fr $(OBJS) $(OBJS:.o=.d) *.exe *.dSYM *.vtk *.o *.d

