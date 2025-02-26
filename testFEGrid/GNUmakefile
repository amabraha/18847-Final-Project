#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..

include $(HOME)/GNUmakefile
CFLAGS += -O3
CFLAGS += -g

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

ShowMesh.exe: GNUmakefile ShowMesh.cpp $(VISIT_OBJS) 
	$(CXX) $(CFLAGS) -o ShowMesh.exe ShowMesh.cpp $(VISIT_OBJS)

testFEGrid.exe: $(OBJS) testFEGrid.cpp $(OBJS)
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o testFEGrid.exe testFEGrid.cpp $(OBJS) 

clean:
	rm -fr $(OBJS) $(OBJS:.o=.d) *.exe *.dSYM *.vtk *.o *.d

