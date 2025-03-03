#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..
DIM = 2
VPATH = . ../VisitWriter ../src/spmat
CFLAGS = -Wall -I. -I ../VisitWriter -I ../src/spmat -std=c++11
CFLAGS += -O3
#CFLAGS += -g
#CXX = g++
CXX = clang++

CPPFLAGS = -D DIM=$(DIM) 

ifeq ($(CXX) , clang++)
CFLAGS += -stdlib=libc++
endif


all: TestJacobi.exe

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d


TestJacobi.exe: TestJacobi.cpp ../src/spmat/SparseMatrix.o ../src/spmat/JacobiSolver.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o TestJacobi.exe TestJacobi.cpp ../src/spmat/SparseMatrix.o ../src/spmat/JacobiSolver.o

run: TestJacobi.exe
	./TestJacobi.exe 10
	./TestJacobi.exe 25
	./TestJacobi.exe 46
	./TestJacobi.exe 150

clean:
	rm -fr *.o *.d *.exe *.dSYM

-include $(OBJS:.o=.d)
