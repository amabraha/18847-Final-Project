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

#ifeq ($(CXX) , clang++)
#CFLAGS += -stdlib=libc++
#endif


all: TestMatrix.exe

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d


TestMatrix.exe: TestMatrix.cpp ../src/spmat/SparseMatrix.o
	$(CXX) $(CPPFLAGS) $(CFLAGS) -o TestMatrix.exe TestMatrix.cpp ../src/spmat/SparseMatrix.o

run: TestMatrix.exe
	./TestMatrix.exe 10
	./TestMatrix.exe 25
	./TestMatrix.exe 46
	./TestMatrix.exe 150

clean:
	rm -fr *.o *.d *.exe *.dSYM

-include $(OBJS:.o=.d)
