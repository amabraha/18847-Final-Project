#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

DIM = 2
FE_SRC = $(HOME)/src/fe
SPMAT_SRC = $(HOME)/src/spmat
VISIT_SRC = $(HOME)/VisitWriter
VPATH = $(HOME) $(VISIT_SRC) $(FE_SRC) $(SPMAT_SRC)
CFLAGS = -Wall -I$(FE_SRC) -I$(VISIT_SRC) -I$(SPMAT_SRC) -std=c++11

#CXX = g++
CXX = clang++
CPPFLAGS = -D DIM=$(DIM) 

#ifeq ($(CXX) , clang++)
#CFLAGS += -stdlib=libc++
#endif

FE_SRCFILES:= $(wildcard $(FE_SRC)/*.cpp)
SPMAT_SRCFILES:= $(wildcard $(SPMAT_SRC)/*.cpp)
FE_OBJS:=$(patsubst %.cpp, %.o, $(FE_SRCFILES))
SPMAT_OBJS:=$(patsubst %.cpp, %.o, $(SPMAT_SRCFILES))
VISIT_OBJS := $(VISIT_SRC)/VisitWriter.o
OBJS = $(FE_OBJS) $(SPMAT_OBJS) $(VISIT_OBJS)

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

