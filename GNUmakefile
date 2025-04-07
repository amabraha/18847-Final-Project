#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

DIM = 2
FE_SRC = $(HOME)/src/fe
SPMAT_SRC = $(HOME)/src/spmat
VISIT_SRC = $(HOME)/VisitWriter
TRIANGLE_SRC = $(HOME)/triangle
VPATH = $(HOME) $(VISIT_SRC) $(FE_SRC) $(SPMAT_SRC) $(TRIANGLE_SRC)
CFLAGS = -Wall -I$(FE_SRC) -I$(VISIT_SRC) -I$(SPMAT_SRC) -I$(TRIANGLE_SRC) -std=c++11

CXX = g++
# CXX = clang++
CPPFLAGS = -D DIM=$(DIM) 

#ifeq ($(CXX) , clang++)
#CFLAGS += -stdlib=libc++
#endif

# Triangle related
TRILIBDEFS = -DTRILIBRARY
CC = gcc

FE_SRCFILES:= $(wildcard $(FE_SRC)/*.cpp)
SPMAT_SRCFILES:= $(wildcard $(SPMAT_SRC)/*.cpp)
FE_OBJS:=$(patsubst %.cpp, %.o, $(FE_SRCFILES))
SPMAT_OBJS:=$(patsubst %.cpp, %.o, $(SPMAT_SRCFILES))
VISIT_OBJS := $(VISIT_SRC)/VisitWriter.o
TRIANGLE_OBJS := $(TRIANGLE_SRC)/triangle.o
OBJS = $(FE_OBJS) $(SPMAT_OBJS) $(VISIT_OBJS) $(TRIANGLE_OBJS)

$(TRIANGLE_SRC)/triangle.o: $(TRIANGLE_SRC)/triangle.c $(TRIANGLE_SRC)/triangle.h
	$(CC) $(TRILIBDEFS) -c -o $(TRIANGLE_SRC)/triangle.o \
		$(TRIANGLE_SRC)/triangle.c

%.o:%.cpp GNUmakefile
	$(CXX) -c $(CPPFLAGS) $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

