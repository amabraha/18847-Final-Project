#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

FE_SRC = $(HOME)/src/fe
SPMAT_SRC = $(HOME)/src/spmat
VISIT_SRC = $(HOME)/VisitWriter
TIME_DEPENDENT_SRC = $(HOME)/src/fe_td
TRIANGLE_SRC = $(HOME)/triangle
VPATH = $(HOME) $(VISIT_SRC) $(FE_SRC) $(SPMAT_SRC) $(TIME_DEPENDENT_SRC) $(TRIANGLE_SRC)
CFLAGS = -Wall -I$(FE_SRC) -I$(VISIT_SRC) -I$(SPMAT_SRC) -I$(TIME_DEPENDENT_SRC) -I$(TRIANGLE_SRC) -std=c++11

#ifeq ($(CXX) , clang++)
#CFLAGS += -stdlib=libc++
#endif

# Triangle related
TRILIBDEFS = -DTRILIBRARY
CC = gcc

FE_SRCFILES:= $(wildcard $(FE_SRC)/*.cpp)
SPMAT_SRCFILES:= $(wildcard $(SPMAT_SRC)/*.cpp)
TIME_DEPENDENT_SRCFILES := $(wildcard $(TIME_DEPENDENT_SRC)/*.cpp)
FE_OBJS:=$(patsubst %.cpp, %.o, $(FE_SRCFILES))
SPMAT_OBJS:=$(patsubst %.cpp, %.o, $(SPMAT_SRCFILES))
TIME_DEPENDENT_OBJS:=$(patsubst %.cpp, %.o, $(TIME_DEPENDENT_SRCFILES))
VISIT_OBJS := $(VISIT_SRC)/VisitWriter.o
TRIANGLE_OBJS := $(TRIANGLE_SRC)/triangle.o
OBJS = $(FE_OBJS) $(SPMAT_OBJS) $(VISIT_OBJS) $(TIME_DEPENDENT_OBJS) $(TRIANGLE_OBJS)

$(TRIANGLE_SRC)/triangle.o: $(TRIANGLE_SRC)/triangle.c $(TRIANGLE_SRC)/triangle.h
	$(CC) $(TRILIBDEFS) -g -c -o $(TRIANGLE_SRC)/triangle.o \
		$(TRIANGLE_SRC)/triangle.c

