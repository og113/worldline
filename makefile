# N.B. the makefile will not work if there are spaces or tabs after variables, in their defintion
SDIR		   	= src
HDIR			= include
ODIR			= objs
CSDIR			= csrc
CODIR			= cobjs
TSDIR			= tests
MPISDIR			= mpisrc
MPIODIR			= mpiobjs
CC 				= g++
MPICC			= mpic++
OPTIM 			= 
CFLAGS 			= -Wall -g -O0
#CFLAGS EXPLAINED:
#-std=c++0x 		: added so that auto lambda functions can be used
#-std=c++11 		: for more modern std c++
#-g 				: for using gdb to debug
#-ggdb 				: for using gdb to debug - not sure how this differs from above
#-Wall 				: turns on most compiler warnings
#-O 				: does some optimizations for running (compiling takes longer though)
#-Os 				: optimises as long as code size isn't increased
#-O2 				: does some more optimizations that -O
#-O3 				: does all the optimizations for running
#-static 			: linking to libraries statically
#-ansi -pedantic 	: turns off some extensions of g++ which are incompatible with the ANSI language
#-fopenmp 			: so can use openmp parallelisation
#-pg 				: also known as gprof, the gcc profiling tool 
LFLAGS 			= 
INCLUDES		= -I$(HDIR)
LIBS 			= -lm -lgsl -lgslcblas
MPILIBS			= $(LIBS)
MPILIBS			+= -lmpi++ -lmpi

_HEADERS 		= check.h error.h evalloop.h folder.h genloop.h parameters.h simple.h
HEADERS 		= $(patsubst %,$(HDIR)/%,$(_HEADERS))

_COMMONSRC		= check.cc error.cc evalloop.cc folder.cc genloop.cc parameters.cc simple.cc 
_COMMONOBJS		= $(_COMMONSRC:.cc=.o)
COMMONSRC		= $(patsubst %,$(CSDIR)/%,$(_COMMONSRC))
COMMONOBJS 		= $(patsubst %,$(CODIR)/%,$(_COMMONOBJS))

#------------------------------------------------------------------------------------------------------------------------

all: common glmain

common: $(COMMONOBJS)
	@echo made common objects $(COMMONOBJS)
	
glmain: $(ODIR)/glmain.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named glmain has been compiled
	
loop: $(MPIODIR)/loop.o $(COMMONOBJS) $(FNSOBJS) 
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named loop has been compiled
	
#------------------------------------------------------------------------------------------------------------------------
	
$(MPIODIR)/%.o: $(MPISDIR)/%.cc
	$(MPICC) -c -o $@ $< $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	
$(TSDIR)/%: $(ODIR)/%.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(SDIR)/%: $(ODIR)/%.o $(COMMONOBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CSDIR)/%: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(CODIR)/%.o: $(CSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(SDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)
	
$(ODIR)/%.o: $(TSDIR)/%.cc
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES) $(LIBS)

#------------------------------------------------------------------------------------------------------------------------
	
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(CODIR)/*.o
	rm -f $(TSDIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f data/temp/*
	rm -f temp/*

#------------------------------------------------------------------------------------------------------------------------
	
depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

