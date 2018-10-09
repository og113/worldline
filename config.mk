# configuration file for worldline programs on laptop

# N.B. the makefile will not work if there are spaces or tabs after variables, in their defintion
SDIR		   	= src
HDIR			= include
ODIR			= objs
LDIR			= libs
CSDIR			= csrc
CODIR			= cobjs
TSDIR			= tests
MPISDIR			= mpisrc
MPIODIR			= mpiobjs
CC 				= g++
MPICC			= mpic++ # other options include mpicxx, mpicxx.mpich
OPTIM 			= 
CFLAGS 			= -Wall -g -O0
MPICFLAGS		= -Wall -g -O3
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
INCLUDES		= -I$(HDIR) -L$(LDIR) -I/usr/include/eigen3 -L/usr/lib/x86_64-linux-gnu/ #-I../c++/eigen-eigen-1306d75b4a21
LIBS 			= -lm -lgsl -lgslcblas -ljbnumlib
MPILIBS			= $(LIBS)
MPILIBS			+= -lmpi # -lmpi++ deprecated
