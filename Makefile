# makefile for compiling programs to do worldline calculations
#------------------------------------------------------------------------------------------------------------------------
# config file

include config.mk

#------------------------------------------------------------------------------------------------------------------------
# local variables

HEADERS 		= $(wildcard $(HDIR)/*.h)
COMMONSRC		= $(wildcard $(CSDIR)/*.cc)
COMMONOBJS 		= $(patsubst $(CSDIR)/%.cc,$(CODIR)/%.o,$(COMMONSRC))
SRC				= $(wildcard $(SDIR)/*.cc)
EXE				= $(patsubst $(SDIR)/%.cc,%,$(SRC))
MPISRC			= $(wildcard $(MPISDIR)/*.cc)
MPIEXE			= $(patsubst $(MPISDIR)/%.cc,%,$(MPISRC))

.PHONY : variables
variables :
	@echo HEADERS: $(HEADERS)
	@echo COMMONSRC: $(COMMONSRC)
	@echo COMMONOBJS: $(COMMONOBJS)
	@echo SRC: $(SRC)
	@echo EXE: $(EXE)
	@echo MPISRC: $(MPISRC)
	@echo MPIEXE: $(MPIEXE)

#------------------------------------------------------------------------------------------------------------------------
# some useful PHONYs

.PHONY: all
all: $(EXE) $(MPIEXE) common

.PHONY: common
common: $(COMMONOBJS)

.PHONY: monte
monte: addVectors binaryToAscii common floop glmain loop loop2 schwingerRate shape

.PHONY: nr
nr: 3dPotentialExtrema addVectors binaryToAscii common dimReduce highTemp highTempParameters nrmain nrmpi perturbativeFiniteTemp shape

#------------------------------------------------------------------------------------------------------------------------
# targets, dependencies and rules for executables

3dPotentialExtrema: $(ODIR)/3dPotentialExtrema.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named 3dPotentialExtrema has been compiled

addVectors: $(ODIR)/addVectors.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named addVectors has been compiled

binaryToAscii: $(ODIR)/binaryToAscii.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named binaryToAscii has been compiled
	
dimReduce: $(ODIR)/dimReduce.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named dimReduce has been compiled
	
dotVectors: $(ODIR)/dotVectors.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named dotVectors has been compiled
	
floop: $(MPIODIR)/floop.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named floop has been compiled
	
gflmain: $(ODIR)/gflmain.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named gflmain has been compiled
	
glmain: $(MPIODIR)/glmain.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named glmain has been compiled
	
highTemp: $(ODIR)/highTemp.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named highTemp has been compiled
	
highTempParameters: $(ODIR)/highTempParameters.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named highTempParameters has been compiled
	
loop: $(MPIODIR)/loop.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named loop has been compiled
	
loop2: $(MPIODIR)/loop2.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named loop2 has been compiled
	
nrmain: $(ODIR)/nrmain.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named nrmain has been compiled
	
nrmpi: $(MPIODIR)/nrmpi.o $(COMMONOBJS) 
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named nrmpi has been compiled
	
perturbativeFiniteTemp: $(ODIR)/perturbativeFiniteTemp.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named perturbativeFiniteTemp has been compiled
	
schwingerRate: $(ODIR)/schwingerRate.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named schwingerRate has been compiled

shape: $(ODIR)/shape.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named shape has been compiled
	
testMpi: $(MPIODIR)/testMpi.o
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named testMpi has been compiled
	
testMpi2: $(MPIODIR)/testMpi2.o
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named testMpi2 has been compiled
	
#------------------------------------------------------------------------------------------------------------------------
# generic rules
	
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
# clean
	
.PHONY: clean
clean:
	rm -f $(ODIR)/*.o
	rm -f $(CODIR)/*.o
	rm -f $(MPIODIR)/*.o
	rm -f $(TSDIR)/*.o
	rm -f *~ core
	rm -f $(HDIR)/*~
	rm -f data/temp/*
	rm -f temp/*

#------------------------------------------------------------------------------------------------------------------------
# makedepend, NOT WORKING

.PHONY: depend	
depend: $(COMMONSRC) $(SRC)
	makedepend -- $(CFLAGS) -- $^ -f- > Makefile.deps

# DO NOT DELETE THIS LINE -- make depend needs it
