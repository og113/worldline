# makefile for compiling programs to do worldline calculations
#------------------------------------------------------------------------------------------------------------------------
# config file

include config.mk

#------------------------------------------------------------------------------------------------------------------------
# local variables

HEADERS 		= $(wildcard $(HDIR)/*.h)
COMMONSRC		= $(wildcard $(CSDIR)/*.cc)
COMMONOBJS 		= $(patsubst $(CSDIR)/%.cc,$(CODIR)/%.o,$(COMMONSRC))

.PHONY : variables
variables :
	@echo HEADERS: $(HEADERS)
	@echo COMMONSRC: $(COMMONSRC)
	@echo COMMONOBJS: $(COMMONOBJS)

#------------------------------------------------------------------------------------------------------------------------
# some useful PHONYs

.PHONY: all
all: binaryToAscii circle common floop glmain loop loop2 nrmain schwingerRate

.PHONY: common
common: $(COMMONOBJS)

.PHONY: monte
monte: binaryToAscii circle common floop glmain loop loop2 schwingerRate

.PHONY: nr
nr: binaryToAscii circle common nrmain

#------------------------------------------------------------------------------------------------------------------------
# targets, dependencies and rules for executables

binaryToAscii: $(ODIR)/binaryToAscii.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named binaryToAscii has been compiled

circle: $(ODIR)/circle.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named circle has been compiled
	
floop: $(MPIODIR)/floop.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named floop has been compiled
	
gflmain: $(ODIR)/gflmain.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named gflmain has been compiled
	
glmain: $(MPIODIR)/glmain.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named glmain has been compiled
	
loop: $(MPIODIR)/loop.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named loop has been compiled
	
loop2: $(MPIODIR)/loop2.o $(COMMONOBJS)
	$(MPICC) -o $@ $^ $(MPICFLAGS) $(INCLUDES) $(MPILIBS)
	@echo Simple compiler named loop2 has been compiled
	
nrmain: $(ODIR)/nrmain.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named nrmain has been compiled
	
schwingerRate: $(ODIR)/schwingerRate.o $(COMMONOBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(INCLUDES) $(LIBS)
	@echo Simple compiler named schwingerRate has been compiled
	
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
# makedepend

.PHONY: depend	
depend: $(SRC)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

