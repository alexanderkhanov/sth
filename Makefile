CPP=g++
LD=g++
SPECIALFLAGS=-O
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs)

CFLAGS = $(SPECIALFLAGS) -I.
LFLAGS = $(SPECIALFLAGS) -L. -lMinuit

RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=d_ana.cpp

%.o: %.cpp matrix.h eval.h lepton_masses.h resval.h minimizer.h output.h
	$(CPP) $(RCXX) -c $<

all: $(SRC:.cpp=.o)
	$(LD) $(SRC:.cpp=.o) $(RLXX) -o d_ana

clean:
	rm -f *~ *.o
	rm -f d_ana
