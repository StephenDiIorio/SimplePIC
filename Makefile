RM=rm -rf
CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic -O3
LDFLAGS=-g -O3
LDLIBS=

TARGET=pic


all: main

SRCDIR=src
BINDIR=bin
OBJDIR=obj

SRCFILES=pic.cpp Simulation.cpp Particle.cpp Species.cpp Field.cpp FFT.cpp ThreeVec.cpp

OBJFILES:=$(SRCFILES:.cpp=.o)

FULLTARGET=$(BINDIR)/$(TARGET)

# Path to look for source files
VPATH=$(SRCDIR):$(OBJDIR)

main: $(FULLTARGET)


# Rule to build the cpp files
%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $(OBJDIR)/$@ $<

$(FULLTARGET): $(OBJFILES)
	mkdir -p $(BINDIR)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

clean:
	$(RM) $(OBJDIR)

cleanall: clean
	$(RM) $(BINDIR)

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

.PHONY: clean cleanall main


# All the dependencies
pic.o: pic.cpp gyro.h uh.h Simulation.h Species.h Particle.h ThreeVec.h Field.h FFT.h
	$(CXX) -c $(CXXFLAGS) -o $(OBJDIR)/$@ $<

Simulation.o: Simulation.cpp Simulation.h Species.h Particle.h ThreeVec.h Field.h FFT.h
Particle.o: Particle.cpp Particle.h ThreeVec.h
Species.o: Species.cpp Species.h Particle.h ThreeVec.h Field.h FFT.h
Field.o: Field.cpp Field.h FFT.h
FFT.o: FFT.cpp FFT.h
ThreeVec.o: ThreeVec.cpp ThreeVec.h
