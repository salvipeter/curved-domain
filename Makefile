all: curved-domain

CXXFLAGS=-g -Wall -pedantic -std=c++17
INCLUDES=-I$(TRANSFINITE_ROOT)/src/geom -I$(HARMONIC_ROOT)/src
LIBS=-L$(TRANSFINITE_ROOT)/debug/geom -lgeom -L$(HARMONIC_ROOT)/target -lharmonic

curved-domain: curved-domain.cc
	g++ $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
