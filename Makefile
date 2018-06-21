all: curved-domain

CXXFLAGS=-g -Wall -pedantic -std=c++17
INCLUDES=-I../transfinite/src/geom -I../rust/harmonic/src
LIBS=-L../transfinite/debug/geom -lgeom -L../rust/harmonic/target -lharmonic

curved-domain: curved-domain.cc
	g++ $(CXXFLAGS) $(INCLUDES) -o $@ $< $(LIBS)
