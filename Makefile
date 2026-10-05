# Simple alternative to CMake.   make            -> build/squirmer_bem + tests
#                                make test       -> run the (quick) test suite
CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall
LIBS     ?= -lopenblas
INC       = -Isrc
HEADERS   = $(wildcard src/*.hpp src/utils/*.hpp)
TESTS     = test_duct_flow test_squirmer test_trajectory

all: build/squirmer_bem $(addprefix build/,$(TESTS))

build/squirmer_bem: src/main.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@ $(LIBS)

build/test_%: tests/test_%.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@ $(LIBS)

build:
	mkdir -p build

test: all
	./build/test_duct_flow && ./build/test_squirmer && ./build/test_trajectory

clean:
	rm -rf build

.PHONY: all test clean
