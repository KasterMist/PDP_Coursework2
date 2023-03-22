SRC = src/simulation_configuration.c src/main.c src/simulation_support.c
LFLAGS=-lm
CFLAGS=-O3 -fopenmp

# .PHONY: archer2 local build
.PHONY: cirrus local build

archer2: CC=cc
archer2: build

cirrus: CC=gcc
cirrus: build

build: 
	$(CC) -o reactor $(SRC) $(CFLAGS) $(LFLAGS)

