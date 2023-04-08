MF=Makefile

CC=mpicc

CFLAGS=	-O3 -fopenmp -lm

SRC = src/simulation_configuration.c src/main.c src/simulation_support.c src/mpi_support.c
LFLAGS= $(CFLAGS)

# .PHONY: archer2 local build
.PHONY: cirrus local build

archer2: CC=mpicc
archer2: build

cirrus: CC=mpicc
cirrus: build

build: 
	$(CC) -o reactor $(SRC) $(CFLAGS) $(LFLAGS)

# MF=	Makefile

# CC=	mpicc
# CFLAGS=	-O3 -Wall -fopenmp

# LFLAGS= $(CFLAGS)

# EXE=	reactor

# INC= \
# 	src/simulation_configuration.h \
# 	src/simulation_support.h

# SRC= \
# 	src/main.c \
# 	src/simulation_configuration.c \
# 	src/simulation_support.c

# archer2: CC=mpicc
# archer2: all

# cirrus: CC=mpicc
# cirrus: all

# #
# # No need to edit below this line
# #

# .SUFFIXES:
# .SUFFIXES: .c .o

# OBJ=	$(SRC:.c=.o)

# .c.o:
# 	$(CC) $(CFLAGS) -c $<

# all:	$(EXE)

# $(OBJ):	$(INC)

# $(EXE):	$(OBJ)
# 	$(CC) $(LFLAGS) -o $@ $(OBJ)

# $(OBJ):	$(MF)

# clean:
# 	rm -f $(EXE) $(OBJ) core