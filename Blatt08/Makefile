# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff partdiff-par-hybrid
OBJS = partdiff.o askparams.o displaymatrix-mpi.o

# Targets ...
all: $(TGTS)

partdiff: $(OBJS) Makefile
    $(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-par-hybrid: partdiff-par-hybrid.o askparams.o displaymatrix-mpi.o Makefile
    $(CC) $(LFLAGS) -fopenmp -o $@ partdiff-par-hybrid.o askparams.o displaymatrix-mpi.o -lm

partdiff.o: partdiff.c Makefile

partdiff-par-hybrid.o: Makefile
    $(CC) -c $(CFLAGS) -fopenmp -o partdiff-par-hybrid.o partdiff.c

askparams.o: askparams.c Makefile

displaymatrix-mpi.o: displaymatrix-mpi.c Makefile

# Rule to create *.o from *.c
%.o: %.c
    $(CC) -c $(CFLAGS) $*.c

clean:
    $(RM) $(OBJS)
    $(RM) $(TGTS)