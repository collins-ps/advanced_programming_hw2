CC=gcc
CFLAGS=-lm -I.
DEPS = timer.h
OBJ_nbody = timer.o nbody.o 
OBJ_nbody_omp = timer.o nbody_omp.o 
OBJ_nbody_pt = timer.o nbody_pt.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

nbody: $(OBJ_nbody)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	\rm -f *.o nbody *~ *#

nbody_pt: $(OBJ_nbody_pt)
	$(CC) -pthread -o $@ $^ $(CFLAGS)

nbody_omp: nbody_omp.c timer.c timer.h
	gcc -fopenmp -o nbody_omp nbody_omp.c timer.c -lm -I.