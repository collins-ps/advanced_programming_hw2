CC=gcc
CFLAGS=-I.
DEPS = timer.h
OBJ = timer.o nbody.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

nbody: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	\rm -f *.o nbody *~ *#
