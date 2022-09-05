CC = gcc

objects = bin_interp.o

bin_interp.a: $(objects)
	ar rcs bin_interp.a $(objects)

bin_interp.o: bin_interp.h
	CC -c -O3 bin_interp.c

clean:
	rm $(objects) bin_interp.a