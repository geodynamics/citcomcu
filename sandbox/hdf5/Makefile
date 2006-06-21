empty: empty.c
	gcc $^ -o $@ -lhdf5

pytables.o: pytables.c pytables.h
	gcc -c pytables.c -o $@

2d-array: 2d-array.c pytables.o
	gcc $^ -o $@ -lhdf5
