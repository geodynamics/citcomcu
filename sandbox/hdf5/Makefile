empty: empty.c
	gcc $^ -o $@ -lhdf5

pytables.o: pytables.c pytables.h
	gcc -c pytables.c -o $@

2d-array: 2d-array.c pytables.o
	gcc $^ -o $@ -lhdf5

3d-array: 3d-array.c pytables.o
	gcc $^ -o $@ -lhdf5

hyperslab: hyperslab.c pytables.o
	mpicc $^ -o $@ -lhdf5
