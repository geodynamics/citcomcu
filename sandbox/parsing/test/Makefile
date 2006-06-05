test_parsing: test_parsing.c ../Parsing.o
	gcc -g $^ -o $@

../Parsing.o: ../Parsing.h ../Parsing.c
	gcc -g -c ../Parsing.c -o ../Parsing.o

clean:
	rm -f ../Parsing.o test_parsing

run:
	./test_parsing input1

run1:
	./test_parsing input1 on

run2:
	./test_parsing input1 off on

.PHONY: clean run run1 run2
