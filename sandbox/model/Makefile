CFLAGS = -g
CPROTO = /usr/bin/cproto
CTAGS = /usr/bin/exuberant-ctags

model.mpi: model.c global_defs.h prototypes.h ../parsing/Parsing.o ctags
	mpicc $(CFLAGS) model.c ../parsing/Parsing.o -o $@

../parsing/Parsing.o: ../parsing/Parsing.c ../parsing/Parsing.h
	mpicc $(CFLAGS) -c ../parsing/Parsing.c -o $@

prototypes.h: model.c
	$(CPROTO) -D__CPROTO__ -q -p -f 3 model.c > prototypes.h 2>/dev/null

ctags:
	$(CTAGS) model.c global_defs.h

clean:
	rm -f model.o ../parsing/Parsing.o
	rm -f prototypes.h
	rm -f tags

.PHONY: clean ctags
