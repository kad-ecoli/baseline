CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static# -lm

all: NWalign SelfScore

SelfScore: SelfScore.cpp basic_fun.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp NWalign.h basic_fun.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: NWalign
	cp NWalign ../bin

clean:
	rm NWalign
