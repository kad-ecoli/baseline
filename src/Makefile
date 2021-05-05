CC=g++
CFLAGS=-O3 -ffast-math
LDFLAGS=-static# -lm

all: NWalign SelfScore

SelfScore: SelfScore.cpp basic_fun.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp NWalign.h basic_fun.h BLOSUM.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: NWalign SelfScore
	cp NWalign SelfScore ../bin

clean:
	rm NWalign SelfScore
