CC=g++ -O3 -Wall
sources=./src/GSGBN.cpp ./src/solve_GSGBN.cpp ./src/solve_GSGBN.h

GSGBN: $(sources)
	$(CC) -o ./bin/GSGBN $(sources) 

