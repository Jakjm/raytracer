CC=gcc

raytracer.exe: main.o matrix2.o Makefile
	$(CC) -O3 -pthread -o raytracer.exe main.o matrix2.o -lm -g -Wall 
matrix2.o : matrix2.c matrix2.h
	$(CC) -c -O3 -pthread matrix2.c -lm -g -Wall
main.o : main.c matrix2.h
	$(CC) -c -O3 -pthread main.c -lm -g -Wall
