CC=gcc

raytracer.exe: main.c matrix2.c doubleMatrix.h
	$(CC) -O3 -pthread -o raytracer.exe main.c matrix2.c -lm -Wall -g
