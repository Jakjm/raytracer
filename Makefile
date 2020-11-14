CC=gcc

raytracer.exe: main.c doubleMatrix.c doubleMatrix.h
	$(CC) -O3 -pthread -o raytracer.exe doubleMatrix.c main.c -lm -Wall -g
