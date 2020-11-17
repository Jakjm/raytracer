CC=gcc

raytracer.exe: main.c matrix2.c matrix2.h
	$(CC) -O3 -pthread -o raytracer.exe main.c matrix2.c -lm -Wall -g
