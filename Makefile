CC=gcc

raytracer.exe: main.c doubleMatrix.c
	$(CC) -O3 -pthread -o raytracer.exe main.c -lm -g
