CC=gcc

a.out: main.c doubleMatrix.c
	$(CC) -O3 -o raytracer.exe main.c -lm
