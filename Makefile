CC=gcc
flags=-g -lm -Wall `pkg-config --cflags --libs MagickWand` -O3
raytracer.exe: main.o matrix2.o Makefile
	$(CC) -pthread -o raytracer.exe main.o matrix2.o $(flags)
matrix2.o : matrix2.c matrix2.h Makefile
	$(CC) -c matrix2.c $(flags)
main.o : main.c matrix2.h Makefile
	$(CC) -c -pthread main.c $(flags) 
clean : 
	rm *.o raytracer.exe
