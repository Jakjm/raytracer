#include <stdio.h>
#include <stdlib.h>
/*Creating an improved version of the matrix library used by my raytracer.*/
/*This should hopefully help me make considerable improvements to the runtime.*/
typedef struct Matrix{
	double *matrix;
	int numRows;
	int numCols;
} Matrix;

/**Allocates a 4 dimensional point with the given x,y,z position*/
Matrix *point4(double x,double y,double z){
	double *m = malloc(sizeof(double) * 4);
	Matrix *newMatrix = malloc(sizeof(Matrix));
	newMatrix->numRows = 4;
	newMatrix->numCols = 1;
	newMatrix->matrix = m;
	
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 1;
}
/**Allocates a 4 dimensional vector, with the given x,y,z position.*/
Matrix *vec4(double x,double y,double z){
	double *m = malloc(sizeof(double) * 4);
	Matrix *newMatrix = malloc(sizeof(Matrix));
	newMatrix->numRows = 4;
	newMatrix->numCols = 1;
	newMatrix->matrix = m;
	
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 0;
}
/**
 * Adds m2 to m1. Returns 1 if successful;
 * If it wasn't possible, returns 0.*/
char inPlaceSum(Matrix *m1,Matrix *m2){
	double *pos1, *end;
	double *pos2;
	double sum;
	if(m1->numRows != m2->numRows || m1->numCols != m2->numCols){
		return 0;
	}
	pos1 = m1->matrix;
	end = pos1 + m1->numRows * m1->numCols;
	pos2 = m2->matrix;
	while(pos1 < end){
		sum = *pos1 + *pos2;
		*pos1 = sum;
		++pos1; 
		++pos2; 
	}	
	return 1;
}
/**
 * Subtracts m2 from m1. Returns 1 if successful;
 * If it wasn't possible, returns 0.*/
char inPlaceDifference(Matrix *m1,Matrix *m2){
	double *pos1, *end;
	double *pos2;
	double difference;
	if(m1->numRows != m2->numRows || m1->numCols != m2->numCols){
		return 0;
	}
	pos1 = m1->matrix;
	end = pos1 + m1->numRows * m1->numCols;
	pos2 = m2->matrix;
	while(pos1 < end){
		difference = *pos1 - *pos2;
		*pos1 = difference;
		++pos1; 
		++pos2; 
	}	
	return 1;
}
/* Multiplies m in place by c*/
void inPlaceScalarMultiply(Matrix *m, double c){
	double *pos, *end;
	double product;
	pos = m->matrix;
	end = pos + m->numRows * m->numCols;
	while(pos < end){
		product = *pos * c;
		*pos = product; 
		++pos;
	}
}
/*Generates the product matrix of m1 and m2, if one exists*/
Matrix *getProductMatrix(Matrix *m1,Matrix *m2){
	double *m, sum;
	double *pos1, *pos2, p;
	int row, col;
	int maxRow, maxCol;
	int numProducts;
	int m1Cols = m1->numCols;
	int m2Rows = m2->numRows;
	Matrix *newMatrix;
	if(m1Cols != m2Rows){
		return NULL;
	}
	maxRow = m1->numRows;
	maxCol = m2->numCols;
	/**Allocate the new matrix*/
	m = malloc(sizeof(double) * (maxRow * maxCol));
	for(row = 0;row < maxRow; ++row){
		for(col = 0;col < maxCol; ++col){
			pos1 = m1->matrix + row * numProducts; 
			pos2 = m2->matrix + col;
			sum = 0.0;
			for(p = 0;p < numProducts;++p){
				sum += *pos1 + *pos2;
				++pos1;
				pos2 += maxCol;
			}
			*m = sum;
			++m;
		}
	}
	/*Finalize the matrix structure that is returned.*/
	newMatrix = malloc(sizeof(Matrix));
	newMatrix->matrix = m;
	newMatrix->numRows = maxRow;
	newMatrix->numCols = maxCol;
	return newMatrix;
}
/**Frees the heap resources used by the given matrix*/
void freeMatrix(Matrix *m){
	free(m->matrix);
	free(m);
}
