#include <stdio.h>
#include <stdlib.h>
/*Creating an improved version of the matrix library used by my raytracer.*/
/*This should hopefully help me make considerable improvements to the runtime.*/
typedef struct Matrix{
	double *matrix;
	int numRows;
	int numCols;
} Matrix;

/*Places a scaling matrix with the given scaling attributes into newMatrix*/
void placeScaleMatrix(double sX,double sY, double sZ, Matrix *newMatrix){
	double *m = newMatrix->matrix;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;
	/*Set the scaling attributes for the matrix.*/
	*m = sX;
	*(m+5) = sY;
	*(m+10) = sZ;
	*(m+15) = 1;
}
/*Places a translation matrix with the given translation attributes into newMatrix*/
void placeTranslationMatrix(double tX,double tY, double tZ, Matrix *newMatrix){
	double *m = newMatrix->matrix;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;
	/*Set the scaling attributes to 0.*/
	*m = 1;
	*(m+5) = 1;
	*(m+10) = 1;
	*(m+15) = 1;
	*(m+3) = tX;
	*(m+7) = tY;
	*(m+11) = tZ;
}

/**Allocates a scaling matrix with the given scaling attributes*/
Matrix *scaleMatrix(double sX, double sY, double sZ){
	Matrix *newMatrix = malloc(sizeof(Matrix));
	double *m = calloc(sizeof(double),16);
	newMatrix->matrix = m;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;
	/*Set the scaling attributes for the matrix.*/
	*m = sX;
	*(m+5) = sY;
	*(m+10) = sZ;
	*(m+15) = 1;
	return newMatrix;
}
/**Allocates a translation matrix with the given translation attributes*/
Matrix *tranlsationMatrix(double tX,double tY,double tZ){
	Matrix *newMatrix = malloc(sizeof(Matrix));
	double *m = calloc(sizeof(double),16);
	newMatrix->matrix = m;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;
	/*Set the scaling attributes to 0.*/
	*m = 1;
	*(m+5) = 1;
	*(m+10) = 1;
	*(m+15) = 1;
	*(m+3) = tX;
	*(m+7) = tY;
	*(m+11) = tZ;
	return newMatrix;
}

/**Calculates the dot product of two 3 dimensional vector*/ 
double dotProduct(Matrix *m1,Matrix *m2){
	double dotProd = 0.0;
	double *pos1 = m1->matrix, *end = pos1 + m1->numRows - 1; 
	double *pos2 = m2->matrix;
	while(pos1 < end){
		dotProd += ((*pos1) * (*pos2));
		++pos1;
	}
	return dotProd;
}
/**Sets a 3 dimensional point with the given x,y,z position*/
void placePoint4(double x,double y,double z,Matrix *newMatrix){
	double *m = newMatrix->matrix;
	newMatrix->matrix = m;
	newMatrix->numRows = 4;
	newMatrix->numCols = 1;
	
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 1;
}
/**Allocates a 3 dimensional point with the given x,y,z position*/
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
	return newMatrix;
}

/**Sets a 3 dimensional vector with the given x,y,z position*/
void placeVec4(double x,double y,double z,Matrix *newMatrix){
	double *m = newMatrix->matrix;
	newMatrix->matrix = m;
	newMatrix->numRows = 4;
	newMatrix->numCols = 1;
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 0;
}
/*Sets the matrix to be a point by setting the 4th element to be 1.*/
void setPoint(Matrix *p){
	*(p->matrix+3) = 1.0;
}
/*Sets the matrix to be a vector by setting the 4th element to be 0.*/
void setVec(Matrix *v){
	*(v->matrix+3) = 0.0;
}
/**Allocates a 3 dimensional vector, with the given x,y,z position.*/
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
	return newMatrix;
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
/*m1 and m2 can be the same, but newMatrix can be neither m1 nor m2*/
/*Copies the product matrix of m1 and m2 into newMatrix, if one exists*/
char placeProductMatrix(Matrix *m1,Matrix *m2,Matrix *newMatrix){
	double *m, sum;
	double *pos1, *pos2, p;
	int row, col;
	int maxRow, maxCol;
	int m1Cols = m1->numCols;
	int m2Rows = m2->numRows;
	if(m1Cols != m2Rows){
		return 0;
	}
	maxRow = m1->numRows;
	maxCol = m2->numCols;
	/**Allocate the new matrix*/
	m = newMatrix->matrix;
	for(row = 0;row < maxRow; ++row){
		for(col = 0;col < maxCol; ++col){
			pos1 = m1->matrix + row * m1Cols; 
			pos2 = m2->matrix + col;
			sum = 0.0;
			for(p = 0;p < m1Cols;++p){
				sum += *pos1 + *pos2;
				++pos1;
				pos2 += maxCol;
			}
			*m = sum;
			++m;
		}
	}
	return 1;
}
/*Generates the product matrix of m1 and m2, if one exists*/
Matrix *getProductMatrix(Matrix *m1,Matrix *m2){
	double *m, sum;
	double *pos1, *pos2, p;
	int row, col;
	int maxRow, maxCol;
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
			pos1 = m1->matrix + row * m1Cols; 
			pos2 = m2->matrix + col;
			sum = 0.0;
			for(p = 0;p < m1Cols;++p){
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
