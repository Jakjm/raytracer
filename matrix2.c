#include <stdio.h>
#include <stdlib.h>
#include "doubleMatrix.h"
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
Matrix *getInverseMatrix(Matrix *m){
	Matrix *inverse = malloc(sizeof(Matrix));
	inverse->matrix = malloc(sizeof(double) * (m->numRows * m->numRows));
	placeInverseMatrix(m,inverse);
	return inverse;
}
/*
 *Produces the inverse of the given matrix.
 */
void placeInverseMatrix(Matrix *m,Matrix *inverseMatrix){
	int row, col, pos, otherRow;
	/*Storing the working positions between two rows of the matrix*/
	double *row1, *rowEnd, *row2;
	double *irow1, *irow2;
	
	double swap;
	double scalar;
	int numRows = m->numRows;
	int numCols = m->numCols;
	/*Put an identity matrix in the place of the inverse matrix.*/
	placeScaleMatrix(1,1,1,inverseMatrix);

	/*Solving the given matrix.*/
	for(row = 0;row < numRows;++row){
		/*Look for nonzero row.*/
		for(col = row;row < numRows; ++row){
			pos = numCols*row + col;
			if(*(m->matrix + pos) != 0.0)break;
		}
		if(row >= numRows){
			return;
		}

		/*Swap the found row and the position where the row should go.*/
		row1 = m->matrix + row * numCols;
		row2 = m->matrix + col * numCols;
		rowEnd = row1 + numRows;

		irow1 = inverseMatrix->matrix + row * numCols;
		irow2 = inverseMatrix->matrix + col * numCols;
		
		/*Scalar divide the row by this value*/
		scalar = 1.0 / *(row1 + col);
		for(;row1 < rowEnd;++row1, ++row2, ++irow1, ++irow2){
			/*Swap the rows in the matrix and apply the scalar multiplication to ensure we have a 1 in that column for this row.*/
			swap = *row1;
			*row1 = *row2;
			*row2 = swap * scalar;
			
			/*Swap the rows in the soon-to-be inverse matrix and apply the scalar multiplication*/
			swap = *irow1;
			*irow1 = *irow2;
			*irow2 = swap * scalar;
		}
		row = col;
		rowEnd = m->matrix + row * numCols + numCols;
		for(otherRow = 0;otherRow < row;++otherRow){
			scalar = *(m->matrix + (numRows * otherRow) + col);
			
			row1 = m->matrix + row * numCols;
			row2 = m->matrix + otherRow * numCols; 
	
			irow1 = inverseMatrix->matrix + numRows*numCols;
			irow2 = inverseMatrix->matrix + numRows*numCols;
			for(;row1 < rowEnd;++row1, ++row2, ++irow1, ++irow2){
				*row2 -= scalar * (*row1);
				*irow2 -= scalar * (*irow1);
			}
		}
		for(otherRow = row + 1;otherRow < numRows;++otherRow){
			scalar = *(m->matrix + (numRows * otherRow) + col);
			
			row1 = m->matrix + row * numCols;
			row2 = m->matrix + otherRow * numCols; 
	
			irow1 = inverseMatrix->matrix + numRows*numCols;
			irow2 = inverseMatrix->matrix + numRows*numCols;
			for(;row1 < rowEnd;++row1, ++row2, ++irow1, ++irow2){
				*row2 -= scalar * (*row1);
				*irow2 -= scalar * (*irow1);
			}
		}
	}
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
Matrix *translationMatrix(double tX,double tY,double tZ){
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
void setPoint4(Matrix *matrix, double x,double y, double z){
	double *m = matrix->matrix;
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 1;
}
void setVec4(Matrix *matrix, double x,double y, double z){
	double *m = matrix->matrix;
	*m = x;
	++m;
	*m = y;
	++m;
	*m = z;
	++m;
	*m = 0;
}
/*Sets the matrix to be a point by setting the 4th element to be 1.*/
void toPoint(Matrix *p){
	*(p->matrix+3) = 1.0;
}
/*Sets the matrix to be a vector by setting the 4th element to be 0.*/
void toVector(Matrix *v){
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
	newMatrix->numRows = maxRow;
	newMatrix->numCols = maxCol;
	/**Allocate the new matrix*/
	m = newMatrix->matrix;
	for(row = 0;row < maxRow; ++row){
		for(col = 0;col < maxCol; ++col){
			pos1 = m1->matrix + row * m1Cols; 
			pos2 = m2->matrix + col;
			sum = 0.0;
			for(p = 0;p < m1Cols;++p){
				sum += (*pos1) * (*pos2);
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
void placeMatrixCopy(Matrix *m1, Matrix *m2){
	int numRows, numCols;
	double *pos1, *pos2, *end;
	pos1 = m1->matrix;
	pos2 = m2->matrix;
	
	numRows = m1->numRows;
	numCols = m1->numCols;
	
	m2->numRows = numRows;
	m2->numCols = numCols;
	end = pos1 + numRows * numCols;
	while(pos1 < end){
		++pos1;
		++pos2;
	}
}
/*This transposes a matrix assuming it's */
char inPlaceTranspose(Matrix *m){
	int row, col;
	double *pos1, *pos2, swap;
	double *mat;
	int numRows = m->numRows;
	if(numRows != m->numCols)return 0;
	for(row = 0;row < numRows;++row){
		for(col = row + 1;col < numRows;++col){
			pos1 = mat + row * numRows + col;
			pos2 = mat + col * numRows + row;
			swap = *pos1;
			*pos1 = *pos2;
			*pos2 = swap;
		}
	}	
	return 1;
}
void placeScalarMultipleMatrix(Matrix *m,Matrix *m2, double scalar){
	int numRows, numCols;
	double *pos1, *pos2, *end;
	
	numRows = m->numRows;
	numCols = m->numCols;
	m2->numRows = numRows;
	m2->numCols = numCols;

	pos1 = m->matrix;
	pos2 = m2->matrix;
	end = m->matrix + numRows * numCols;
	while(pos1 < end){
		*pos2 = *pos1 * scalar;
		++pos1;
		++pos2;
	}
}
Matrix *getScalarMultipleMatrix(Matrix *m,double scalar){
	Matrix *productMatrix = malloc(sizeof(Matrix)); 
	productMatrix->matrix = malloc(sizeof(double) * m->numRows * m->numCols);
	placeScalarMultipleMatrix(m,productMatrix,scalar);
	return productMatrix;
}
Matrix *matrixCopy(Matrix *m){
	Matrix *copyMatrix = malloc(sizeof(Matrix));
	copyMatrix = malloc(sizeof(Matrix) * m->numRows * m->numCols);
	placeMatrixCopy(m,copyMatrix);
	return copyMatrix;
}
void printMatrix(Matrix *m){
	int i;
	for(i = 0;i < m->numRows * m->numCols;++i){
		printf("| %.2f |%s",m->matrix[i],((i + 1) % m->numCols) == 0 ? "\n" : " ");
	}
}
int main2(){
	Matrix m, m2, m3, m4, m5, *mptr;
	double list[16], list2[16], list3[16], list4[4], list5[4]; 
	m.matrix = list;
	m2.matrix = list2;
	m3.matrix = list3;
	m4.matrix = list4;
	m5.matrix = list5;
	
	placeTranslationMatrix(1,1,1,&m);
	printMatrix(&m);
	placeScaleMatrix(2,2,2,&m2);
	
	printf("\n");
	printMatrix(&m2);
	
	/*Test product matrix function*/
	placeProductMatrix(&m2,&m,&m3);
	
	printf("\n");
	printMatrix(&m3);

	mptr = getInverseMatrix(&m3);
	printf("\n");
	printMatrix(mptr);

	placePoint4(5,5,5, &m4);
	printf("\n");
	printMatrix(&m4);
	
	printf("\n");
	printMatrix(&m3);

	placeProductMatrix(&m3,&m4,&m5);
	printf("\n");
	printMatrix(&m5);
	return 0;
}
