#include <stdio.h>
#include <stdlib.h>
#include "matrix2.h"
/*Creating an improved version of the matrix library used by my raytracer.*/
/*This should hopefully help me make considerable improvements to the runtime.*/

/*Places a scaling matrix with the given scaling attributes into newMatrix*/
void placeScaleMatrix(double sX,double sY, double sZ, Matrix *newMatrix){
	double *m = newMatrix->matrix;
	double *mPos, *endM;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;

	/*Zero the matrix*/
	endM = m + 16;
	for(mPos = m;mPos < endM;++mPos){
		*mPos = 0.0;
	}
	/*Set the scaling attributes for the matrix.*/
	*m = sX;
	*(m+5) = sY;
	*(m+10) = sZ;
	*(m+15) = 1;
}
/*Places a translation matrix with the given translation attributes into newMatrix*/
void placeTranslationMatrix(double tX,double tY, double tZ, Matrix *newMatrix){
	double *m = newMatrix->matrix;
	double *mPos, *endM;
	newMatrix->numRows = 4;
	newMatrix->numCols = 4;
	
	/*Zero the matrix.*/
	endM = m + 16;
	for(mPos = m;mPos < endM;++mPos){
		*mPos = 0.0;
	}
	/*Set the matrix to an identity matrix, then add the translation attributes.*/
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
void placeInverseMatrix(Matrix *matrix, Matrix *inverseMatrix){
	Matrix m;
	double copyBuf[matrix->numRows * matrix->numCols];
	int rowCount = matrix->numRows;
	int row, col, otherRow;
	int pos;
	double swap;
	double val;
	double *rowPtr, *otherRowPtr, *rowEnd;
	double *IrowPtr, *IotherRowPtr;
	m.matrix = copyBuf;
	placeMatrixCopy(matrix,&m);

	placeScaleMatrix(1,1,1,inverseMatrix);
	/*Using gaussian elimination to solve a copy to the first argument matrix.
	 *Applying those same operations to an identity matrix to produce the inverse matrix.*/
	for(row = 0;row < rowCount;++row){
		col = row;
		while(row < rowCount){
			pos = rowCount * row + col;
			val = *(m.matrix + pos);
			if(val != 0.0)break;
			++row;
		}
		if(row == rowCount)return;

		val = 1.0 / val;

		/*Swap row and col, and scalar divide the row so that col,col = 1*/
		rowPtr = m.matrix + rowCount * col;
		otherRowPtr = m.matrix + rowCount * row;
		rowEnd = rowPtr + rowCount;


		IrowPtr = inverseMatrix->matrix + rowCount * col;
		IotherRowPtr = inverseMatrix->matrix + rowCount * row;
		
		while(rowPtr < rowEnd){
			swap = *rowPtr;
			*rowPtr = (*otherRowPtr) * val;
			*otherRowPtr = swap;

			swap = *IrowPtr;
			*IrowPtr = (*IotherRowPtr) * val;
			*IotherRowPtr = swap;

			++rowPtr; ++otherRowPtr; ++IrowPtr; ++IotherRowPtr;
		}
		row = col;

		/*Make sure every other row has a zero for this column.*/
		for(otherRow = 0;otherRow < col;++otherRow){
			val = *(m.matrix + otherRow * rowCount + col);
			rowPtr = m.matrix + rowCount * col;
			otherRowPtr = m.matrix + rowCount * otherRow;
		
			rowEnd = rowPtr + rowCount;
			IrowPtr = inverseMatrix->matrix + rowCount * col;
			IotherRowPtr = inverseMatrix->matrix + rowCount * otherRow;

			while(rowPtr < rowEnd){
				*otherRowPtr -= ((*rowPtr) * val);
				*IotherRowPtr -= ((*IrowPtr) * val);
				++rowPtr; ++otherRowPtr; ++IrowPtr; ++IotherRowPtr;
			}
		}

		for(otherRow = col+1;otherRow < rowCount;++otherRow){
			val = *(m.matrix + otherRow * rowCount + col);
			
			rowPtr = m.matrix + rowCount * col;
			otherRowPtr = m.matrix + rowCount * otherRow;

			IrowPtr = inverseMatrix->matrix + rowCount * col;
			IotherRowPtr = inverseMatrix->matrix + rowCount * otherRow;
			
			while(rowPtr < rowEnd){
				*otherRowPtr -= ((*rowPtr) * val);
				*IotherRowPtr -= ((*IrowPtr) * val);
				++rowPtr; ++otherRowPtr; ++IrowPtr; ++IotherRowPtr;
			}
				
		}
	}
}	

/**Allocates a scaling matrix with the given scaling attributes*/
Matrix *scaleMatrix(double sX, double sY, double sZ){
	Matrix *newMatrix = malloc(sizeof(Matrix));
	double *m = malloc(sizeof(double) * 16);
	newMatrix->matrix = m;
	placeScaleMatrix(sX,sY,sZ,newMatrix);
	return newMatrix;
}
/**Allocates a translation matrix with the given translation attributes*/
Matrix *translationMatrix(double tX,double tY,double tZ){
	Matrix *newMatrix = malloc(sizeof(Matrix));
	newMatrix->matrix = malloc(sizeof(double) * 16);
	placeTranslationMatrix(tX,tY,tZ,newMatrix);
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
	Matrix *newMatrix = malloc(sizeof(Matrix));
	double *m = malloc(sizeof(double) * 4);
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
	Matrix *newMatrix = malloc(sizeof(Matrix));
	double *m = malloc(sizeof(double) * 4);
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
		sum = (*pos1) + (*pos2);
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
		difference = (*pos1) - (*pos2);
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
	double *pos1, *pos2;
	int row, col, p;
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
	/*Set the elements of the new matrix*/
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
	/*Finalize the matrix structure that is returned.*/
	Matrix *newMatrix = malloc(sizeof(Matrix));
	newMatrix->matrix = malloc(sizeof(double) * 16);
	placeProductMatrix(m1,m2,newMatrix);
	return newMatrix;
}
/**Frees the heap resources used by the given matrix*/
void freeMatrix(Matrix *m){
	free(m->matrix);
	free(m);
}
/*Copies matrix 1 into matrix 2*/
void placeMatrixCopy(Matrix *m1, Matrix *m2){
	int numRows, numCols;
	double *pos1, *pos2, *end;
	
	numRows = m1->numRows;
	numCols = m1->numCols;
	m2->numRows = numRows;
	m2->numCols = numCols;
	
	pos1 = m1->matrix;
	pos2 = m2->matrix;
	end = pos1 + numRows * numCols;
	while(pos1 < end){
		*pos2 = *pos1;
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
	int numCols = m->numCols;
	if(numRows != numCols)return 0;
	mat = m->matrix;
	for(row = 0;row < numRows;++row){
		for(col = row;col < numRows;++col){
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
	end = pos1 + numRows * numCols;
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
	copyMatrix->matrix = malloc(sizeof(double) * m->numRows * m->numCols);
	placeMatrixCopy(m,copyMatrix);
	return copyMatrix;
}
void printMatrix(Matrix *m){
	int i;
	for(i = 0;i < m->numRows * m->numCols;++i){
		printf("| %.2f |%s",m->matrix[i],((i + 1) % m->numCols) == 0 ? "\n" : " ");
	}
	printf("\n");
}
int test(){
	//Matrix *m = translationMatrix(12,40,20);
	printf("Running test...\n");
	return 0;
}
int main(){
	Matrix m1, m2, m3, m4, *m5;
	Matrix *tmp;
	double buf[16], buf2[16], buf3[16], buf4[16];
	m1.matrix = buf;
	m2.matrix = buf2;
	m3.matrix = buf3;
	m4.matrix = buf4;
	placeTranslationMatrix(1,2,3,&m1);
	placeScaleMatrix(1,2,3,&m2);
	placeProductMatrix(&m1,&m2,&m3);
	placeProductMatrix(&m2,&m1,&m4);
	printMatrix(&m3);
	printMatrix(&m4);

	m5 = getProductMatrix(&m1,&m2);
	printf("m5:\n");
	printMatrix(m5);
	
	tmp = getInverseMatrix(m5);
	printf("inverse\n");
	printMatrix(tmp);

	printf("Product\n");
	tmp = getProductMatrix(tmp,m5);
	printMatrix(tmp);
	test();
	return 0;
}
