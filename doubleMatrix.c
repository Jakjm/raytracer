#include <stdlib.h>
#include <stdio.h>
#include "doubleMatrix.h"
//typedef struct Matrix{
//    double **matrix;
//    int rows;
//    int cols;
//} Matrix;
/**
 *Frees the memory used by this matrix so it can be reused.
 *@param Matrix *matrix- a pointer to the matrix that should be freed. 
 */
void freeMatrix(Matrix *matrix){
    int currentRow;
    for(currentRow = 0;currentRow < matrix->rows;++currentRow){
		free(matrix->matrix[currentRow]);
		matrix->matrix[currentRow] = NULL;
    }
    free(matrix->matrix);
	matrix->matrix = NULL;
    free(matrix);
}

//Creates a matrix with given dimensions with zeroed elements. 
Matrix *zeroMatrix(int rows,int cols){
	return createMatrix(rows,cols);
}
/** Creates a matrix in the dimensions specified. 
 * @param rows - the number of rows to place in the matrix.
 * @param cols - the number of elements within each row of the matrix.
 * @return - a pointer to the new matrix created.
 */
Matrix *createMatrix(int rows,int cols){
    /* Allocating space for the structure */
    Matrix *newMatrix = malloc(sizeof(Matrix));
    newMatrix->rows = rows;
    newMatrix->cols = cols;
    /* Now allocating space for the actual matrix. */
    newMatrix->matrix = malloc(sizeof(double*)*rows);

    int currentRow;

    /*Looping through each of the rows, and allocating a column of doubles. */
    for(currentRow = 0;currentRow < rows;++currentRow){
	newMatrix->matrix[currentRow] = calloc(cols,sizeof(double));
    }
    return newMatrix;
}
Matrix *getMatrix(double **m,int rows,int cols){
    Matrix *newMatrix = malloc(sizeof(Matrix));
    newMatrix->matrix = m;
    newMatrix->rows = rows;
    newMatrix->cols = cols;
    return newMatrix;
}
Matrix *identityMatrix(int);
Matrix *translationMatrix(double tX,double tY,double tZ){
    Matrix *m = identityMatrix(4);
    m->matrix[0][3] = tX;
    m->matrix[1][3] = tY;
    m->matrix[2][3] = tZ;
    return m;
}

//Produces a Scale affine-transformation matrix. 
Matrix *scaleMatrix(double sX,double sY,double sZ){
	Matrix *m = zeroMatrix(4,4);
	m->matrix[0][0] = sX;
	m->matrix[1][1] = sY;
	m->matrix[2][2] = sZ;
	m->matrix[3][3] = 1;
	return m;
}
//4 Dimensional vector creation method.
Matrix *vec4(double x,double y,double z){
	int r;
	Matrix *n = malloc(sizeof(Matrix));
	
	n->matrix = malloc(sizeof(double *) * 4);
	for(r = 0;r < 4;++r){
		n->matrix[r] = malloc(sizeof(double));
	}
	n->matrix[0][0] = x;
	n->matrix[1][0] = y;
	n->matrix[2][0] = z;
	n->matrix[3][0] = 0;
	n->rows = 4;
	n->cols = 1;
	return n;
}
//4 Dimensional point creation method. 
Matrix *point4(double x,double y,double z){
	Matrix *m = vec4(x,y,z);
	m->matrix[3][0] = 1;
	return m;
}
//Dot product calculation,
//Where both of the matricies passed are vectors. 
double dotProduct(Matrix *m1,Matrix *m2){
	double prod = 0;
	int i;
	//Be aware: we are skipping the bottom row since,
	//for our purposes, the fourth dimension should be ignored. 
	for(i = 0;i < m1->rows - 1;++i){
		prod += m1->matrix[i][0] * m2->matrix[i][0];
	}
	return prod;
}



/**
 *Creates a deep copy of the given matrix.
 *@param Matrix * - a pointer to the matrix that should be copied.
 *@return - a pointer to the copy of the parameter matrix.
 */
extern inline Matrix *matrixCopy(Matrix *matrix){
    double **newMatrix;

    /*Allocating memory to create row pointers for new matrix to point to. */
    newMatrix = malloc(sizeof(double*)*matrix->rows);

    /*Allocating some variables for the copying loops*/
    int currentRow;
    int currentCol;
    /* Looping through to allocate memory for each row pointer */
    for(currentRow = 0;currentRow < matrix->rows;++currentRow){
	newMatrix[currentRow] = malloc(sizeof(double) * matrix->cols);
	/* Now that the row pointer is allocated, we can copy the element of the original matrix onto the new one. */
	for(currentCol = 0; currentCol < matrix->cols;++currentCol){
	    newMatrix[currentRow][currentCol] = matrix->matrix[currentRow][currentCol];
	}
    }
    return getMatrix(newMatrix,matrix->rows,matrix->cols);
}

void printMatrix(Matrix *matrix){
	if(matrix == NULL){
		printf("Matrix is null!\n");
		return;
	}
    int row,col;
    printf("Matrix contents \n");
    for(row = 0; row < matrix->rows;row++){
	printf("| ");    
        for(col = 0; col < matrix->cols;col++){
           printf("% .3f  %c",matrix->matrix[row][col], (matrix->cols > matrix->rows && col == (matrix->rows - 1)) ? '|' : '\0');
        }
        printf("|\n");
    }
    printf("\n");
}

extern inline void swapRows(Matrix *matrix,int rowOne,int rowTwo){
     double *temp = matrix->matrix[rowOne];
     matrix->matrix[rowOne] = matrix->matrix[rowTwo];
     matrix->matrix[rowTwo] = temp; 
}

extern inline void scalarMultiplyRow(Matrix *matrix,int row,double scalar){
    int currentCol;
    for(currentCol = 0;currentCol < matrix->cols;++currentCol){
	matrix->matrix[row][currentCol] *= scalar;
    }
}
extern inline void addRow(Matrix *matrix,int sumRow,int addedRow,double scalar){
    int currentCol;
    for(currentCol = 0; currentCol < matrix->cols;++currentCol){
	matrix->matrix[sumRow][currentCol] += scalar * matrix->matrix[addedRow][currentCol];
    }
}
extern inline void subtractRow(Matrix *matrix,int differenceRow,int subtractedRow,double scalar){
    addRow(matrix,differenceRow,subtractedRow,-scalar);
}
//Returns the sum of two matricies. 
Matrix *getSumMatrix(Matrix *mOne,Matrix *mTwo){
	if(mOne->cols != mTwo->cols || mOne->rows != mTwo->rows)return NULL;
	int r,c;
	Matrix *sumMatrix = createMatrix(mOne->rows,mOne->cols);
	for(r = 0;r < mOne->rows;++r){
		for(c = 0;c < mOne->cols;++c){
			sumMatrix->matrix[r][c] = mOne->matrix[r][c] + mTwo->matrix[r][c];
		}
	}
	return sumMatrix;
}
/**
 *Generates the product of the given Matricies. Neither are modified.
 *@param Matrix *mOne - a pointer to the first matrix being multiplied.
 *@param Matrix *mTwo - a pointer to the second matrix being multiplied.
 *@return Matrix * - a pointer to the product matrix generated from the multiplied matricies. 
 */
Matrix *getProductMatrix(Matrix *mOne,Matrix*mTwo){
    /*If the num cols of the first matrix don't match the rows of the second matrix, then multiplication is undefined.*/
    if(mOne->cols != mTwo->rows)return NULL;

    /*Allocating memory for the product matrix */
    Matrix *productMatrix = createMatrix(mOne->rows,mTwo->cols);

    /*Allocating memory for keeping track of the row, column,multiplication position,and the sum of multiplication.*/
    int row,col,currentPosition;
    double sum;
    /*For each position inside the product matrix... */
    for(row = 0;row < productMatrix->rows;++row){
	for(col = 0;col < productMatrix->cols;++col){
	    /*Generating the sum by adding the products of the row values of matrix one by the col values of matrix two */
	    sum = 0;
            for(currentPosition = 0;currentPosition < mOne->cols;++currentPosition){
                sum += mOne->matrix[row][currentPosition] * mTwo->matrix[currentPosition][col];
	    }

	    /*Now that the sum has been found, setting the product matrix element to the sum.*/
	    productMatrix->matrix[row][col] = sum;
	}
    }
    return productMatrix;
}

/**
 *Produces a pointer to a newly formed identity matrix of the given size.
 *@param size - the number of rows and columns the Matrix should have.
 *@return Matrix * - a pointer to the newly formed identity matrix.
 */
extern inline Matrix *identityMatrix(int size){
    /**Creating the matrix */
    Matrix *matrix = createMatrix(size,size);
    int currentPos;
    /** Filling in the matrix with ones along the main diagonal */
    for(currentPos = 0;currentPos < matrix->rows;++currentPos){
	matrix->matrix[currentPos][currentPos] = 1;
    }
    return matrix;
}

/**
 *Produces a pointer to a scalar multiplied version of the given matrix, without modifying the original. 
 *@param Matrix *matrix - a pointer to the matrix that the new matrix is based off of. It is not modified.
 *@return Matrix * - a pointer to the new, scalar multiplied matrix. 
 */
Matrix *getScalarMultipleMatrix(Matrix *matrix,double scalar){
    Matrix *newMatrix = matrixCopy(matrix);
    int currentRow;
    for(currentRow = 0;currentRow < matrix->rows;currentRow++){
	scalarMultiplyRow(newMatrix,currentRow,scalar);
    }
    return newMatrix;
}
/**
 *Transposes the given matrix in place. The matrix is modified to reflect the change.
 *@param Matrix *matrix - a pointer to the matrix which should be transposed. 
 *@return 0 if the transpose was not possible, 1 if the transpose was successfully completed.
 */
char inPlaceTranspose(Matrix *matrix){
    if(matrix->rows != matrix->cols)return 0;
    int row,col;
    double temp;
    for(row = 0;row < matrix->rows - 1;++row){
	for(col = row + 1;col < matrix->cols;col++){
	    temp = matrix->matrix[row][col];
	    matrix->matrix[row][col] = matrix->matrix[col][row];
	    matrix->matrix[col][row] = temp;
	}
    }
    return 1;
}
//Performs a summation of two matricies in the place of m1.
char inPlaceSum(Matrix *m1,Matrix *m2){
	if(m1->rows != m2->rows || m1->cols != m2->cols)return 0;
	for(int i = 0;i < m1->rows;++i){
		for(int b = 0;b < m1->cols;++b){
			m1->matrix[i][b] += m2->matrix[i][b];
		}
	}
	return 1;
}
//Computes the difference of two matricies in the place of m1.
//In other words, m1 becomes m1 - m2; 
char inPlaceDifference(Matrix *m1,Matrix *m2){
	if(m1->rows != m2->rows || m1->cols != m2->cols)return 0;
	for(int i = 0;i < m1->rows;++i){
		for(int b = 0;b < m1->cols;++b){
			m1->matrix[i][b] -= m2->matrix[i][b];
		}
	}
	return 1;
}
void toVector(Matrix *m){
	m->matrix[3][0] = 0;
}
void inPlaceScalarMultiply(Matrix *m1,double scalar){
	for(int i = 0;i < m1->rows;++i){
		for(int b = 0;b < m1->cols;++b){
			m1->matrix[i][b] *= scalar;
		}
	}
}
/**
 * Computes the transpose matrix of the given matrix.
 * @param Matrix *matrix - a pointer to the matrix that the transpose should be made for. The structure is not modified.
 * @return Matrix * - a pointer to the newly formed transpose matrix.
 */
Matrix *getTransposeMatrix(Matrix *matrix){
    /* Allocating memory for the adjugate matrix, and variables for the row and col within matrix. */
    int currentRow,currentCol;
    Matrix *transpose = createMatrix(matrix->cols,matrix->rows);
    /* Looping through and copying the values from matrix to the new adjugate. */
    for(currentRow = 0;currentRow < matrix->rows;++currentRow){
	for(currentCol = 0;currentCol < matrix->cols;++currentCol){
	    transpose->matrix[currentCol][currentRow] = matrix->matrix[currentRow][currentCol];
	}
    }
    return transpose;
}
char inPlaceInverse(Matrix *m){
	if(m->rows != m->cols)return 0;
	Matrix *inv = getInverseMatrix(m);
	for(int r = 0;r < m->rows;++r){
		for(int c = 0;c < m->cols;++c){
			m->matrix[r][c] = inv->matrix[r][c];
		}
	}
	freeMatrix(inv);
	return 1;
}
/**
 *Computes an inverse matrix to the one given.
 *@param Matrix *matrix - a pointer to the matrix that the inverse should be made for.
 *@return Matrix * - a pointer to the newly formed inverse matrix, or NULL if the inverse DNE.
 */
extern inline Matrix *getInverseMatrix(Matrix *matrix){
    	/** If the matrix is not square, then inversion is impossible, so we return null. */
    	if(matrix->rows != matrix->cols)return NULL;

    	Matrix *workMatrix = createMatrix(matrix->rows,matrix->cols * 2);
    	int row, col;
    	for(row = 0;row < matrix->rows;++row){
		//Copying matrix into left side of working matrix
		for(col = 0; col < matrix->cols; ++col){
			workMatrix->matrix[row][col] = matrix->matrix[row][col];
		}
		//Putting identity matrix on the other side. 
		for(col = matrix->cols; col < workMatrix->cols; ++col){
			workMatrix->matrix[row][col] = ((row + matrix->cols) == col);
		}
	}
	solveMatrix(workMatrix);

    	Matrix *inverseMatrix = createMatrix(matrix->rows,matrix->cols);
	for(row = 0;row < inverseMatrix->rows;++row){
		for(col = 0;col < inverseMatrix->cols;++col){
			inverseMatrix->matrix[row][col] = workMatrix->matrix[row][col + inverseMatrix->cols];
		}
	}
    	/*Since we no longer need the working matrix, it is freed. */
    	freeMatrix(workMatrix);
    	return inverseMatrix;
}
void solveMatrix(Matrix *m){
	int currentRow = 0;
	int currentCol = 0;
	int remainingRow;
	double scalar;
	for(currentRow = 0;currentRow < m->rows;++currentRow){
		currentCol = currentRow;
		while(currentRow < m->rows){
			if(m->matrix[currentRow][currentCol] != 0){
				break;
			}
			++currentRow;
		}
		if(currentRow >= m->rows)return; 

		swapRows(m, currentRow, currentCol);
		currentRow = currentCol;
		scalar = 1.0 / m->matrix[currentRow][currentCol];
		scalarMultiplyRow(m,currentRow,scalar);

		for(remainingRow = currentRow + 1;remainingRow < m->rows;++remainingRow){
			scalar = m->matrix[remainingRow][currentCol] / m->matrix[currentRow][currentCol];
			subtractRow(m, remainingRow, currentRow, scalar);	
		
		}
		for(remainingRow = currentRow - 1;remainingRow >= 0;--remainingRow){
			scalar = m->matrix[remainingRow][currentCol];
			subtractRow(m,remainingRow,currentRow,scalar);
		}
	}
}
