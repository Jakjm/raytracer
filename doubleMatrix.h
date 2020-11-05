#ifndef DOUBLE_MATRIX
#define DOUBLE_MATRIX
typedef struct Matrix Matrix;
typedef struct Matrix{
    double **matrix;
    int rows;
    int cols;
} Matrix;
extern inline void swapRows(Matrix *matrix,int rowOne,int rowTwo);
extern inline void addRow(Matrix *matrix,int sumRow,int addedRow,double scalar);
extern inline void subtractRow(Matrix *matrix,int differenceRow,int subtractedRow,double scalar);
extern inline void scalarMultiplyRow(Matrix *matrix,int row,double scalar);
extern inline Matrix *createMatrix(int rows, int cols);
extern inline void solveMatrix(Matrix *m);
char inPlaceTranspose(Matrix *m);
char inPlaceSum(Matrix *,Matrix *);
char inPlaceDifference(Matrix *,Matrix*);
char inPlaceInverse(Matrix *);
void freeMatrix(Matrix *m);
void toVector(Matrix *m);
void inPlaceScalarMultiply(Matrix *,double);
void printMatrix(Matrix *m);
//Matrix operations, sum, product, inverse
extern inline Matrix *getProductMatrix(Matrix *m1,Matrix *m2);
extern inline Matrix *getInverseMatrix(Matrix *m1);
extern inline Matrix *getSumMatrix(Matrix *m1,Matrix *m2);
//Matrix generation.
extern inline Matrix *matrixCopy(Matrix *);
extern inline Matrix *point4(double x,double y,double z);
extern inline Matrix *vec4(double x,double y,double z);
extern inline Matrix *identityMatrix(int size);
extern inline Matrix *getScalarMultipleMatrix(Matrix *,double);

/**Set a vector to the following dimensions*/
void setVec4(Matrix *,double, double, double);
void setPoint4(Matrix *,double, double, double);

//Affine transformation matricies. 
extern inline Matrix *scaleMatrix(double sX,double sY,double sZ);
extern inline Matrix *translationMatrix(double,double,double);

//Vector functions
double dotProduct(Matrix *,Matrix*);
#endif

