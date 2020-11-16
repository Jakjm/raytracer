#ifndef DOUBLE_MATRIX
#define DOUBLE_MATRIX 1
typedef struct Matrix Matrix;
void swapRows(Matrix *,int,int);
void addRow(Matrix *,int ,int ,double );
void subtractRow(Matrix *matrix,int differenceRow,int subtractedRow,double scalar);
void scalarMultiplyRow(Matrix *matrix,int row,double scalar);
Matrix *createMatrix(int rows, int cols);
void solveMatrix(Matrix *m);
char inPlaceTranspose(Matrix *m);
char inPlaceSum(Matrix *,Matrix *);
char inPlaceDifference(Matrix *,Matrix*);
char inPlaceInverse(Matrix *);
void freeMatrix(Matrix *m);
void toVector(Matrix *m);
void inPlaceScalarMultiply(Matrix *,double);
void printMatrix(Matrix *);
//Matrix operations, sum, product, inverse
char placeProductMatrix(Matrix *,Matrix *,Matrix *);
Matrix *getProductMatrix(Matrix *m1,Matrix *m2);
void placeInverseMatrix(Matrix *, Matrix *);
Matrix *getInverseMatrix(Matrix *m1);
Matrix *getSumMatrix(Matrix *m1,Matrix *m2);
//Matrix generation.
Matrix *matrixCopy(Matrix *);
Matrix *point4(double x,double y,double z);
Matrix *vec4(double x,double y,double z);
Matrix *identityMatrix(int size);
Matrix *getScalarMultipleMatrix(Matrix *,double);

/**Set a vector to the following dimensions*/
void setVec4(Matrix *,double, double, double);
void setPoint4(Matrix *,double, double, double);

//Affine transformation matricies. 
Matrix *scaleMatrix(double sX,double sY,double sZ);
Matrix *translationMatrix(double,double,double);
void placeScaleMatrix(double, double, double, Matrix *);
void placeTranslationMatrix(double,double,double,Matrix *);
//Vector functions
double dotProduct(Matrix *,Matrix*);
#endif

