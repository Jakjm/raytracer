#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/**Updating to allow the use of multithreading.*/
#include <pthread.h>
#include <time.h>
#include <math.h>
#include "matrix2.h"

//@Author Jordan Malek
#define MINIMUM_T 0.0000000001
#define NUM_BOUNCES 3
#define MAX_THREADS 16
//Warning the compiler that I'll be defining these structs at some point.
typedef struct sphere sphere;
typedef struct light light;
typedef struct cube cube;

//Program variables....
//Variable that says whether debug info should be printed or not.
int verbose = 0;

/*Antialias option*/
int antialias = 0;

/**Values of the plane variables to be read by the parser. **/
float near, left, right, bottom, top;

/**Resolution of the image**/
int cols, rows;

//Spheres in the scene
int numSpheres = 0;

//Lights in the scene
int numLights = 0;

//cubes in the scene.
int numCubes = 0;
light **lightList;
cube **cubeList;
sphere **sphereList;
//Back color
double r, g, b;

//Ambient light
double aR, aG, aB;
char *outputFile;
unsigned char *byteBuffer;

typedef struct cube {
	double posX, posY, posZ;
	double scaleX, scaleY, scaleZ;
	double rX, rY;
	double r, g, b;
	double kAmb, kDif, kSpec, kR;
	int specExp;
	Matrix *inverseMatrix;
	Matrix *inverseTranspose;
} cube;
//Defining a structure for spheres. 
typedef struct sphere{
	//The name of the sphere. 
	char *name;
	//The position of the sphere in space.
	double posX, posY, posZ;
	//The scaling of the sphere. 
	double scaleX, scaleY, scaleZ;
	//Colors of the sphere, between 0 and 1. 
	double r, g, b;
	//Coefficients of reflectiveness
	double kAmb, kDif, kSpec, kR;
	//Specular brightness exponent
       	int specExp; 	

	Matrix *inverseMatrix;
	Matrix *inverseTranspose; 
} sphere;
Matrix *getSphereMatrix(sphere *);
//Function for parsing the sphere from the text file.
sphere *readSphere(FILE *fp){
	int result;
	Matrix *sMatrix;
	sphere *oval = malloc(sizeof(sphere));
	oval->name = malloc(sizeof(char) * 30);

	//Reading the name, position and scale. Breaking if there is a read error.
	result = fscanf(fp, " %29s %lf %lf %lf %lf %lf %lf ", oval->name, &oval->posX, &oval->posY, &oval->posZ, &oval->scaleX, &oval->scaleY, &oval->scaleZ);
	if(result != 7){
		free(oval->name);
		free(oval);
		return NULL;
	}
	
	//Reading the color and the lighting parameters. Breaking if there is a read error. 
	result = fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %d ",&oval->r,&oval->g,&oval->b,&oval->kAmb,&oval->kDif,&oval->kSpec,&oval->kR, &oval->specExp);
	if(result != 8){
		free(oval->name);
		free(oval);
		return NULL;
	}
	
	sMatrix = getSphereMatrix(oval);
	oval->inverseMatrix = getInverseMatrix(sMatrix);
	oval->inverseTranspose = matrixCopy(oval->inverseMatrix);
	inPlaceTranspose(oval->inverseTranspose);
	free(sMatrix);
	return oval;
}
/*
 *Frees the heap space used by the given sphere. 
 */
void freeSphere(sphere *s){
	free(s->name);
	freeMatrix(s->inverseMatrix);
	freeMatrix(s->inverseTranspose);
	free(s);
}
//Function for printing the parsed sphere. 
void printSphere(sphere *oval){
	printf("Sphere Name: %s Position: (%.1lf %.1lf %.1lf) Scale: (%.1lf %.1lf %.1lf)\n", oval->name, oval->posX, oval->posY, oval->posZ, oval->scaleX, oval->scaleY, oval->scaleZ);
	printf("\t Color: (%.1lf %.1lf %.1lf) Lighting Coefficients: Ambient: %.2lf Diffuse: %.2lf Specular: %.2lf Reflection: %lf Brightness: %d\n",oval->r,oval->g,oval->b,oval->kAmb,oval->kDif,oval->kSpec,oval->kR, oval->specExp);
}

//Defining a structure for lights. 
typedef struct light{
	char *name;
	double posX, posY, posZ;
	double iR, iG, iB;
	Matrix *lightPoint;
}light;

//Function for parsing the light from the text file.
light *readLight(FILE *fp){
	int result;
	light *lamp = malloc(sizeof(light));
	lamp->name = malloc(sizeof(char) * 30);

	//Reading the name and the position of the light. Breaking if reading fails. 
	result = fscanf(fp," %29s %lf %lf %lf ",lamp->name,&lamp->posX,&lamp->posY,&lamp->posZ);
	if(result != 4){
		free(lamp->name);
		free(lamp);
		return NULL;
	}

	//Reading the intensity values of the light, and breaking if there is an error with the read. 
	result = fscanf(fp," %lf %lf %lf ",&lamp->iR,&lamp->iG,&lamp->iB);
	if(result != 3){
		free(lamp->name);
		free(lamp);
		return NULL;
	}

	lamp->lightPoint = point4(lamp->posX,lamp->posY,lamp->posZ);
	return lamp;
}
/*
 *Frees the heap space used by the given light.
 */
void freeLight(light *l){
	free(l->name);
	freeMatrix(l->lightPoint);
	free(l);
}

//Function for printing the parsed light. 
void printLight(light *lamp){
	printf("Light name: %s Position(%.1lf %.1lf %.1lf)\n",lamp->name,lamp->posX,lamp->posY,lamp->posZ);
	printf("\t Intensity red: %.1lf green: %.1lf blue :%.1lf\n",lamp->iR,lamp->iG,lamp->iB);
}
Matrix *getCubeMatrix(cube *);
cube *readCube(FILE *fp){
	int result;
	Matrix *matrix;
	cube *c = malloc(sizeof(cube));
	result = fscanf(fp," %lf %lf %lf %lf %lf %lf ",&c->posX,&c->posY,&c->posZ,&c->scaleX,&c->scaleY,&c->scaleZ);
	if(result != 6){
		free(c);
		return NULL;
	}
	result = fscanf(fp," %lf %lf %lf %lf %lf ",&c->rX,&c->rY,&c->r,&c->g,&c->b);
	if(result != 5){
		free(c);
		return NULL;
	}
	result = fscanf(fp," %lf %lf %lf %lf %d ",&c->kAmb,&c->kDif,&c->kSpec,&c->kR,&c->specExp);
	if(result != 5){
		free(c);
		return NULL;
	}
	matrix = getCubeMatrix(c);
	c->inverseMatrix = getInverseMatrix(matrix);
	c->inverseTranspose = matrixCopy(c->inverseMatrix);
	inPlaceTranspose(c->inverseTranspose);
	free(matrix);
	return c;
}
void printCube(cube *c){
	printf("Cube Position(%.1lf %.1lf %.1lf)\n",c->posX,c->posY,c->posZ);
}
void freeCube(cube *c){
	freeMatrix(c->inverseMatrix);
	freeMatrix(c->inverseTranspose);
	free(c);
}
//Void for producing smaller versions of the lists, if possible.
void trimLists(){
	int i;
	if(numSpheres < 15){
		sphere **dummy = malloc(sizeof(sphere *) * numSpheres);
		for(i = 0;i < numSpheres;++i){
			dummy[i] = sphereList[i];
		}
		free(sphereList);
		sphereList = dummy;
	}
	if(numLights < 10){
		light **dummy2 = malloc(sizeof(light *) * numLights);
		for(i = 0;i < numLights;++i){
			dummy2[i] = lightList[i];
		}
		free(lightList);
		lightList = dummy2;
	}
	if(numCubes < 15){
		cube **dummy3 = malloc(sizeof(cube*) * numCubes);
		for(i = 0;i < numCubes;++i){
			dummy3[i] = cubeList[i];
		}
		free(cubeList);
		cubeList = dummy3;
	}
}
//Printing the results of the parsing of the file. 
void printParsedFile(char* fileName){
	int i;
	printf("~~~~~ Printing processed inputs ~~~~~\n\n");
	printf("Input file: %s\n\n",fileName); //Filename
	printf("Viewing Planes: Near %.2f Left %.2f Right %.2f Bottom %.2f Top %.2f\n",near,left,right,bottom,top); //Viewing planes
	printf("Picture Resolution (Width x Height): %dx%d \n",cols, rows); //Resolution
	//Printing the spheres that have been read. 
	if(numSpheres > 0){
		printf("\n");
		for(i = 0;i < numSpheres;++i){
			printSphere(sphereList[i]);
		}
		printf("\n");
	}
	else printf("\n");
	//Printing the lights that have been read. 
	if(numLights > 0){
		printf("\n");
		for(i = 0;i < numLights;++i){
			printLight(lightList[i]);
		}
		printf("\n");
	}
	if(numCubes > 0){
		printf("\n");
		for(i = 0;i < numCubes;++i){
			printCube(cubeList[i]);
		}
		printf("\n");
	}
	printf("Back Color: (%.1lf %.1lf %.1lf)\n",r,g,b); //Color of the background
	printf("Ambient Light: red: %.1lf green: %.1lf blue: %.1lf\n",aR,aG,aB); //Ambient light
	printf("\nOutput File: %s\n\n",outputFile); //Name of the output file.
	printf("~~~~~ Completed reading file ~~~~~\n");
}
//Second parser, if only to comply with the 'testParser' test case.
//Reads the necessary inputs from the file specified by the string.
//The inputs are stored in global variables. 

//Returns -1 if the parse was unsucessful.
int parseFile(char *fileName,long *parseTime){
	/*Measuring time used to read the file.*/
	clock_t endTime, startTime = clock();
	//Flags for whether the viewing planes have been parsed or not.
	int readNear, readLeft,readRight,readTop,readBottom;
	readNear = 0;readLeft = 0;readRight = 0;readTop = 0;readBottom = 0;
	
	//Flags for whether the other necessary inputs have been parsed or not.
	int readRes, readBack, readAmbient, readOutput;
	readRes = 0;readBack = 0;readAmbient = 0;readOutput = 0;


	//Initializing lists
	sphereList = malloc(sizeof(sphere *) * 15);
	lightList = malloc(sizeof(light *) * 10);
	cubeList = malloc(sizeof(cube *) * 15);

	//Opening file, returning with an error if file does not open
	FILE *file = fopen(fileName,"r");
	if(file == NULL)return -1;
	
	int result;
	char input[30];
	result = fscanf(file," %29s ",input);
	while(result == 1){
		//Read the near plane. 
		if(strcmp(input,"NEAR") == 0 && readNear == 0){
			result = fscanf(file," %f ",&near);
			if(result != 1)return -1;
			readNear = 1;
		}
		//Read the left plane.
		else if(strcmp(input,"LEFT") == 0 && readLeft == 0){
			result = fscanf(file," %f ",&left);
			if(result != 1)return -1;
			readLeft = 1;
		}
		//Read the right plane.
		else if(strcmp(input,"RIGHT") == 0 && readRight == 0){
			result = fscanf(file," %f ",&right);
			if(result != 1)return -1;
			readRight = 1;
		}
		//Read the bottom plane
		else if(strcmp(input,"BOTTOM") == 0 && readBottom == 0){
			result = fscanf(file," %f ",&bottom);
			if(result != 1)return -1;
			readBottom = 1;
		}
		//Read the top plane
		else if(strcmp(input,"TOP") == 0 && readTop == 0){
			result = fscanf(file," %f ",&top);
			if(result != 1)return -1;
			readTop = 1;
		}
		//Read screen resolution
		else if(strcmp(input,"RES") == 0 && readRes == 0){
			result = fscanf(file," %d %d ",&cols,&rows);
			if(result != 2)return -1;
			readRes = 1;
		}
		//Read a sphere, up to 15 times.
		else if(strcmp(input,"SPHERE") == 0){
			if(numSpheres >= 15){
				fprintf(stderr,"::Too many spheres. \n::There is a limit of 15 based on the given specifications.\n\n");
				return -1;
			}
			sphere *s = readSphere(file);
			if(s == NULL)return -1;
			sphereList[numSpheres] = s;
			++numSpheres; 
		}
		else if(strcmp(input,"CUBE") == 0){
			if(numCubes >= 15){
				fprintf(stderr,"::Too many cubes. \n::There is a limit of 15 based on the given specifications.\n\n");
				return -1;
			}
			cube *c = readCube(file);
			if(c == NULL)return -1;
			cubeList[numCubes] = c;
			++numCubes;
		}
		//Read a light, up to 10 times. 
		else if(strcmp(input,"LIGHT") == 0){
			if(numLights >= 10){
				fprintf(stderr,"::Too many lights. \n::There is a limit of 10 lights for the given specifications.\n\n");
				return -1;
			}
			light *l = readLight(file);
			if(l == NULL)return -1;
			lightList[numLights] = l;
			++numLights;
		}
		//Read the background color
		else if(strcmp(input,"BACK") == 0 && readBack == 0){
			result = fscanf(file," %lf %lf %lf ",&r,&g,&b);
			if(result != 3)return -1;
			readBack = 1;
		}
		//Read ambient light
		else if(strcmp(input,"AMBIENT") == 0 && readAmbient == 0){
			result = fscanf(file," %lf %lf %lf ",&aR,&aG,&aB);
			if(result != 3)return -1;
			readAmbient = 1;
		}
		//Read output filename. 
		else if(strcmp(input,"OUTPUT") == 0 && readOutput == 0){
			outputFile = malloc(sizeof(char) * 50);
			result = fscanf(file," %49s ",outputFile);
			if(result != 1)return -1;
			readOutput = 1;
		}
		//Allow the use of comments in file, preceeded by #.
		else if(strchr(input,'#')){
			if(!strchr(input,'\n')){
				int tmp = fgetc(file);
				while(tmp != '\n' && tmp != EOF){
					tmp = fgetc(file);
				}
			}
		}
		//If one of these input descriptors was not read, the file format is not valid. Exit. 
		else{
			return -1;
		}
		result = fscanf(file," %29s ",input);
	}

	//Cleaning up the extra memory used within this function
	trimLists();
	if(readNear == 0 || readTop == 0 || readBottom == 0 || readLeft == 0 || readRight == 0){
		return -1;
	}
	if(readRes == 0 || readAmbient == 0 || readBack == 0 || readOutput == 0){
		return -1;
	}
	result = fclose(file);
	
	endTime = clock();
	*parseTime = (endTime - startTime);
	return 0;		
}
//Function to convert rgb floats into an rgb color. 
int convertIntColor(double red, double green, double blue){
	if(red > 1.0 || green > 1.0 || blue > 1.0)printf("COLOR OVERFLOW %lf %lf %lf\n",red,green,blue);
	int color = (int)(red * 255) << 16;
       	color += (int)(green * 255) << 8;
	color += (int)(blue * 255); 
	return color; 	
}

//Function to create the image array for the ray tracer. 
void createImageArray(){
	byteBuffer = malloc(sizeof(unsigned char) * (rows * cols * 3));
}
//Function to save the image
//Copied from the website. 
void save_image(int Width, int Height, char* fname,unsigned char* pixels,long *saveTime) {
	FILE *fp;
	const int maxVal=255;
	clock_t endTime, startTime = clock();
	if(verbose)printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"wb");
	if(!fp){
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P6\n%d %d\n%d\n",Width,Height,maxVal);
	fwrite(pixels,3,Width*Height,fp);
	fclose(fp);
	endTime = clock();
	*saveTime = (endTime - startTime);
}
Matrix *getCubeMatrix(cube *c){
	Matrix rX; double rotationXBuffer[16]; rX.matrix = rotationXBuffer;
	Matrix rY; double rotationYBuffer[16]; rY.matrix = rotationYBuffer;
	Matrix scale; double scaleBuffer[16]; scale.matrix = scaleBuffer;
	Matrix translation; double translationBuffer[16]; translation.matrix = translationBuffer;
	Matrix product; double productBuffer[16];product.matrix = productBuffer;
	Matrix product2; double product2Buffer[16];product2.matrix = product2Buffer;

	/**double theta*/
	placeRotation_X_Matrix(&rX,c->rX);
	placeRotation_Y_Matrix(&rY,c->rY);
	placeTranslationMatrix(c->posX,c->posY,c->posZ,&translation);
	placeScaleMatrix(c->scaleX,c->scaleY,c->scaleZ,&scale);
	/*Product = rotationY * rotationX.*/
	placeProductMatrix(&rY,&rX,&product);
	/*Product2 = rotationY * rotationX * scale */
	placeProductMatrix(&product,&scale,&product2);
	/*Product3 = translation * scale * rY * rX*/
	return getProductMatrix(&translation,&product2);
}
//Gets the transformation matrix for the current sphere. 
Matrix *getSphereMatrix(sphere *s){
	Matrix *scale = scaleMatrix(s->scaleX,s->scaleY,s->scaleZ);
	Matrix *t = translationMatrix(s->posX,s->posY,s->posZ);
	Matrix *product = getProductMatrix(t,scale);
	freeMatrix(scale);
	freeMatrix(t);
	return product;
}
double computeTToSphere(Matrix*ray,Matrix *origin,sphere *s,double minimum);

double computeTToCube(Matrix*ray,Matrix *origin,Matrix **normal,cube *s,double minimum);
/**Checks if there is a sphere or a cube in the way of this shadow ray on its way to a light source.
 * Used for checking if this spot should get some extra lighting.
 * */
int existsCollision(Matrix *origin,Matrix *ray){
	sphere *s;
	cube *c;
	double t;
	sphere **list = sphereList;
	sphere **listEnd = sphereList + numSpheres;
	while(list < listEnd){	
		s = *list;
		t = computeTToSphere(ray,origin,s,MINIMUM_T);
		//If there is indeed an n such that 0 <= t <= 1, there is a collision between a sphere and the shadow ray.
		if(t >= MINIMUM_T && t <= 1.0 - MINIMUM_T){
			return 1;
		}
		++list;
	}
	cube **list2 = cubeList;
	cube **list2End = cubeList + numCubes;
	while(list2 < list2End){
		c = *list2;
		t = computeTToCube(ray,origin,NULL,c,MINIMUM_T);
		if(t >= MINIMUM_T && t <= 1.0 - MINIMUM_T){
			return 1;
		}
		++list2;	
	}
	return 0;
}

//Computing the light color at the given point of collision.
//the colPoint is the point of collision with the surface. 
//origin is the observing point - but this changes as rays are traced.
//returns the sum of all the light colors. 
void computeLightColor(Matrix *colPoint,Matrix *origin,Matrix *normal,double r, double g, double b, double kDif, double kSpec, int specExp, double *red,double *green,double *blue){
	double cR = 0.0,cG = 0.0,cB = 0.0;
	int i;
	light *l;
	double dot;
	//The ray from the collision point to the light source.
	Matrix shadowRay;
	double rayBuf[4]; shadowRay.matrix = rayBuf;
	double normalLengthSQ = dotProduct(normal,normal);
	double normalLength = sqrt(normalLengthSQ);
	for(i = 0;i < numLights;++i){
		l = lightList[i];

		//Computing a vector from the surface of the sphere to the light. 
		placeMatrixCopy(l->lightPoint,&shadowRay);
		inPlaceDifference(&shadowRay,colPoint);

		/*Can skip this light if the dot product between the shadow ray and the normal vector is < 0*/
		dot = dotProduct(&shadowRay,normal);
		if(dot < 0.0)continue;
		//Now, need to ensure there aren't any spheres in the way.
		//In other words, need to check if there is a sphere in between the colPoint and the lightPoint
		if(!existsCollision(colPoint,&shadowRay)){
			//Need to compute light, as long as there was a true collision.
			//Getting the dot product of the normal and the light ray. 
			
			//Color placeholders. 
			double specR,specG,specB;
			double difR,difG,difB;
			//The "r" vector;the light ray after it has been reflected on the surface using the normal.
			Matrix ref;
			double refBuf[4]; ref.matrix = refBuf;
			//The "v" vector;the viewing ray from the eye to the collision point. 
			Matrix viewRay;
			double viewRayBuf[4];viewRay.matrix = viewRayBuf;
			/*Projectiion of the ray to calculate reflection*/
			Matrix projection;
			double projectionBuf[4];projection.matrix = projectionBuf;
				
			double shininess;
			double spec,diff;

			//Diminishing the dot product between the ray.
			dot /= sqrt(dotProduct(&shadowRay,&shadowRay));
			dot /= normalLength;

			//Computing the diffuse light. 
			diff = dot * kDif;
			difR = diff * l->iR * r;
		       	difG = diff * l->iG * g;
			difB = diff * l->iB * b;

			/*Calculating the reflection of the negative shadow ray off the surface*/
			/*Ref = negative shadow ray. In other words vector from light to colPoint*/
			placeScalarMultipleMatrix(&shadowRay,&ref,-1);

			/*Twice the projection of the ref onto normal gives us the amount of bounce.*/
			dot = 2*dotProduct(&ref,normal) / normalLengthSQ;
			placeScalarMultipleMatrix(normal,&projection,dot);
			
			/*Subtract the bounce from ref.*/
			inPlaceDifference(&ref,&projection);

			//Computing the viewing ray... 
			//The viewing ray should be the vector from the collision point to the eye. 
			placeScalarMultipleMatrix(colPoint,&viewRay,-1);
			inPlaceSum(&viewRay,origin);
			toVector(&viewRay);

			//Computing the shininess coefficient. 
			shininess = dotProduct(&ref,&viewRay);
			//If the shininess is zero, the reflected light ray and the view ray do not coincide. 
			//So there is no specular light. 
			if(shininess >= 0 && kSpec > 0.0){
				//shininess = fabs(shininess);
				shininess /= sqrt(dotProduct(&ref,&ref));
				shininess /= sqrt(dotProduct(&viewRay,&viewRay));
				shininess = pow(shininess,specExp);
		
				//Computing the specular light.
				spec = shininess * kSpec;
				specR = spec * l->iR;
				specG = spec * l->iG;
				specB = spec * l->iB;
			}
			else{
				specR = 0;
				specG = 0;
				specB = 0;
			}
			cR += (difR + specR);
			cG += (difG + specG);
			cB += (difB + specB);
		}
	}
	*red = cR;
	*green = cG;
	*blue = cB;
}
//Traces the ray to the closest sphere, if possible, then computes the color.
void traceRay(Matrix *ray,Matrix *origin,int bounceCount,double *red,double *green,double *blue){
	double cR, cG, cB;
	sphere *s = NULL;
	double t; double lowestT = INFINITY;
	int i;
	cube *c = NULL;
	Matrix *potentialNormal;
	Matrix *normalPrime = NULL;
	Matrix normal; double normalBuf[4]; normal.matrix = normalBuf;
	//Computing the closest sphere intersection with the ray. 
	double minimum = MINIMUM_T;
	if(bounceCount == NUM_BOUNCES){/*If bounceCount = NUM_BOUNCES, collision must be after near plane*/
		minimum = near + MINIMUM_T;
	}
	for(i = 0;i < numSpheres;++i){
		//Need to compute the normal here...
		t = computeTToSphere(ray,origin,sphereList[i],minimum);
		//If t is smaller than the lowest T so far, we take it. 
		if(t > minimum && t < lowestT){
			lowestT = t;
			s = sphereList[i];
		}
	}
	for(i = 0;i < numCubes;++i){
		t = computeTToCube(ray,origin,&potentialNormal,cubeList[i],minimum);
		if(t > minimum && t < lowestT){
			lowestT = t;
			s = NULL;
			c = cubeList[i];
			normalPrime = potentialNormal;
		}
	}
	/*Set t to the lowest t found that is still greater than the minimum*/
	t = lowestT;
	if(s != NULL || c != NULL){
		double kAmb;
		double shapeRed, shapeGreen, shapeBlue, kDif, kSpec, kRef;
		int specExp;
		double lightR, lightG, lightB;
		/*Calculate the collision point...*/
		Matrix colPoint; double colPointBuf[4]; colPoint.matrix = colPointBuf;
		
		/*Col point = origin + ray * t*/
		placeScalarMultipleMatrix(ray,&colPoint,t);
		inPlaceSum(&colPoint,origin);

		//There's a collision if we get inside this statement....
		//We'll start making recursive calls...
		if(s != NULL){
			Matrix rayPrime; double rayPrimeBuf[4]; rayPrime.matrix = rayPrimeBuf;
			Matrix originPrime; double originPrimeBuf[4]; originPrime.matrix = originPrimeBuf;
			Matrix colPointPrime; double colPointPrimeBuf[4]; colPointPrime.matrix = colPointPrimeBuf;
			
			//Variables for making the color calculations of the pixel.
			shapeRed = s->r;
			shapeGreen = s->g;
			shapeBlue = s->b;
			kDif = s->kDif;
			kSpec = s->kSpec;
			kAmb = s->kAmb;
			kRef = s->kR;
			specExp = s->specExp;
			
			//Building the normal vector...
			placeProductMatrix(s->inverseMatrix,ray,&rayPrime);
			placeProductMatrix(s->inverseMatrix,origin,&originPrime);

			placeScalarMultipleMatrix(&rayPrime,&colPointPrime,t);
			inPlaceSum(&colPointPrime,&originPrime);
			toVector(&colPointPrime); /*Vectorize the colPointPrime before applying the inverse transpose to it*/

			
			/*If the distance from the eye to the collision point is longer than the distance from the eye to the center of the sphere,
			 *then the collision happened at the back of the sphere, not the front*/		
			if(dotProduct(&rayPrime,&rayPrime) > dotProduct(&originPrime,&originPrime)){
				inPlaceScalarMultiply(&colPointPrime,-1);
			}
			normalPrime = &colPointPrime;
			/*Apply the inverse transpose of the sphere to get the canonical normal.*/
			placeProductMatrix(s->inverseTranspose,&colPointPrime,&normal);
			toVector(&normal);
		}
		/*Set lighting parameters and colour for cubes.*/
		else if(c != NULL){
			shapeRed = c->r;
			shapeGreen = c->g;
			shapeBlue = c->b;
			kAmb = c->kAmb;
			kRef = c->kR;
			kDif = c->kDif;
			kSpec = c->kSpec;
			specExp = c->specExp;

			/**Apply the inverse transpose of the cube to the normal with respect to the cube to get the true normal vector*/ 
			placeProductMatrix(c->inverseTranspose,normalPrime,&normal);
			toVector(&normal); /*Need to set the normal to a vector...*/
		}
		/*Calculating the ambient light*/
		cR = kAmb * aR * shapeRed;
		cG = kAmb * aG * shapeGreen;
		cB = kAmb * aB * shapeBlue;
		/*Calculating the diffuse and speciular light*/
		computeLightColor(&colPoint,origin,&normal,shapeRed,shapeGreen,shapeBlue,kDif,kSpec,specExp,&lightR,&lightG,&lightB);
		cR += lightR;
		cG += lightG;
		cB += lightB;
		
		/*Calculating the reflective light*/
		if(bounceCount > 0 && kRef > 0.0){
				double refR, refG, refB;
				
				Matrix reflectedRay;
				double reflectedRayBuffer[4];reflectedRay.matrix = reflectedRayBuffer;
				Matrix projection;
				double projectionBuffer[4];projection.matrix = projectionBuffer;
				--bounceCount;
				
				/*Calculating the projection of the ray onto the normal vector. Subtracting twice the normal projection to get the reflection.*/
				placeScalarMultipleMatrix(&normal,&projection,2 * ((dotProduct(ray,&normal) / dotProduct(&normal,&normal))));
				/*Copy ray into reflected ray*/
				placeMatrixCopy(ray,&reflectedRay);
				/*Subtract projection from reflected ray*/
				inPlaceDifference(&reflectedRay,&projection);
				
				traceRay(&reflectedRay,&colPoint,bounceCount,&refR,&refG,&refB);
				cR += (kRef * refR);
				cG += (kRef * refG);
				cB += (kRef * refB);
		}
	}
	else{
		//If this is a bounced ray, return black if there is no collision.
		if(bounceCount < NUM_BOUNCES){
			cR = 0.0;
			cG = 0.0;
			cB = 0.0;
		}
		//Otherwise, return the background color. 
		else{
			cR = r;
			cG = g;
			cB = b;
		}
	}
	*red = cR;
	*green = cG;
	*blue = cB;
}
Matrix **cubeMatricies;
int numCubeMatricies = 6;
void makeCubeMatricies(){
	Matrix **mList;
	cubeMatricies = malloc(sizeof(Matrix*) * numCubeMatricies);
	mList = cubeMatricies;
	*mList = vec4(-1.0,0.0,0.0);
	++mList;
	*mList = vec4(1.0,0.0,0.0);
	++mList;
	*mList = vec4(0.0,-1.0,0.0);
	++mList;
	*mList = vec4(0.0,1.0,0.0);
	++mList;
	*mList = vec4(0.0,0.0,-1.0);
	++mList;
	*mList = vec4(0.0,0.0,1.0);
}
double computeTToCube(Matrix *ray,Matrix *origin,Matrix **n, cube *c,double minimum){
	/*Allocate matrix for placing the product of m and ray, and m and origin.*/
	Matrix rayCP, originCP, surface, *normal;
	Matrix colPoint;
	double colPtBuf[4];colPoint.matrix = colPtBuf;
	double rayBuf[4], originBuf[4], surBuf[4];
	double originProj, rayProj, surProj, distance;
	double minT, t;
	rayCP.matrix = rayBuf;
	originCP.matrix = originBuf;
	surface.matrix = surBuf;
	int i;
	int a, b;
	minT = -1;
	
	placeProductMatrix(c->inverseMatrix,origin,&originCP);
	placeProductMatrix(c->inverseMatrix,ray,&rayCP);
	for(i = 0;i < numCubeMatricies;++i){
		normal = cubeMatricies[i];
		placeMatrixCopy(normal,&surface);
		toPoint(&surface);

		rayProj = dotProduct(&rayCP,normal);
		originProj = dotProduct(&originCP,normal);
		surProj = dotProduct(&surface,normal);
		distance = surProj - originProj;
		if(rayProj == 0.0){
			continue;
		}
		else if(rayProj > 0.0){
			//Normal needs to be flipped in this case, since ray and normal are in the same direction.
			if((i % 2) == 0){
				normal = cubeMatricies[i+1];
			}
			else{
				normal = cubeMatricies[i-1];
			}
		}
		//distance = originProj - surProj;
		t = distance / rayProj;
		if(t < minimum){
			continue;
		}

		placeScalarMultipleMatrix(&rayCP,&colPoint,t);
		inPlaceSum(&colPoint,&originCP);
		if(i < 2){
			a = 1;
			b = 2;
		}
		else if(i < 4){
			a = 0;
			b = 2;
		}
		else{
			a = 0;
			b = 1;
		}
		if((colPtBuf[a] >= -1.0 && colPtBuf[a] <= 1.0) && (colPtBuf[b] >= -1.0 && colPtBuf[b] <= 1.0)){
			if(minT > 0){
				if(t < minT){
					minT = t;
					if(n != NULL)*n = normal;
				}
			}
			else{
				minT = t;
				if(n != NULL)*n = normal;
			}
		}
	}
	return minT;
}
//Traces the given ray to the given sphere,
//from the given starting point. 
double computeTToSphere(Matrix *ray,Matrix *origin,sphere *s,double minimum){
	double a,b,c;
	double t = -1;
	double det;
	double originDistSQ;

	/*Allocate matrix for placing the product of m and ray, and m and origin.*/
	Matrix rayPrime, originPrime;
	double rayBuf[4], originBuf[4];
	rayPrime.matrix = rayBuf;
	originPrime.matrix = originBuf;

	//Need to find the distance to the sphere, if a collision between the ray and the sphere exists. 
	//Obtains the matrix of the sphere.
	
	//Applying the matrix to the ray. 
	placeProductMatrix(s->inverseMatrix,ray,&rayPrime);
	//Applying the matrix to the origin of the vector.
	placeProductMatrix(s->inverseMatrix,origin,&originPrime);
	
	a = dotProduct(&rayPrime,&rayPrime); /*Length of ray^2*/
	b = dotProduct(&originPrime,&rayPrime);/*Length of origin along ray*/
	originDistSQ = dotProduct(&originPrime,&originPrime); /*Length of origin wrt spehre ^ 2*/
	c = originDistSQ - 1; /*Length of origin with respect to sphere ^ 2 - 1*/

	det = b * b - a * c;
	//If there is a collision between the sphere and the ray...
	if(det >= 0.0){
		//Compute the earlier collision with t > 0.
		double rootDet = sqrt(det);
		double reciprocalA = 1.0 / a;
		double tOne;
		double tTwo;
		tOne = (-b - rootDet) * reciprocalA;
		tTwo = (-b + rootDet) * reciprocalA;
		if(tOne > minimum){
			t = tOne;
		}
		else if(tTwo > minimum){
			t = tTwo;
		}
	}
	return t;
}


void pack(int *arr, int start, int end){
	arr[0] = start;
	arr[1] = end;
}
void unpack(int *arr, int *start,int *end){
	*start = arr[0];
	*end = arr[1];
}
/**This method is used to handle the raytracing for a particular portion of the scene.*/
void *computePixelThread(void *range){
	//long mask = 0XFFFFFFFF;
	int rowStart, rowEnd;
	unpack(range,&rowStart,&rowEnd);
	/*The location of the pixel at row 0, column 0 in camera coordinates*/
	double zeroX = ((left * cols) * 0.5) + 0.5;
	double zeroY = ((bottom * rows) * 0.5) + 0.5;

	double planeX = (right - left) / cols;
	double planeY = -((top - bottom) / rows);

	Matrix eye;
	double eyeBuf[4];eye.matrix = eyeBuf;
	Matrix ray;
	double rayBuf[4];ray.matrix = rayBuf;

	int x, y;
	double rayX, rayY, rayZ;
	double tR, tG, tB;
	double cR, cG, cB;
	double antiX, antiY;
	double antiCoefficient = 1.0 / 9.0;
	unsigned char *buffer = byteBuffer + (rowStart * cols * 3);
	placePoint4(&eye,0,0,0);
	placeVec4(&ray,0,0,0);

	/**Render the pixels that have been assigned to this thread.*/
	for(y = rowStart; y < rowEnd; ++y){
		for(x = 0;x < cols; ++x){
			tR = 0.0; tG = 0.0; tB = 0.0;
			
			/*Calculate the ray for the particular pixel we're rendering.*/
			rayX = (x + zeroX) * planeX;
			rayY = (y + zeroY) * planeY;
			rayZ = -near;
			setVec4(&ray,rayX,rayY,rayZ);
			setPoint4(&eye,0.0,0.0,0.0);

			//Computing the pixel color.
			traceRay(&ray,&eye,NUM_BOUNCES,&cR,&cG,&cB);
			
			//Clamping the color, if the color has exceeded one. 
			if(cR > 1.0)cR = 1.0;
			if(cG > 1.0)cG = 1.0;
			if(cB > 1.0)cB = 1.0;
			
			tR += cR; 
			tG += cG;
			tB += cB;
			if(antialias){
				for(antiX = -0.5;antiX <= 0.52;antiX += 0.5){
					for(antiY = -0.5;antiY <= 0.52;antiY += 0.5){
						/*Skip the middle ray since it has already been considered.*/
						if(antiX >= -0.1 && antiX < 0.1 && antiY >= -0.1 && antiY <= 0.1)continue;
						rayX = (x + zeroX + antiX) * planeX;
						rayY = (y + zeroY + antiY) * planeY;
						rayZ = -near;

						setVec4(&ray,rayX,rayY,rayZ);
						setPoint4(&eye,0.0,0.0,0.0);
						traceRay(&ray,&eye,NUM_BOUNCES,&cR,&cG,&cB);
						//Clamping the color, if the color has exceeded one. 
						if(cR > 1.0)cR = 1.0;
						if(cG > 1.0)cG = 1.0;
						if(cB > 1.0)cB = 1.0;
						tR += cR; 
						tG += cG;
						tB += cB;
					}
				}
				/*Average out colours of samples.*/
				tR *= antiCoefficient;
				tG *= antiCoefficient;
				tB *= antiCoefficient;
			}
			
			
			

			/*Put the colour for this pixel in the buffer.*/
			*buffer = (unsigned char)(tR * 255);
			++buffer;
			*buffer = (unsigned char)(tG * 255);
			++buffer;
			*buffer = (unsigned char)(tB * 255);
			++buffer;
		}
	}
	return NULL;
}
/*This method distributes the workload to render the image along multiple threads.*/
void computePixels2(int threadCount,long *renderTime){
	clock_t endTime, startTime = clock();
	int thread;
	/*The i-th thread renders rowStart to rowEnd (not including the row-endth row).*/
	int rowStart, rowEnd;
	
	/*A particular thread will render rowHeight * cols pixels.*/ 
	int rowHeight = rows / threadCount;

	//End offset
	int offset = rows % threadCount;
	

	/**
	 * EG. if using 4 threads
	 * # # # # # < - done by main thread
	 * # # # # # < - done by helper thread 1
	 * # # # # # < - done by helper thread 2
	 * # # # # # < - done by helper thread 3 
	 */

	pthread_t threads[MAX_THREADS];
	int range[2];
	/*The work will be divided equally among the main thread (the currently active ones)
	 *and the helper threads.*/
	/*The helper threads will handle lower slices of rows of the image.*/
	for(thread = 1;thread < threadCount;++thread){
		rowStart = (thread * rowHeight) + offset;
		rowEnd = rowStart + rowHeight;
	 	pack(range,rowStart,rowEnd);
		
		//Create a thread to handle the workload
		pthread_create(&threads[thread - 1],NULL,computePixelThread,(void*)range);
	}

	/*The main thread will handle the first slice of rows*/
	rowStart = 0;
	rowEnd = rowHeight + offset;
	pack(range,rowStart,rowEnd);
	computePixelThread((void*)range);

	/*Wait to join each thread.*/
	for(int thread = 1;thread < threadCount;++thread){
		pthread_join(threads[thread - 1],NULL);
	}

	endTime = clock();
	*renderTime = (endTime - startTime); 
}
//Freeing the sphereList and the spheres it contains
//and the light list and the lights contained. 
void freeLists(){
	int i;
	for(i = 0;i < numSpheres;++i){
		freeSphere(sphereList[i]);
	}
	free(sphereList);
	for(i = 0;i < numLights;++i){
		freeLight(lightList[i]);
	}
	free(lightList);
	for(i = 0;i < numCubes;++i){
		freeCube(cubeList[i]);
	}
	free(cubeList);
	for(i = 0;i < numCubeMatricies;++i){
		freeMatrix(cubeMatricies[i]);
	}
	free(cubeMatricies);
}
//Main function of the program. 
int main(int argc,char **argv){
	char **endArgs = argv + argc;
	char *arg;

	char *filename = NULL;
	int numThreads = 1;
	long parseTime, renderTime, saveTime;
	//Checking if the program lacks arguments. 
	if(argc < 2){
		fprintf(stderr,"Please retry running the program with the correct number of arguments\n");
		return 1;
	}
	++argv;
	while(argv < endArgs){
		arg = *argv;
		/**If argument matches filename*/
		if(strcmp(arg,"-v") == 0){
			verbose = 1;
		}
		else if(sscanf(arg,"-t%d",&numThreads) == 1){
			/**Verify that the number of requested threads is valid...*/
			if(numThreads <= 0 || numThreads > MAX_THREADS){
				fprintf(stderr,"Incorrect number of threads\n");
				return 1;
			}
		}
		else if(strcmp(arg,"-aa") == 0){
			antialias = 1;
		}
		else if(filename == NULL){
			filename = arg;
			/**Check if the argument is a correct filename.*/
			int result = parseFile(filename,&parseTime);
			if(result == -1){
				fprintf(stderr,"Parsing the file failed.\n::Please ensure you have provided a valid txt file in the specified format.\n");
				return 1;
			}
		}
		++argv;
	}
	/*If no file was read, exit.*/
	if(filename == NULL){
		fprintf(stderr,"No input file specified. Please specify a text file describing a scene to render.\n");
		return 1;
	}
	if(verbose)printParsedFile(filename);
	makeCubeMatricies();
	createImageArray();
	computePixels2(numThreads,&renderTime);
	//Freeing the lists of lights and spheres. 
	freeLists();
	
	//Saving the image file:
	save_image(cols,rows,outputFile,byteBuffer,&saveTime);

	/*Print runtimes for different segments of raytracer, if the option was specified.*/
	if(verbose){
		printf("Summary of runtime\n");
		printf("Parsing: %.3lf ms\n", parseTime * (1000.0 / CLOCKS_PER_SEC));	
		printf("Rendering: %.3lf ms\n",renderTime * (1000.0 / CLOCKS_PER_SEC));
		printf("Save time: %.3lf ms\n",saveTime * (1000.0 / CLOCKS_PER_SEC));
	}
	
	//Freeing the remaining used resources. 
	free(byteBuffer);
	free(outputFile);
	return 0;
}
