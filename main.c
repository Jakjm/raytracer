#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/**Updating to allow the use of multithreading.*/
#include <pthread.h>
#include <time.h>
#include "matrix2.h"

//@Author Jordan Malek

#define MINIMUM_T 0.0000000001
#define NUM_BOUNCES 3
//Warning the compiler that I'll be defining these structs at some point.
typedef struct sphere sphere;
typedef struct light light;
typedef struct cube cube;

//Program variables....
//Variable that says whether debug info should be printed or not.
int verbose = 0;

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

	Matrix *matrix;
	Matrix *inverseMatrix;
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

	Matrix *matrix;
	Matrix *inverseMatrix;
	Matrix *inverseTranspose; 
} sphere;
Matrix *getSphereMatrix(sphere *);
//Function for parsing the sphere from the text file.
sphere *readSphere(FILE *fp){
	int result;
	sphere *oval = malloc(sizeof(sphere));
	oval->name = malloc(sizeof(char) * 30);

	//Reading the name, position and scale. Breaking if there is a read error.
	result = fscanf(fp, " %29s %lf %lf %lf %lf %lf %lf ", oval->name, &oval->posX, &oval->posY, &oval->posZ, &oval->scaleX, &oval->scaleY, &oval->scaleZ);
	if(result != 7)return NULL;
	
	//Reading the color and the lighting parameters. Breaking if there is a read error. 
	result = fscanf(fp, " %lf %lf %lf %lf %lf %lf %lf %d ",&oval->r,&oval->g,&oval->b,&oval->kAmb,&oval->kDif,&oval->kSpec,&oval->kR, &oval->specExp);
	if(result != 8)return NULL;
	
	oval->matrix = getSphereMatrix(oval);
	oval->inverseMatrix = getInverseMatrix(oval->matrix);
	oval->inverseTranspose = matrixCopy(oval->inverseMatrix);
	inPlaceTranspose(oval->inverseTranspose);
	
	return oval;
}
/*
 *Frees the heap space used by the given sphere. 
 */
void freeSphere(sphere *s){
	free(s->name);
	freeMatrix(s->matrix);
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
	cube *c = malloc(sizeof(cube));
	result = fscanf(fp," %lf %lf %lf %lf %lf %lf ",&c->posX,&c->posY,&c->posZ,&c->scaleX,&c->scaleY,&c->scaleZ);
	if(result != 6){
		free(c);
		return NULL;
	}
	c->matrix = getCubeMatrix(c);
	c->inverseMatrix = getInverseMatrix(c->matrix);
	return c;
}
void printCube(cube *c){
	printf("Cube Position(%.1lf %.1lf %.1lf)\n",c->posX,c->posY,c->posZ);
}
void freeCube(cube *c){
	freeMatrix(c->matrix);
	freeMatrix(c->inverseMatrix);
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
int parseFile(char *fileName){
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
	char *input = malloc(sizeof(char) * 30);
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
		//Allow the use of comments file.
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
	free(input);
	trimLists();
	if(readNear == 0 || readTop == 0 || readBottom == 0 || readLeft == 0 || readRight == 0){
		return -1;
	}
	if(readRes == 0 || readAmbient == 0 || readBack == 0 || readOutput == 0){
		return -1;
	}
	result = fclose(file);

	endTime = clock();
	if(verbose)printf("Parsing: clock time: %ld clocks per second: %ld runtime: %.3lfms\n",endTime - startTime,CLOCKS_PER_SEC,(endTime - startTime) / (CLOCKS_PER_SEC / 1000.0));
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
void save_image(int Width, int Height, char* fname,unsigned char* pixels) {
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
	if(verbose)printf("Writing: clock time: %ld clocks per sec: %ld runtime: %.3lfms\n",endTime - startTime,CLOCKS_PER_SEC,(endTime - startTime) / (CLOCKS_PER_SEC / 1000.0));
}
Matrix *getCubeMatrix(cube *c){
	Matrix scale; double scaleBuffer[16]; scale.matrix = scaleBuffer;
	Matrix translation; double translationBuffer[16]; translation.matrix = translationBuffer;
	placeTranslationMatrix(c->posX,c->posY,c->posZ,&translation);
	placeScaleMatrix(c->scaleX,c->scaleY,c->scaleZ,&scale);
	return getProductMatrix(&translation,&scale);
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

double computeTToCube(Matrix*ray,Matrix *origin,cube *s,double minimum);
/**Checks if there is a sphere in the way of this shadow ray on its way to a light source.
 * Used for checking if this spot should get some extra lighting.
 * */
int existsCollision(Matrix *origin,Matrix *ray){
	sphere *s;
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
	return 0;
}

//Computing the light color at the given point of collision.
//the colPoint is the point of collision with the surface. 
//the normal is the normal vector of collision with the surface
//s is the sphere. 
//TODO: this method can be easily adapted to both spheres and cubes by just inputting the lighting parameters instead of the kind of matrix.
//origin is the observing point - but this changes as rays are traced.
//returns the sum of all the light colors. 
void computeLightColor(Matrix *colPoint,Matrix *origin,Matrix *normal,sphere *s, double *red,double *green,double *blue){
	double cR = 0.0,cG = 0.0,cB = 0.0;
	int i;
	light *l;

	//The ray from the collision point to the light source.
	Matrix ray;
	double rayBuf[4]; ray.matrix = rayBuf;
	for(i = 0;i < numLights;++i){
		l = lightList[i];

		//Computing a vector from the surface of the sphere to the light. 
		placeMatrixCopy(l->lightPoint,&ray);
		inPlaceDifference(&ray,colPoint);

		//Now, need to ensure there aren't any spheres in the way.
		//In other words, need to check if there is a sphere in between the colPoint and the lightPoint
		if(!existsCollision(colPoint,&ray)){
			//Need to compute light, as long as there was a true collision.
			//Getting the dot product of the normal and the light ray. 
			double dot = dotProduct(&ray,normal);
			
			if(dot >= 0){
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
				dot /= sqrt(dotProduct(&ray,&ray));
				dot /= sqrt(dotProduct(normal,normal));

				//Computing the diffuse light. 
				diff = dot * s->kDif;
				difR = diff * l->iR * s->r;
			       	difG = diff * l->iG * s->g;
				difB = diff * l->iB * s->b;

				//Computing the reflection of the light ray on the surface. 
				placeScalarMultipleMatrix(&ray,&ref,-1);
				dot = 2*dotProduct(&ray,normal) / dotProduct(normal,normal);
				placeScalarMultipleMatrix(normal,&projection,dot);
				inPlaceSum(&ref,&projection);

				//Computing the viewing ray... 
				//The viewing ray should be the vector from the collision point to the eye. 
				placeScalarMultipleMatrix(colPoint,&viewRay,-1);
				inPlaceSum(&viewRay,origin);
				toVector(&viewRay);

				//Computing the shininess coefficient. 
				shininess = dotProduct(&ref,&viewRay);
				//If the shininess is zero, the reflected light ray and the view ray do not coincide. 
				//So there is no specular light. 
				if(shininess >= 0 && s->kSpec > 0.0){
				//shininess = fabs(shininess);
					shininess /= sqrt(dotProduct(&ref,&ref));
					shininess /= sqrt(dotProduct(&viewRay,&viewRay));
					shininess = pow(shininess,s->specExp);
			
				//Computing the specular light.
					spec = shininess * s->kSpec;
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
	}
	*red = cR;
	*green = cG;
	*blue = cB;
}
//Traces the ray to the closest sphere, if possible, then computes the color.
void traceRay(Matrix *ray,Matrix *origin,int bounceCount,double *red,double *green,double *blue){
	double cR, cG, cB;
	sphere *s = NULL;
	double t;
	double lowestT = 2020202020202020;
	int i;
	cube *c = NULL;
	//Computing the closest sphere intersection with the ray. 
	double minimum = MINIMUM_T;
	if(bounceCount == NUM_BOUNCES){
		minimum = 1 + MINIMUM_T;
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
		t = computeTToCube(ray,origin,cubeList[i],minimum);
		if(t > minimum && t < lowestT){
			lowestT = t;
			s = NULL;
			c = cubeList[i];
		}
	}
	t = lowestT;
	//There's a collision if we get inside this statement....
	//We'll start making recursive calls...
	if(s != NULL){
		//Variables for making the color calculations of the pixel.
		Matrix normal;double normalBuf[4]; normal.matrix = normalBuf;
		Matrix colPoint; double colPointBuf[4]; colPoint.matrix = colPointBuf;
		Matrix rayPrime; double rayPrimeBuf[4]; rayPrime.matrix = rayPrimeBuf;
		Matrix originPrime; double originPrimeBuf[4]; originPrime.matrix = originPrimeBuf;
		
		double lightR, lightG, lightB;
		
		placeScalarMultipleMatrix(ray,&colPoint,t);
		inPlaceSum(&colPoint,origin);
		//Computing the ambient light.
		cR = s->kAmb * aR * s->r;
		cG = s->kAmb * aG * s->g;
		cB = s->kAmb * aB * s->b;	
		
		//Building the normal vector...
		placeProductMatrix(s->inverseMatrix,ray,&rayPrime);
		placeProductMatrix(s->inverseMatrix,origin,&originPrime);
	
		//Getting the sum of the ray prime and the origin prime, then using the inverse transpose to generate the actual normal vector. 
		//Getting the collision point with respect to the canonical sphere - origin plus t times the ray. 
		inPlaceScalarMultiply(&rayPrime,t);
		inPlaceSum(&rayPrime,&originPrime);
		toVector(&rayPrime); //We want the vector with respect to the origin - but the origin is 0,0,0 so we can just knock off the point value. 
		//The normal is inverse transpose applied to the collision point minus (0,0,0) as a vector - but subtracing (0,0,0) is redundant so we skip that.  
		placeProductMatrix(s->inverseTranspose,&rayPrime,&normal);
		toVector(&normal);
		//If the ray from the origin to collision point is longer than the vector from the origin to the center of the sphere,
		//The normal should be flipped. 
		inPlaceDifference(&rayPrime,&originPrime);
		toVector(&rayPrime);  
		if(dotProduct(&rayPrime,&rayPrime) > dotProduct(&originPrime,&originPrime)){
			inPlaceScalarMultiply(&normal,-1);
		}

		//Light collision methods go here.
		computeLightColor(&colPoint,origin,&normal,s,&lightR,&lightG,&lightB);
		cR += lightR;
		cG += lightG;
		cB += lightB;

		//If the ray should be reflected, reflect it. Here goes nothing.
		if(bounceCount > 0 && s->kR > 0.0){
			double refR, refG, refB;

			Matrix reflectedRay;
			double reflectedRayBuffer[4];reflectedRay.matrix = reflectedRayBuffer;
			Matrix projection;
			double projectionBuffer[4];projection.matrix = projectionBuffer;
			--bounceCount;

			//calculating the projection of the ray onto the normal. 
			placeScalarMultipleMatrix(&normal,&projection,2 * ((dotProduct(ray,&normal) / dotProduct(&normal,&normal))));

			//Calculating the reflected ray. 
			placeMatrixCopy(ray,&reflectedRay);
			//Subtract 2 times the normal from the reflected ray.
			inPlaceDifference(&reflectedRay,&projection);

			traceRay(&reflectedRay,&colPoint,bounceCount,&refR,&refG,&refB);
			cR += (s->kR * refR);
			cG += (s->kR * refG);
			cB += (s->kR * refB);
		}
	}
	/*Do stuff for cubes...*/
	else if(c != NULL){
		cR = 0.0;
		cG = 0.0;
		cB = 0.0;
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
//TODO.... work in progress....
double computeTToCube(Matrix *ray,Matrix *origin,cube *c,double minimum){
	/*Allocate matrix for placing the product of m and ray, and m and origin.*/
	Matrix rayCP, originCP, *normal, surface;
	Matrix colPoint;
	double colPtBuf[4];colPoint.matrix = colPtBuf;
	double rayBuf[4], originBuf[4], surBuf[4];
	double originProj, rayProj, surProj;
	double minT, t;
	rayCP.matrix = rayBuf;
	originCP.matrix = originBuf;
	surface.matrix = surBuf;
	int i;
	minT = -1;
	
	placeProductMatrix(c->inverseMatrix,origin,&originCP);
	placeProductMatrix(c->inverseMatrix,ray,&rayCP);
	for(i = 0;i < numCubeMatricies;++i){
		normal = cubeMatricies[i];
		placeMatrixCopy(normal,&surface);
		toPoint(&surface);

		originProj = dotProduct(&originCP,normal);
		rayProj = dotProduct(&rayCP,normal);
		surProj = dotProduct(&surface,normal);
		t = (surProj - originProj) / rayProj;
		if(t > minimum){
			int a, b;
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
					if(minT < t){
						minT = t;
					}
				}
				else{
					minT = t;
				}
			}
		}
	}
	/*printMatrix(&colPoint)*/
	return minT;
}
//Traces the given ray to the given sphere,
//from the given starting point. 
double computeTToSphere(Matrix *ray,Matrix *origin,sphere *s,double minimum){
	double a,b,c;
	double t = -1;
	double det;

	/*Allocate matrix for placing the product of m and ray, and m and origin.*/
	Matrix rayCP, originCP;
	double rayBuf[4], originBuf[4];
	rayCP.matrix = rayBuf;
	originCP.matrix = originBuf;

	//Need to find the distance to the sphere, if a collision between the ray and the sphere exists. 
	//Obtains the matrix of the sphere.
	
	//Applying the matrix to the ray. 
	placeProductMatrix(s->inverseMatrix,ray,&rayCP);
	//Applying the matrix to the origin of the vector.
	placeProductMatrix(s->inverseMatrix,origin,&originCP);
	
	a = dotProduct(&rayCP,&rayCP);
	b = dotProduct(&originCP,&rayCP);
	c = dotProduct(&originCP,&originCP) - 1;

	det = b * b - a * c;
	//If there is a collision between the sphere and the ray...
	if(det >= 0.0){
		//Compute the earlier collision with t > 0.
		double rootDet = sqrt(det);
		double tOne;
		double tTwo;
		tOne = (-b - rootDet) / a;
		tTwo = (-b + rootDet) / a;

		//Computing which one should be the t. 
		
		if(tOne > minimum || tTwo > minimum){
			if(tOne <= tTwo){
				if(tOne >= minimum){
					t = tOne;
				}
				else{
					t = tTwo;
				}
			}
			else{
				if(tTwo >= minimum){
					t = tTwo;
				}
				else{
					t = tOne;
				}
			}
		}
	}
	return t;
}
/**This method is used to handle the raytracing for a particular portion of the scene.*/
void *computePixelThread(void *encoding){
	long encodedRowParms = (long)encoding;
	long mask = 0XFFFFFFFF;
	int rowStart = encodedRowParms >> 32;
	int rowEnd = (encodedRowParms & mask);

	/*The location of pixel 0,0 in space.*/
	double zeroX = (-cols / 2.0) + 0.5;
	double zeroY = (-rows / 2.0) + 0.5;

	double planeX = (right - left) / cols;
	double planeY = (bottom - top) / rows;

	Matrix eye;
	double eyeBuf[4];eye.matrix = eyeBuf;
	Matrix ray;
	double rayBuf[4];ray.matrix = rayBuf;

	int x, y;
	double rayX, rayY, rayZ;
	double cR, cG, cB;

	unsigned char *buffer = byteBuffer + (rowStart * cols * 3);
	placePoint4(&eye,0,0,0);
	placeVec4(&ray,0,0,0);
	/**Render the pixels that have been assigned to this thread.*/
	for(y = rowStart; y < rowEnd; ++y){
		for(x = 0;x < cols; ++x){
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

			/*Put the colour for this pixel in the buffer.*/
			*buffer = (unsigned char)(cR * 255);
			++buffer;
			*buffer = (unsigned char)(cG * 255);
			++buffer;
			*buffer = (unsigned char)(cB * 255);
			++buffer;
		}
	}
	return NULL;
}
void computePixels2(int threadCount){
	clock_t endTime, startTime = clock();
	int thread;
	/*The i-th thread renders rowStart to rowEnd (not including the row-endth row).*/
	int rowStart, rowEnd;
	/*A particular thread will render rowHeight * cols pixels.*/ 
	int rowHeight = rows / threadCount;
	long encoding;
	
	/**Default to 1 if threadCount invalid*/
	if(threadCount < 1 || threadCount > 16)threadCount = 1;
	pthread_t threads[threadCount];
	
	for(thread = 0;thread < threadCount;++thread){
		rowStart = thread * rowHeight;
		rowEnd = rowStart + rowHeight;
		encoding = rowStart;
		encoding = encoding << 32;
		encoding += rowEnd;
		pthread_create(&threads[thread],NULL,computePixelThread,(void*)encoding);
	}
	/*Wait to join each thread.*/
	for(int thread = 0;thread < threadCount;++thread){
		pthread_join(threads[thread],NULL);
	}

	endTime = clock();
	if(verbose)printf("Raytracing: clock time: %ld clocks per second: %ld runtime: %.3lfms\n",endTime - startTime,CLOCKS_PER_SEC,(endTime - startTime) / (CLOCKS_PER_SEC / 1000.0));
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
			if(numThreads < 0 || numThreads > 16){
				fprintf(stderr,"Incorrect number of threads\n");
				return 1;
			}
		}
		else{
			filename = arg;
			/**Check if the argument is a correct filename.*/
			int result = parseFile(filename);
			if(result == -1){
				fprintf(stderr,"::Parsing the file failed.\n::Please ensure you have provided a valid txt file in the format specified by the assignment.\n");
				return 1;
			}
			if(verbose)printParsedFile(filename);
		}
		++argv;
	}
	makeCubeMatricies();
	createImageArray();
	computePixels2(numThreads);
	//Freeing the lists of lights and spheres. 
	freeLists();
	
	//Saving the image file:
	save_image(cols,rows,outputFile,byteBuffer);

	//Freeing the remaining used resources. 
	free(byteBuffer);
	free(outputFile);
	return 0;
}
