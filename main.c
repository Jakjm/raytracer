#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "doubleMatrix.h"
#include "doubleMatrix.c"

//@Author Jordan Malek

#define MINIMUM_T 0.0000000001
#define NUM_BOUNCES 3
//Warning the compiler that I'll be defining these structs at some point.
typedef struct sphere sphere;
typedef struct light light;

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
light **lightList;
sphere **sphereList;
//Back color
double r, g, b;

//Ambient light
double aR, aG, aB;
char *outputFile;
int **imageArray;
char *byteBuffer;

//Defining a structure for color. 
typedef struct color{
	double r;
	double g;
	double b;
} color;

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
} sphere;

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

	return oval;
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
}light;

//Function for parsing the light from the text file.
light *readLight(FILE *fp){
	int result;
	light *lamp = malloc(sizeof(light));
	lamp->name = malloc(sizeof(char) * 30);

	//Reading the name and the position of the light. Breaking if reading fails. 
	result = fscanf(fp," %29s %lf %lf %lf ",lamp->name,&lamp->posX,&lamp->posY,&lamp->posZ);
	if(result != 4)return NULL;

	//Reading the intensity values of the light, and breaking if there is an error with the read. 
	result = fscanf(fp," %lf %lf %lf ",&lamp->iR,&lamp->iG,&lamp->iB);
	if(result != 3)return NULL;

	return lamp;
}

//Function for printing the parsed light. 
void printLight(light *lamp){
	printf("Light name: %s Position(%.1lf %.1lf %.1lf)\n",lamp->name,lamp->posX,lamp->posY,lamp->posZ);
	printf("\t Intensity red: %.1lf green: %.1lf blue :%.1lf\n",lamp->iR,lamp->iG,lamp->iB);
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
	printf("Back Color: (%.1lf %.1lf %.1lf)\n",r,g,b); //Color of the background
	printf("Ambient Light: red: %.1lf green: %.1lf blue: %.1lf\n",aR,aG,aB); //Ambient light
	printf("\nOutput File: %s\n\n",outputFile); //Name of the output file.
	printf("~~~~~ Completed reading file ~~~~~\n\n");
}
//Second parser, if only to comply with the 'testParser' test case.
//Reads the necessary inputs from the file specified by the string.
//The inputs are stored in global variables. 

//Returns -1 if the parse was unsucessful.
int parse2(char *fileName){
	//Flags for whether the viewing planes have been parsed or not.
	int readNear, readLeft,readRight,readTop,readBottom;
	readNear = 0;readLeft = 0;readRight = 0;readTop = 0;readBottom = 0;
	
	//Flags for whether the other necessary inputs have been parsed or not.
	int readRes, readBack, readAmbient, readOutput;
	readRes = 0;readBack = 0;readAmbient = 0;readOutput = 0;


	//Initializing lists
	sphereList = malloc(sizeof(sphere *) * 15);
	lightList = malloc(sizeof(light *) * 10);
	
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
		//Read a light, up to 10 times. 
		else if(strcmp(input,"LIGHT") == 0){
			if(numLights >= 10){
				fprintf(stderr,"::Too many lights. \n::There is a limit of 10 lights for this assignment.\n\n");
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
		//If one of these input descriptors was not read, the file format is not valid. Exit. 
		else{
			return -1;
		}
		result = fscanf(file," %s ",input);
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
	return 0;		
}
//Function to convert rgb floats into an rgb color. 
int convertIntColor(double red, double green, double blue){
	if(red > 1.0 || green > 1.0 || blue > 1.0){
		printf("COLOR OVERFLOW %lf %lf %lf\n",red,green,blue);
	}
	int color = (int)(red * 255) << 16;
       	color += (int)(green * 255) << 8;
	color += (int)(blue * 255); 
	return color; 	
}

//Function to create the image array for the ray tracer. 
void createImageArray(){
	int row, col;
	imageArray = malloc(sizeof(int *) * rows);
	for(row = 0;row < rows;++row){
		imageArray[row] = malloc(sizeof(int) * cols); 
	}
}
//Function for dumping the completed image into a byte buffer.
//the 2D integer array for the image is freed. 
char *convertImageToChars(){
	char* byteBuffer = malloc(sizeof(char) * (rows * cols * 3));
	int x, y;
	int index;
	int color;
	//Buffering the image into a buffer array of 24 bit colors. 
	for(y = 0;y < rows;++y){
		for(x = 0;x < cols;++x){
			index = y * cols * 3 + x * 3;
			color = imageArray[y][x];
			byteBuffer[index] = (color >> 16) % 256; //Red
			byteBuffer[index + 1] = (color >> 8) % 256; //Green 
			byteBuffer[index + 2] = color % 256; //Blue
		}
		//Freeing the current row of the image buffer since we don't need it anymore.
		free(imageArray[y]);
	}
	//Freeing the image array pointer since we're done with it. 
	free(imageArray);
	return byteBuffer;
}
//Function to save the image
//Copied from the website. 
void save_image(int Width, int Height, char* fname,unsigned char* pixels) {
	FILE *fp;
	const int maxVal=255; 
  
	if(verbose)printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"wb");
	if(!fp){
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);

	for(int j = 0; j < Height; j++) {
		fwrite(&pixels[j*Width*3], 3,Width,fp);
	}
	fclose(fp);
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

//Checks if there exists a collision
//in the way of the path of this ray from the origin (the collision point) to light. 
int existsCollision(Matrix *origin,Matrix *ray){
	sphere *s;
	double t;
	for(int i = 0;i < numSpheres;++i){
		s = sphereList[i];
		t = computeTToSphere(ray,origin,s,MINIMUM_T);
		//If there is indeed an n such that 0 <= t <= 1, there is a collision between a sphere and the shadow ray.
		if(t >= MINIMUM_T && t <= 1.0 - MINIMUM_T){
			return 1;
		}
	}
	return 0;
}

//Computing the light color at the given point of collision.
//the colPoint is the point of collision with the surface. 
//the normal is the normal vector of collision with the surface
//s is the sphere. 
//origin is the observing point - but this changes as rays are traced.
//returns the sum of all the light colors. 
color *computeLightColor(Matrix *colPoint,Matrix *origin,Matrix *normal,sphere *s){
	//Initializing color to zero.
	color *c = malloc(sizeof(color));
	c->r = 0.0;
	c->g = 0.0;
	c->b = 0.0;

	light *l;

	//The ray from the collision point to the light source.
	Matrix *ray;
	//The point of the light source.
	Matrix *lightPoint;
	for(int i = 0;i < numLights;++i){
		l = lightList[i];
		lightPoint = point4(l->posX,l->posY,l->posZ);

		//Computing a vector from the surface of the sphere to the light. 
		ray = matrixCopy(lightPoint);
		inPlaceDifference(ray,colPoint);

		//Now, need to ensure there aren't any spheres in the way.
		//In other words, need to check if there is a sphere in between the colPoint and the lightPoint
		if(!existsCollision(colPoint,ray)){
			//Need to compute light, as long as there was a true collision.
			//Getting the dot product of the normal and the light ray. 
			double dot = dotProduct(ray,normal);
			
			if(dot >= 0){
				//Color placeholders. 
				double specR,specG,specB;
				double difR,difG,difB;
				//The "r" vector;the light ray after it has been reflected on the surface using the normal.
				Matrix *ref;
				//The "v" vector;the viewing ray from the eye to the collision point. 
				Matrix *viewRay;
				Matrix *projection;
				double shininess;
				double spec,diff;

				//Diminishing the dot product between the ray.
				dot /= sqrt(dotProduct(ray,ray));
				dot /= sqrt(dotProduct(normal,normal));

				//Computing the diffuse light. 
				diff = dot * s->kDif;
				difR = diff * l->iR * s->r;
			       	difG = diff * l->iG * s->g;
				difB = diff * l->iB * s->b;

				//Computing the reflection of the light ray on the surface. 
				ref = getScalarMultipleMatrix(ray,-1);
				dot = 2*dotProduct(ray,normal) / dotProduct(normal,normal);
				projection = getScalarMultipleMatrix(normal,dot);
				inPlaceSum(ref,projection);

				//Computing the viewing ray... 
				//The viewing ray should be the vector from the collision point to the eye. 
				viewRay = getScalarMultipleMatrix(colPoint,-1);
				inPlaceSum(viewRay,origin);
				toVector(viewRay);

				//Computing the shininess coefficient. 
				shininess = dotProduct(ref,viewRay);
				//If the shininess is zero, the reflected light ray and the view ray do not coincide. 
				//So there is no specular light. 
				if(shininess >= 0 && s->kSpec > 0.0){
				//shininess = fabs(shininess);
					shininess /= sqrt(dotProduct(ref,ref));
					shininess /= sqrt(dotProduct(viewRay,viewRay));
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
				c->r += (difR + specR);
				c->g += (difG + specG);
				c->b += (difB + specB);
				
				freeMatrix(projection);
				freeMatrix(viewRay);
				freeMatrix(ref);
			}
		}
		freeMatrix(ray);
		freeMatrix(lightPoint);
	}
	
	return c;
}
//Traces the ray to the closest sphere, if possible, then computes the color.
color *traceRay(Matrix *ray,Matrix *origin,int bounceCount){
	color *c = malloc(sizeof(color));
	sphere *s = NULL;
	double t;
	double lowestT = 2020202020202020;
	int col;
	color *lightColor;

	//Computing the closest sphere intersection with the ray. 
	double minimum = MINIMUM_T;
	if(bounceCount == NUM_BOUNCES){
		minimum = 1 + MINIMUM_T;
	}
	for(int i = 0;i < numSpheres;++i){
		//Need to compute the normal here...
		t = computeTToSphere(ray,origin,sphereList[i],minimum);
		//If t is smaller than the lowest T so far, we take it. 
		if(t > minimum && t < lowestT){
			lowestT = t;
			s = sphereList[i];
		}
	}
	t = lowestT;
	//There's a collision if we get inside this statement....
	//We'll start making recursive calls...
	if(s != NULL){
		//Variables for making the color calculations of the pixel.
		Matrix *normal = NULL;
		Matrix *tmp = NULL;
		Matrix *inverse = NULL;
		Matrix *inverseTranspose = NULL;
		Matrix *colPoint = getScalarMultipleMatrix(ray,t);
		inPlaceSum(colPoint,origin);
		Matrix *rayPrime = NULL;
		Matrix *originPrime = NULL;
		//Computing the ambient light.
		c->r = s->kAmb * aR * s->r;
		c->g = s->kAmb * aG * s->g;
		c->b = s->kAmb * aB * s->b;	
		
		//Computing the inverse matrix of the sphere's transformation matrix. 
		inverse = getSphereMatrix(s);
		tmp = inverse;
		inverse = getInverseMatrix(inverse);
		freeMatrix(tmp);
		
		//Building the normal vector...
		rayPrime = getProductMatrix(inverse,ray);
		originPrime = getProductMatrix(inverse,origin);
	
		//Getting the sum of the ray prime and the origin prime, then using the inverse transpose to generate the actual normal vector. 
		//Getting the collision point with respect to the canonical sphere - origin plus t times the ray. 
		inPlaceScalarMultiply(rayPrime,t);
		inPlaceSum(rayPrime,originPrime);
		toVector(rayPrime); //We want the vector with respect to the origin - but the origin is 0,0,0 so we can just knock off the point value. 
		inverseTranspose = matrixCopy(inverse);
		inPlaceTranspose(inverseTranspose);
		//The normal is inverse transpose applied to the collision point minus (0,0,0) as a vector - but subtracing (0,0,0) is redundant so we skip that.  
		normal = getProductMatrix(inverseTranspose,rayPrime);
		toVector(normal);
		//If the ray from the origin to collision point is longer than the vector from the origin to the center of the sphere,
		//The normal should be flipped. 
		inPlaceDifference(rayPrime,originPrime);
		rayPrime->matrix[3][0] = 0.0;
		if(dotProduct(rayPrime,rayPrime) > dotProduct(originPrime,originPrime)){
			inPlaceScalarMultiply(normal,-1);
		}

		//Light collision methods go here.
		lightColor = computeLightColor(colPoint,origin,normal,s);
		c->r += lightColor->r;
		c->g += lightColor->g;
		c->b += lightColor->b;
		

		//If the ray should be reflected, reflect it. Here goes nothing.
		if(bounceCount > 0 && s->kR > 0){
			Matrix *reflectedRay;
			Matrix *projection;
			color *reflectedColor;
			--bounceCount;

			//calculating the projection of the ray onto the normal. 
			projection = getScalarMultipleMatrix(normal,2*((dotProduct(ray,normal) / dotProduct(normal,normal))));

			//Calculating the reflected ray. 
			reflectedRay = matrixCopy(ray);
			inPlaceDifference(reflectedRay,projection);

			reflectedColor = traceRay(reflectedRay,colPoint,bounceCount);
			c->r += (s->kR * reflectedColor->r);
			c->g += (s->kR * reflectedColor->g);
			c->b += (s->kR * reflectedColor->b);
			
			//Freeing the resources used here. 
			freeMatrix(projection);
			free(reflectedColor);
			freeMatrix(reflectedRay);
		}

		freeMatrix(inverse);
		freeMatrix(rayPrime);
		freeMatrix(originPrime);
		freeMatrix(inverseTranspose);
		free(lightColor);
		freeMatrix(colPoint);
		freeMatrix(normal);
	}

	else{
		//If this is a bounced ray, return black if there is no collision.
		if(bounceCount < NUM_BOUNCES){
			c->r = 0;
			c->g = 0;
			c->b = 0;
		}
		//Otherwise, return the background color. 
		else{
			c->r = r;
			c->g = g;
			c->b = b;
		}
	}
	return c;
}
//Traces the given ray to the given sphere,
//from the given starting point. 
double computeTToSphere(Matrix *ray,Matrix *origin,sphere *s,double minimum){
	double a,b,c;
	double t = -1;
	double det;
	//Need to find the distance to the sphere, if one exists. 

	//Obtains the matrix of the sphere.
	Matrix *m = getSphereMatrix(s);

	//Getting the inverse matrix of the transformation, m.
	inPlaceInverse(m);

	//Applying the matrix to the ray. 
	ray = getProductMatrix(m,ray);

	//Applying the matrix to the origin of the vector.
	origin = getProductMatrix(m,origin);

	//Freeing the inverse matrix, since we're done with it. 
	freeMatrix(m);
	
	a = dotProduct(ray,ray);
	b = dotProduct(origin,ray);
	c = dotProduct(origin,origin) - 1;

	det = b * b - a * c;
	//If there is a collision between the sphere and the ray...
	freeMatrix(ray);
	freeMatrix(origin);
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
//Starting to compute the color of each pixel,
//and placing that pixel in the array of pixels. 
void computePixels(){
	//The ray should go directly through the middle of each pixel. 
	double zeroX = (-cols / 2.0) + 0.5;
	double zeroY = (-rows / 2.0) + 0.5;
	//Used for converting the ray's x and y values into viewing coordinates. 
	double planeX = (right - left) / (cols);
	double planeY = (bottom - top) / (rows);

	int x, y;

	//Adjusting the pixel at (x,y) from the top left corner. 
	//For each row of pixels...
	Matrix *eye;
	Matrix *ray;
	double rayX, rayY, rayZ;
	for(y = 0;y < rows;++y){
	//for(y = rows / 2;y < rows / 2 + 1;++y){
		//For each pixel within the row... 
		for(x = 0;x < cols;++x){
		//for(x = cols / 2;x < cols / 2 + 1;++x){
			//Getting the ray vector of the current pixel.
			//This is now the vector from the eye into the pixel. 
			rayX = (x + zeroX) * planeX;
			rayY = (y + zeroY) * planeY;
			rayZ = -near;
			//Converting the color that has been traced into an integer. 
			ray = vec4(rayX,rayY,rayZ);
			eye = point4(0,0,0);

			//Computing the pixel color.
			color *pxColor = traceRay(ray,eye,NUM_BOUNCES);
			
			//Clamping the color, if the color has exceeded one. 
			if(pxColor->r > 1){
				pxColor->r = 1;
			}
			if(pxColor->g > 1){
				pxColor->g = 1;
			}
			if(pxColor->b > 1){
				pxColor->b = 1;
			}
			imageArray[y][x] = convertIntColor(pxColor->r,pxColor->g,pxColor->b);

			freeMatrix(eye);
			eye = NULL;
			freeMatrix(ray);
			ray = NULL;
			free(pxColor);
			pxColor = NULL;
		}
	}
}
//Freeing the sphereList and the spheres it contains
//and the light list and the lights contained. 
void freeLists(){
	int i;
	for(i = 0;i < numSpheres;++i){
		free(sphereList[i]);
	}
	free(sphereList);
	for(i = 0;i < numLights;++i){
		free(lightList[i]);
	}
	free(lightList);
}
//Main function of the program. 
int main(int argc,char **argv){
	//Checking if the program lacks arguments. 
	if(argc < 2){
		printf("Please retry running the program with the correct number of arguments\n");
		return 1;
	}
	
	//If the parsing option appears first, swap the command line arguments so that the actual file appears first. 
	if(argc > 2 && strcmp(argv[1],"-v") == 0){
		char *swap = argv[1];
		argv[1] = argv[2];
		argv[2] = swap;
	}

	//Parsing the file and checking for success. 
	int result = parse2(argv[1]);
	//Stop program with an error if the parsing failed. 
	if(result == -1){
		fprintf(stderr,"::Parsing the file failed.\n::Please ensure you have provided a valid txt file in the format specified by the assignment.\n");
		return 1;
	}

	//Printing the output of the parser, if requested. 
	if(argc > 2 && (strcmp(argv[2],"1") == 0 || strcmp(argv[2],"-v") == 0)){
		verbose = 1;
		printParsedFile(argv[1]);
	}
	createImageArray();

	computePixels();
	//Freeing the lists of lights and spheres. 
	freeLists();
	//Saving the image file:
	byteBuffer = convertImageToChars();
	save_image(cols,rows,outputFile,byteBuffer);

	//Freeing the remaining used resources. 
	free(byteBuffer);
	free(outputFile);
	return 0;
}
