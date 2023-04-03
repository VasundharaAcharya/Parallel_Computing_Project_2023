#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "/content/clockcycle.h"

#define X 1049088     // X Total Number of pixels in the image, to be given as an input. 
#define Y 3           // No change required basically the RBG components of an image
#define K 4         // NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 10  // NUMBER OF ITERATIONS in k K-medoids

double *num;          // pointer to all the data to be stored
double *medoids;      // pointer to the data of medoids
int *idx;             // pointer to data  which stores the index of the medoid nearest to each pixel

typedef unsigned long long ticks;


//This function appears to work fine
void findclosestmedoids(double *num, double *medoids, int* idx){
    
    int i, j, l, min_ind; 
    double sum, dist[K],min_dist;
    for (i=0;i<X;i++){
        for (j=0;j<K;j++){
            sum=0;
            for (l=0;l<Y;l++){
              //manhattan distance
                sum = sum + fabs(*(num+i*Y+l)-*(medoids+j*Y+l));
            }
            dist[j]=sum;
        }
        min_dist=dist[0];
        min_ind=0;
        for (j=0; j<K; j++){
            if (dist[j]<min_dist) {
                min_dist=dist[j];
                min_ind=j;
            }
        }
        *(idx+i)=min_ind;
    }
}





/*It iterates over each medoid, and then for each medoid it iterates over all the data points that belong to that medoid. For each combination of data points that belong to that medoid, it calculates the sum of Manhattan distances between them. The data point that has the smallest sum of Manhattan distances to all the other points in the same medoid is chosen as the new medoid for that cluster.

The new medoid is then assigned to the corresponding index in the medoid array.*/

//the computeMedoids function randomly samples a subset of points from each cluster and calculates the distances only for those points. 
//Currently, the computeMedoid function calculates the distance between each point and every other point in the cluster. This results in a large number of distance calculations, which can be computationally expensive. One way to reduce the number of distance calculations 
//is to randomly sample a subset of points from the cluster and calculate the distances only for those points.
void computeMedoids(double* num, int* idx, double* medoids) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;
    int sample_size = 1000; // set the sample size to 10
    for (i = 0; i < K; i++) {
        min_distance = DBL_MAX;
        // randomly sample a subset of points from the cluster
        int* sample = (int*)malloc(sample_size * sizeof(int));
        for (j = 0; j < sample_size; j++) {
            sample[j] = -1;
            while (sample[j] == -1 || idx[sample[j]] != i) {
                sample[j] = rand() % X;
            }
        }
        // calculate distances only for the sampled points
        for (j = 0; j < sample_size; j++) {
            sum = 0.0;
            for (k = 0; k < sample_size; k++) {
                distance = 0.0;
                for (l = 0; l < Y; l++) {
                  //manhattan distance
                    distance += fabs(*(num + sample[j] * Y + l) - *(num + sample[k] * Y + l));
                }
                sum += distance;
            }
            if (sum < min_distance) {
                min_distance = sum;
                for (m = 0; m < Y; m++) {
                    temp = *(num + sample[j] * Y + m);
                    *(medoids + i * Y + m) = temp;
                }
            }
        }
        free(sample);
    }
}


int main() {
  //some print are added to see if its running
    //printf("here");
    // Define variables to keep track of time
    unsigned long long start = 0;
    unsigned long long finish = 0;

    FILE *fp, *fw;
    double num1;
    int i, j, k, rnd_num;
    num = (double*) calloc(X * Y, sizeof(double));
    medoids = (double*) calloc(K * Y, sizeof(double));
    idx = (int*) calloc(X, sizeof(int));

    // Upper and lower bounds for the random numbers 
    int lower = 0;
    int upper = X - 1;
    
    srand(time(0));
    
    // File open for reading
    fp = fopen("input.txt", "r");
    if (fp == NULL) {
        printf("Exiting no file with such name \n");
        exit(-1);
    }

   for (i=0;i<X;i++){
		for (j=0;j<Y;j++){
			fscanf(fp,"%lf", &num1);
			*(num+i*Y+j)=num1;
		}
	}
	fclose(fp);

	start = clock_now();


	//Creation of random number and start with random medoids
	for (i = 0; i < K; i++) {

			rnd_num = (rand()%(upper-lower + 1)) + lower;
			//printf("%d ", rnd_num);  
			for (j=0;j<Y;j++){ 
        		*(medoids+i*Y+j) = *(num+rnd_num*Y+j); 
        	} 
       // printf("thy");
    }
  

  
  for (i=0; i<MAX_ITERS; i++){

	findclosestmedoids((double*)num, (double *)medoids, &idx[0]);
 

  //call to computeMedoids works to be okay
	computeMedoids((double *)num, &idx[0], (double *)medoids);

	}


	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("Medoids are tt %lf  ",*(medoids+i*Y+j));

			}
		printf("\n");
	}
	



  //here instead of this we need to cluster the genes.

	//Repalcement of each pixel in the image by the medoid it is closest to. 

	for (i=0; i<X;i++){
		//printf("%d==%d\n",i+1, idx[i]+1);

		for (k=0;k<K;k++){

			if (idx[i]==k){

					for (j=0;j<Y;j++){			
						*(num+i*Y+j)=*(medoids+k*Y+j);
					}
			}
				
		}

	}

	// KMeans algorithm completed
  printf("kmedoids completed");

	// Close timer and print total time taken
	finish = clock_now();
	printf("Total time taken to run K-Medoids on %d pixel image with %d clusters is %e seconds.\n", X, K, (finish-start)/512000000.0f);


	//Writing the data to the output file
	fw=fopen("output.txt","w");
	
	for(i=0; i<X;i++){
	
		for(j=0; j<Y;j++){

				fprintf(fw,"%lf  ",*(num+i*Y+j));
			
			}
		fprintf(fw, "\n");
		//printf("\n");
	}


//printf("thy");
	free (num);
	free (medoids);
	free(idx);
}