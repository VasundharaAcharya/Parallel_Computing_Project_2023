
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "clockcycle.h"

#define X 7129     // X Total Number of genes to be given as an input. 
#define Y 34          // Represents the sample genes
#define K 15 // NUMBER OF CLUSTERS TO DIVIDE THE DATA INTO
#define MAX_ITERS 10  // NUMBER OF ITERATIONS in k K-medoids
#define DBL_MAX 1e9
double *num;          // pointer to all the data to be stored
double *medoids;      // pointer to the data of medoids
int *idx;             // pointer to data  which stores the index of the medoid nearest to each pixel

typedef unsigned long long ticks;





void findclosestmedoids(double *num, double *medoids, int* idx) {
    int i, j, l;
    double sum, dist[K], min_dist;

    for (i = 0; i < X; i++) {
        min_dist = INFINITY;

        for (j = 0; j < K; j++) {
            sum = 0;

            for (l = 0; l < Y; l++) {
              //Manhattan distance
                sum += fabs(num[i * Y + l] - medoids[j * Y + l]);
            }

            dist[j] = sum;

            if (dist[j] < min_dist) {
                min_dist = dist[j];
                idx[i] = j;
            }
        }
    }
}



void computeMedoids(double* num, int* idx, double* medoids) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;
    int sample_size = 10; // set the sample size to 1000
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
        // I have used the manhattan distance, we can experiment with others as well
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



double euclidean_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i ++ ) {
        distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(distance);
}


// Calculate the average distance for the data point with index 'idx' in 'data'
//  within its cluster
double average_cluster_distance(double **data, int *labels, int numdata, int dim, int idx) {
    double distance_sum = 0.0;
    int count = 0;
    for (int i = 0; i < numdata; i ++ ) {
        if (labels[i] != labels[idx]) continue;
        distance_sum += euclidean_distance(data[idx], data[i], dim);
        count ++ ;
    }
    return distance_sum / count;
}


#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

// Calculate the average nearest cluster distance for the data point with index 'idx' in 'data'
//  to its nearest cluster
// argument 'k' is the number of clusters
double average_nearest_cluster_distance(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double cluster_distance[k];
    memset(cluster_distance, 0.0, sizeof cluster_distance);
    int count[k];
    memset(count, 0, sizeof count);
    for (int i = 0; i < numdata; i ++ ) {
        cluster_distance[labels[i]] += euclidean_distance(data[idx], data[i], dim);
        count[labels[i]] ++ ;
    }
    double res = INFINITY;
    for (int i = 0; i < k; i ++ ) {
        if (i == labels[idx]) continue;
        res = MIN(res, cluster_distance[i] / count[i]);
    }
    return res;
}


#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif

// Calculate silhouette coefficient for the data point with index 'idx' in 'data'
double silhouette_coefficient_single_point(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double a = average_cluster_distance(data, labels, numdata, dim, idx);
    double b = average_nearest_cluster_distance(data, labels, numdata, k, dim, idx);
    return (b - a) / MAX(a, b);
}

// Calculate average silhouette coefficient
double silhouette_coefficient(double **data, int *labels, int numdata, int k, int dim) {
    double totalsc = 0.0;
    for (int i = 0; i<numdata; i ++ ) {
        if(i%1000==0)
            printf("%d\n",i);
        totalsc += silhouette_coefficient_single_point(data, labels, numdata, k, dim, i);
    }
    return totalsc / numdata;
}


int main() {
    // Define variables to keep track of time
    unsigned long long start = 0;
    unsigned long long finish = 0;

    FILE *fp, *fw;
    double num1;
    int i, j, k, rnd_num;

 int n =X, fs = Y, k1=K, *label,lab;

    double **data,temp;



    num = (double*) calloc(X * Y, sizeof(double));
    medoids = (double*) calloc(K * Y, sizeof(double));
    idx = (int*) calloc(X, sizeof(int));


    // Upper and lower bounds for the random numbers 
    int lower = 0;
    int upper = X - 1;
    
    srand(time(0));
    
    // File open for reading
    fp = fopen("training.txt", "r");
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
		
			for (j=0;j<Y;j++){ 
        		*(medoids+i*Y+j) = *(num+rnd_num*Y+j); 
        	} 
      
    }
  

  
  for (i=0; i<MAX_ITERS; i++){

	findclosestmedoids((double*)num, (double *)medoids, &idx[0]);
 


	computeMedoids((double *)num, &idx[0], (double *)medoids);

	}


	for(i=0; i<K;i++){
			for(j=0; j<Y;j++){

				printf("Medoids are tt %lf  ",*(medoids+i*Y+j));

			}
		printf("\n");
	}
	








	// kmedoids algorithm completed
  printf("kmedoids completed");

	// Close timer and print total time taken
	finish = clock_now();
	printf("Total time taken to run K-Medoids on %d genes with %d clusters is %e seconds.\n", X, K, (finish-start)/512000000.0f);





// Write clustered gene data to file
fw=fopen("output_gene_data.txt","w");
for(i=0; i<X;i++){
    for(j=0; j<Y;j++){
        fprintf(fw,"%lf  ",*(num+i*Y+j));
    }
    fprintf(fw, "\n");
}
fclose(fw);



// Write gene medoids to file
fw=fopen("output_gene_medoids.txt","w");
for(i=0; i<K;i++){
    for(j=0; j<Y;j++){
        fprintf(fw,"%lf  ",*(medoids+i*Y+j));
    }
    fprintf(fw, "\n");
}
fclose(fw);



fw = fopen("gene_med_ind.txt","r");
for (i = 0; i < X; i++) {
    // printf("%d ", idx[i]);
    fprintf(fw,"%d  ",idx[i]);

}


// start of sill coeff




    data = (double**)malloc(n*sizeof(double*));

    for(int i=0;i<n;i++)
    {
        data[i]=(double*)malloc(fs*sizeof(double));
    }

    label = (int*)malloc(n*sizeof(int));


    fp = fopen("output_gene_data.txt", "r");

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<fs;j++)
        {
            fscanf(fp,"%lf",&temp);
            data[i][j]=temp;
        }
    }

    fclose(fp);

    fp = fopen("gene_med_ind.txt", "r");

    for(int i=0;i<n;i++)
    {

        fscanf(fp,"%d",&lab);
        label[i]=lab;
    }
    fclose(fp);

    printf("sillhouette coefficient:%lf\n",silhouette_coefficient(data,label,n,k1,fs));



// end of sill coeff



	free (num);
	free (medoids);
	free(idx);
}

