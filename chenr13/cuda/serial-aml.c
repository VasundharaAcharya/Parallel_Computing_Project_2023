/*AML DATASET: This is the serial version of the code that uses the euclidean distance as the distance metric. 
Everything else is similar to the other code. We also have the Silhoutte score computation logic in the same code.
You can vary the K value to obtain results of different cluster values.*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "clockcycle.h"

#define INF 0x3f3f3f3f

#define Genes 7129     // X Total Number of genes to be given as an input. 
#define Samples 34          // Represents the sample genes
// Change the value of K to obtain the results with different clusters
#define K 8 // Number of clusters

const int blockSize = 32;
const int Iteration = 1000;

int *cluster_idx;             
double *gene_data;          
double *medoids;      

#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif


void init();
void readdata();
void initmedoids();
void finalize();

// clustering algorithm
double manhattan_distance(double *p1, double *p2, int dim);
double euclidean_distance(double *p1, double *p2, int dim);
void findclosestmedoids(double *gene_data, double *medoids, int* idx);
void computeMedoids(double* gene_data, int* idx, double* medoids);
double average_cluster_distance(double **data, int *labels, int numdata, int dim, int idx);
double average_nearest_cluster_distance(double **data, int *labels, int numdata, int k, int dim, int idx);
double silhouette_coefficient_single_point(double **data, int *labels, int numdata, int k, int dim, int idx);
double silhouette_coefficient(double **data, int *labels, int numdata, int k, int dim);


int main() {
    // Define variables to keep track of time
    unsigned long long start = 0;
    unsigned long long finish = 0;

    int i, j;
    int *label, lab;
    double **data, temp;

    init();

    readdata();

    //Start the clockcycle
    start = clock_now();

    initmedoids();

    for (i = 0; i < Iteration; i ++ ) {
        findclosestmedoids((double*)gene_data, (double *)medoids, &cluster_idx[0]);
        computeMedoids((double *)gene_data, &cluster_idx[0], (double *)medoids);
    }

    for (i = 0; i < K; i ++ ) {
        for (j = 0; j < Samples; j ++ ) {
            printf("Final Medoids are  %lf  ", *(medoids + i * Samples + j));
        }
        printf("\n");
    }


    printf("kmedoids completed with clustering of the gene samples \n");

    FILE *fp, *fw;
    // Write clustered gene data to file
    int flag1 = 0;
    char filename3[100];
    sprintf(filename3, "output_gene_data_%d_euclidean.txt", K);
    fw = fopen(filename3, "w");
    for (i = 0; i < Genes; i ++ ) {
        for (j = 0; j < Samples; j ++ ) {
            flag1 = fprintf(fw, "%lf  ", *(gene_data + i * Samples + j));
        }
        fprintf(fw, "\n");
    }
    fclose(fw);

    if (flag1) { printf("Successfully wrote the clustered gene data to the file \n"); }

    // Write gene medoids to file
    int flag = 0;
    char filename2[100];
    sprintf(filename2, "output_gene_medoids_%d_clusters_euclidean.txt", K);
    fw = fopen(filename2, "w");
    for (i = 0; i < K; i ++ ) {
        for (j = 0; j < Samples; j ++ ) {
            flag = fprintf(fw, "%lf  ", *(medoids + i * Samples + j));
        }
        fprintf(fw, "\n");
    }
    fclose(fw);

    if (flag) { printf("Successfully wrote the gene medoids to the file \n"); }


    //Write the cluster assignments to the file
    int flag2 = 0;
    char filename[100];
    sprintf(filename, "gene_med_ind_%d_clusters_euclidean.txt", K);
    fw = fopen(filename, "w");
    for (i = 0; i < Genes; i ++ ) {
        flag2 = fprintf(fw, "%d  ", cluster_idx[i]);
    }
    fclose(fw);
    if(flag2) { printf("Successfully wrote the cluster assignments to the file \n"); }

    finish = clock_now();

    unsigned long long starting = 0;
    unsigned long long finishing = 0;

    // Start of silhoutte coefficient. We are using this as we cannot plot these huge data.
    // To validate the results of the clustering algorithm, we use this. 

    starting = clock_now();
    data = (double**)calloc(Genes, sizeof(double*));

    for (int i = 0; i < Genes; i ++ ) {
        data[i] = (double*)calloc(Samples, sizeof(double));
    }

    label = (int*)calloc(Genes, sizeof(int));

    fp = fopen(filename3, "r");
    if (fp == NULL) {
        printf("Error, the requested file does not exist\n");   
        exit(1);             
    }
    for (int i = 0; i < Genes; i ++ )
    {
        for(int j = 0; j < Samples; j ++ )
        {
            fscanf(fp, "%lf", &temp);
            data[i][j] = temp;
        }
    }

    fclose(fp);

    fp = fopen(filename, "r");
    if(fp == NULL) {
        printf("Error, the requested file does not exist");   
        exit(1);             
    }
    for(int i = 0; i < Genes; i ++ ) {
        fscanf(fp,"%d", &lab);
        label[i] = lab;
    }
    fclose(fp);

    // end of sill coeff
    printf("sillhouette coefficient:%lf\n", silhouette_coefficient(data, label, Genes, K, Samples));

    // End the timer here
    finishing = clock_now();
    printf("Total clockcycles taken to run K-Medoids on %d genes with %d clusters is %lld cycles.\n", Genes, K, (finish - start));
    printf("Total time taken to run K-Medoids on %d genes with %d clusters is %e seconds.\n", Genes, K, (finish - start) / 512000000.0f);
    printf("Total clockcycles taken to compute silhouette score on %d genes with %d clusters is %lld cycles.\n", Genes, K, (finishing - starting));
    printf("Total time taken to compute silhouette score on %d genes with %d clusters is %e seconds.\n", Genes, K, (finishing - starting) / 512000000.0f);

    finalize();

    return EXIT_SUCCESS;
}


void init() { // initialization
    srand(time(0));

    gene_data = calloc(Genes * Samples, sizeof(double));
    medoids = calloc(K * Samples, sizeof(double));
    cluster_idx = calloc(Genes, sizeof(int));
}


void readdata() { // read data from "training.txt" in serial
    FILE *fp;
    double num1;
    fp = fopen("training.txt", "r");
    if (fp == NULL) {
        printf("The requested input file does not exist. \n");
        exit(1);
    }
    for (int i = 0; i < Genes; i ++ ){
        for (int j = 0; j < Samples; j ++ ){
            fscanf(fp, "%lf", &num1);
            *(gene_data + i * Samples + j) = num1;
        }
    }
    fclose(fp);
}


void initmedoids() {
    //Creation of random number and start with random medoids
    int rnd_num;
    for (int i = 0; i < K; i ++ ) {
        rnd_num = rand() % Genes;
        for (int j = 0; j < Samples; j ++ ){ 
            *(medoids + i * Samples + j) = *(gene_data + rnd_num * Samples + j); 
        }
    }
    // Print initial medoids
    // for (i = 0; i < K; i++) {
    //     printf("Initial Medoid %d: ", i+1);
    //     for (j = 0; j < Samples; j++) {
    //         printf("%f ", *(medoids + i * Samples + j));
    //     }
    //     printf("\n");
    // }
}


void finalize() {
    free(gene_data);
    free(medoids);
    free(cluster_idx);
}



// Specification of the distance metrics for the task of gene clustering.
// We use Manhattan Distance here
double manhattan_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i ++ ) {
        distance += fabs(p1[i] - p2[i]) + fabs(p1[i] - p2[i]);
    }
    return distance;
}


double euclidean_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i ++ ) {
        distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(distance);
}


// This function is used to compute the closest medoids 
void findclosestmedoids(double *gene_data, double *medoids, int* idx) {
    int i, j, l;
    double sum, dist[K], min_dist;

    for (i = 0; i < Genes; i ++ ) {
        min_dist = INF;

        for (j = 0; j < K; j ++ ) {
            sum = 0;

            // compute Euclidean distance
            for (l = 0; l < Samples; l ++ ) {
                double diff = gene_data[i * Samples + l] - medoids[j * Samples + l];
                sum += diff * diff;
            }

            dist[j] = sqrt(sum);

            if (dist[j] < min_dist) {
                min_dist = dist[j];
                idx[i] = j;
            }
        }
    }
}


void computeMedoids(double* gene_data, int* idx, double* medoids) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;

    // set the sample size as needed here
    int sample_size = 10; 
    for (i = 0; i < K; i ++ ) {
        min_distance = 1e9;
        // randomly sample a subset of points from the cluster
        int* sample = calloc(sample_size, sizeof(int));
        for (j = 0; j < sample_size; j ++ ) {
            sample[j] = -1;
            while (sample[j] == -1 || idx[sample[j]] != i) {
                sample[j] = rand() % Genes;
            }
        }
        // calculate distances only for the sampled points
        // I have used the manhattan distance, we can experiment with others as well
        for (j = 0; j < sample_size; j ++ ) {
            sum = 0.0;
            for (k = 0; k < sample_size; k ++ ) {
                distance = 0.0;
                // compute Euclidean distance
                for (l = 0; l < Samples; l ++ ) {
                    double diff = *(gene_data + sample[j] * Samples + l) - *(gene_data + sample[k] * Samples + l);
                    distance += diff * diff;
                }
                sum = sqrt(distance);
            }
            if (sum < min_distance) {
                min_distance = sum;
                for (m = 0; m < Samples; m++) {
                    temp = *(gene_data + sample[j] * Samples + m);
                    *(medoids + i * Samples + m) = temp;
                }
            }
        }
        free(sample);
    }
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


// Calculate the average nearest cluster distance for the data point with index 'idx' in 'data'
//  to its nearest cluster
// argument 'k' is the number of clusters
double average_nearest_cluster_distance(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double cluster_distance[k];
    //To fill the values to cluster_distance with zero
    memset(cluster_distance, 0.0, sizeof cluster_distance);
    int count[k];
    memset(count, 0, sizeof count);
    for (int i = 0; i < numdata; i ++ ) {
        cluster_distance[labels[i]] +=euclidean_distance(data[idx], data[i], dim);
        count[labels[i]] ++ ;
    }
    double res = INF;
    for (int i = 0; i < k; i ++ ) {
        if (i == labels[idx]) continue;
        res = MIN(res, cluster_distance[i] / count[i]);
    }
    return res;
}


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
        totalsc += silhouette_coefficient_single_point(data, labels, numdata, k, dim, i);
    }
    return totalsc / numdata;
}

