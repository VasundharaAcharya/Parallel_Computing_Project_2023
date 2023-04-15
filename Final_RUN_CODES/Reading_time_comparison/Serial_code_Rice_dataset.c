/* This is the time needed to read the rice dataset. Serial code version */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "clockcycle.h"

#define Genes 10000   // X Total Number of genes to be given as an input. 
#define Samples 35         // Represents the sample genes
//Change the value of K to obtain the results with different clusters
#define K 2 // Number of clusters
//Initializations
int *cluster_idx;             
double *gene_data;          
double *medoids;      


//Define a macro to compute the minimum
#ifndef MIN1
#define MIN1(x, y) ((x < y) ? 1: 0)
#endif


//Macro
#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

//Macro to compute the maximum
#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif

#ifdef INFINITY
/* INFINITY is supported */
#endif




//Main function starts here
int main() {
    // Define variables to keep track of time
    unsigned long long start = 0;
    unsigned long long finish = 0;

//All the initializations
    FILE *fp, *fw;
    double num1;
    int i, j, k, rnd_num;

    int n =Genes, fs = Samples, k1=K, *label,lab;

    double **data,temp;



    gene_data = (double*) calloc(Genes * Samples, sizeof(double));
    medoids = (double*) calloc(K * Samples, sizeof(double));
    cluster_idx = (int*) calloc(Genes, sizeof(int));



    
    srand(time(0));


    start=clock_now();
    // File open for reading
    fp = fopen("train_rice.txt", "r");
    if (fp == NULL) {
        printf("The requested input file does not exist. \n");
        exit(-1);
    }






   for (i=0;i<Genes;i++){
		for (j=0;j<Samples;j++){
			fscanf(fp,"%lf", &num1);
			*(gene_data+i*Samples+j)=num1;
		}
	}




	fclose(fp);




   finish = clock_now();

 
   printf("Total clockcycles taken to read file serially with %d genes with %d clusters is %lld cycles.\n", Genes, K, (finish-start));
   printf("Total time taken to read file serially with %d genes with %d clusters is %e seconds.\n", Genes, K, (finish-start)/512000000.0f);
 

	free (gene_data);
	free (medoids);
	free(cluster_idx);
}
