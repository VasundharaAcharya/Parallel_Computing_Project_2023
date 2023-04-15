#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "clockcycle.h"

#define Genes 7000   // X Total Number of genes to be given as an input. 
#define Samples 10          // Represents the sample genes
//Change the value of K to obtain the results with different clusters
#define K 10 // Number of clusters
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

// Specification of the distance metrics for the task of gene clustering.
// We use Manhattan Distance here
double manhattan_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i ++ ) {
        distance += abs(p1[i] - p2[i]) + abs(p1[i] - p2[i]);
    }
    return distance;
}




// This function is used to compute the closest medoids 

void findclosestmedoids(double *gene_data, double *medoids, int* idx) {
    int i, j, l;
    double sum, dist[K], min_dist;

    for (i = 0; i < Genes; i++) {
        min_dist = INFINITY;

        for (j = 0; j < K; j++) {
            sum = 0;

            for (l = 0; l < Samples; l++) {
              //Manhattan distance
                sum += fabs(gene_data[i * Samples + l] - medoids[j * Samples + l]);
            }

            dist[j] = sum;

            if (MIN1(dist[j],min_dist)) {
                min_dist = dist[j];
                idx[i] = j;
            }
        }
    }
}




void computeMedoids(double* gene_data, int* idx, double* medoids) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;

    for (i = 0; i < K; i++) {
        min_distance =1e9;
        // calculate distances for all points in the cluster
        for (j = 0; j < Genes; j++) {
            if (idx[j] == i) {
                sum = 0.0;
                for (k = 0; k < Genes; k++) {
                    if (idx[k] == i) {
                        distance = 0.0;
                        for (l = 0; l < Samples; l++) {
                            //manhattan distance
                            double diff = *(gene_data + j * Samples + l) - *(gene_data + k * Samples + l);
                            distance += fabs(diff);
                        }
                        sum += distance;
                    }
                }
                if (MIN1(sum,min_distance)) {
                    min_distance = sum;
                    for (m = 0; m < Samples; m++) {
                        temp = *(gene_data + j * Samples + m);
                        *(medoids + i * Samples + m) = temp;
                    }
                }
            }
        }
       
    }
}




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


    
    // File open for reading
    fp = fopen("training.txt", "r");
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


//Start the clockcycle

	start = clock_now();


	//Creation of random number and start with random medoids
	for (i = 0; i < K; i++) {

			rnd_num = (rand()%n);
		
			for (j=0;j<Samples;j++){ 
        		*(medoids+i*Samples+j) = *(gene_data+rnd_num*Samples+j); 
        	} 
      
    }
  

  

      i = 0;
      //Change the number of iterations here if required
      while (i < 1000) {
        findclosestmedoids((double*)gene_data, (double *)medoids, &cluster_idx[0]);
      


        computeMedoids((double *)gene_data, &cluster_idx[0], (double *)medoids);
          i++;
      }

        for(i=0; i<K;i++){
            for(j=0; j<Samples;j++){

              printf("Medoids are  %lf  ",*(medoids+i*Samples+j));

            }
          printf("\n");
        }
	


  printf("kmedoids completed with clustering of the gene samples \n");








    // Write clustered gene data to file
    int flag1=0;
    char filename3[100];
    sprintf(filename3, "output_gene_data_%d_euclidean.txt", K);
    fw=fopen(filename3,"w");
    for(i=0; i<Genes;i++){
        for(j=0; j<Samples;j++){
          flag1= fprintf(fw,"%lf  ",*(gene_data+i*Samples+j));
        }
        fprintf(fw, "\n");
    }
    fclose(fw);


    if(flag1)
    {

      printf("Successfully wrote the clustered gene data to the file \n");
    }



    // Write gene medoids to file
    int flag=0;
    char filename2[100];
    sprintf(filename2, "output_gene_medoids_%d_clusters_euclidean.txt", K);
    fw=fopen(filename2,"w");
    for(i=0; i<K;i++){
        for(j=0; j<Samples;j++){
            flag=fprintf(fw,"%lf  ",*(medoids+i*Samples+j));
        }
        fprintf(fw, "\n");
    }
    fclose(fw);




    if(flag)
    {

      printf("Successfully wrote the gene medoids to the file \n");
    }


    //Write the cluster assignments to the file



    int flag2=0;
    char filename[100];
    sprintf(filename, "gene_med_ind_%d_clusters_euclidean.txt", K);
    fw = fopen(filename,"w");
    for (i = 0; i < Genes; i++) {
    
        flag2=fprintf(fw,"%d  ",cluster_idx[i]);

    }

  fclose(fw);
  if(flag2)
  {

    printf("Successfully wrote the cluster assignments to the file \n");
  }

   finish = clock_now();

 
   printf("Total clockcycles taken to run K-Medoids on %d genes with %d clusters is %lld cycles.\n", Genes, K, (finish-start));
   printf("Total time taken to run K-Medoids on %d genes with %d clusters is %e seconds.\n", Genes, K, (finish-start)/512000000.0f);
 

	free (gene_data);
	free (medoids);
	free(cluster_idx);
}
