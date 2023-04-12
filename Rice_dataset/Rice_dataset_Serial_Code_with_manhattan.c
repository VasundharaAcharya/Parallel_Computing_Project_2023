

/* This is the rice dataset. So we change the number of genes to 10k and the samples to 35. As the dataset size increases in this dataset, we 
use the precomputed pairwise distance to achieve the speedup in this case */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "clockcycle.h"

#define Genes 10000     // X Total Number of genes to be given as an input. 
#define Samples 35     // Represents the sample genes
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

//Especially for this dataset, as the size increased, after talking to the professor, we decided to use the precomputed distances

void computeMedoids(double* gene_data, int* idx, double* medoids) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;
    
    // This is the only extra step we do for this dataset, we precompute the distances to speedup as the dataset is huge
    double** pairwise_distances = (double**)malloc(Genes * sizeof(double*));
    for (i = 0; i <Genes; i++) {
        pairwise_distances[i] = (double*)malloc(Genes * sizeof(double));
    }
    for (i = 0; i < Genes; i++) {
        for (j = i + 1; j <Genes; j++) {
            sum = 0.0;
            for (k = 0; k < Samples; k++) {
                sum += fabs(*(gene_data + i * Samples + k) - *(gene_data + j * Samples + k));
            }
            pairwise_distances[i][j] = pairwise_distances[j][i] = sum;
        }
    }


//This part is same as the previous codes
  // set the sample size as needed here
    int sample_size = 10; 
    for (i = 0; i < K; i++) {
        //Initialize to some large numbers
        min_distance =1e9;
        // randomly sample a subset of points from the cluster
        int* sample = (int*)malloc(sample_size * sizeof(int));
        // I used this code to use the index only once

        for (j = 0; j < sample_size; j++) {
            sample[j] = -1;
            while (sample[j] == -1 || idx[sample[j]] != i) {
                //Sample the genes randomly
                sample[j] = rand() % Genes;
            }
        }
        // calculate distances only for the sampled points
        // I have used the manhattan distance, we can experiment with others as well
        for (j = 0; j < sample_size; j++) {
            sum = 0.0;
            for (k = 0; k < sample_size; k++) {
                distance = 0.0;
                for (l = 0; l < Samples; l++) {
                  //manhattan distance
                    distance += fabs(*(gene_data + sample[j] * Samples + l) - *(gene_data + sample[k] * Samples + l));
                }
                sum += distance;
            }
            if (MIN1(sum,min_distance)) {
                min_distance = sum;
                for (m = 0; m < Samples; m++) {
                    temp = *(gene_data + sample[j] * Samples + m);
                    *(medoids + i * Samples + m) = temp;
                }
            }
        }
        free(sample);
    }



    // free memory allocated for pairwise distances
    for (i = 0; i <Genes; i++) {
        free(pairwise_distances[i]);
    }
    
}





// Calculate the average distance for the data point with index 'idx' in 'data'
//  within its cluster
double average_cluster_distance(double **data, int *labels, int numdata, int dim, int idx) {
    double distance_sum = 0.0;
    int count = 0;
    for (int i = 0; i < numdata; i ++ ) {
        if (labels[i] != labels[idx]) continue;
        distance_sum += manhattan_distance(data[idx], data[i], dim);
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
        cluster_distance[labels[i]] +=manhattan_distance(data[idx], data[i], dim);
        count[labels[i]] ++ ;
    }
    double res = INFINITY;
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
      while (i < 5) {
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
    sprintf(filename3, "output_gene_data_%d.txt", K);
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
    sprintf(filename2, "output_gene_medoids_%d_clusters.txt", K);
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
    sprintf(filename, "gene_med_ind_%d_clusters.txt", K);
    fw = fopen(filename,"w");
    for (i = 0; i < Genes; i++) {
    
        flag2=fprintf(fw,"%d  ",cluster_idx[i]);

    }

  fclose(fw);
  if(flag2)
  {

    printf("Successfully wrote the cluster assignments to the file \n");
  }




  // Start of silhoutte coefficient. We are using this as we cannot plot these huge data.
  //To validate the results of the clustering algorithm, we use this. 




    data = (double**)malloc(n*sizeof(double*));

    for(int i=0;i<n;i++)
    {
        data[i]=(double*)malloc(fs*sizeof(double));
    }

    label = (int*)malloc(n*sizeof(int));




    fp = fopen(filename3, "r");
    if(fp == NULL)
    {
        printf("Error, the requested file does not exist");   
        exit(1);             
    }
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<fs;j++)
        {
            fscanf(fp,"%lf",&temp);
            data[i][j]=temp;
        }
    }

    fclose(fp);





   fp = fopen(filename, "r");
   if(fp == NULL)
   {
      printf("Error, the requested file does not exist");   
      exit(1);             
   }
    for(int i=0;i<n;i++)
    {

        fscanf(fp,"%d",&lab);
        label[i]=lab;
    }
    fclose(fp);

    printf("sillhouette coefficient:%lf\n",silhouette_coefficient(data,label,n,k1,fs));



// end of sill coeff

	// End the timer here
  finish = clock_now();
  printf("Total clockcycles taken to run K-Medoids and compute silhouette score on %d genes in rice dataset with %d clusters is %lld cycles.\n", Genes, K, (finish-start));
  printf("Total time taken to run K-Medoids and compute silhouette score on %d genes in rice dataset with %d clusters is %e seconds.\n", Genes, K, (finish-start)/512000000.0f);

	free (gene_data);
	free (medoids);
	free(cluster_idx);
}
