
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdbool.h>
#include<time.h>
#include<mpi.h>
#include<string.h>
#include<math.h>
#include "clockcycle.h"
#define Genes 7129    // X Total Number of genes to be given as an input. 
#define Samples 34       // Represents the sample genes
//Change the value of K to obtain the results with different clusters
#define K 3 // Number of clusters
//Initializations
int *cluster_idx;             
double *gene_data;          
double *medoids;   					 // pointer to data  which stores the index of the centroid nearest to each pixel


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


/*INFINITY is a macro constant defined in the <math.h> library, and we use it to perform mathematical comparisons*/
#ifdef INFINITY
/* INFINITY is supported */
#endif


//Finding the closeset medoids
//This function works totally fine
void findclosestmedoids(double *num, double *medoids, int *idx, int rank, int size,int process_job, int si,int ei) {
    int i, j, l, for_i;
    double sum, dist[K], min_dist, local_min_dist;




    // Broadcast the medoids to all processes
    MPI_Bcast(medoids, K * Samples, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Find the closest medoid for each local data point
    for (for_i = 0; for_i <process_job; for_i++) {
        i = si + for_i;
        local_min_dist = INFINITY;

        for (j = 0; j < K; j++) {
            sum = 0;

            for (l = 0; l < Samples; l++) {
             
                sum += fabs(num[i * Samples + l] - medoids[j * Samples + l]);
            }

            dist[j] = sum;

            if (MIN1(dist[j],local_min_dist)) {
                local_min_dist = dist[j];
                idx[i] = j;
            }
        }
    }



    // Reduce the local min distances to find the global min distance
    //The reason why we use
    MPI_Allreduce(&local_min_dist, &min_dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // Update the indices for data points with the global min distance
    for (i = si; i < ei; i++) {
        if (dist[idx[i]] == min_dist) {
            idx[i] = j;
        }
    }
}





//This code is same as the serial code, but the difference is the splitting of the clusters among processes
void computeMedoids(double* gene_data, int* idx, double* medoids, int rank, int size) {
    int i, j, k, l, m;
    double min_distance, distance, sum, temp;

    // using entire gene data instead of sampling
    int sample_size = Genes; 
   
    // Divide the K clusters across the processes
    int medoids_per_proc = K / size;
    int remainder = K % size;
    int start_idx = rank * medoids_per_proc;
    int end_idx = start_idx + medoids_per_proc;
   
    if (rank == size - 1) {
        end_idx += remainder;
    }

    //printf("Process %d: start_idx = %d, end_idx = %d\n", rank, start_idx, end_idx);

    // Seed the random number generator with the rank
    srand(time(0) + rank);

    for (i = start_idx; i < end_idx; i++) {
        min_distance = 1e9;

        // Calculate distances for all points in the cluster
        for (j = 0; j < Genes; j++) {
            if (idx[j] == i) {
                sum = 0.0;
                for (k = 0; k < Genes; k++) {
                    if (idx[k] == i) {
                        distance = 0.0;
                        for (l = 0; l < Samples; l++) {
                            // manhattan distance
                            double diff = *(gene_data + j * Samples + l) - *(gene_data + k * Samples + l);
                            distance += fabs(diff);
                        }
                        sum += distance;
                    }
                }
//I am using the macro i defined for this
                if (MIN1(sum, min_distance)) {
                    min_distance = sum;
                    for (m = 0; m < Samples; m++) {
                        temp = *(gene_data + j * Samples + m);
                        *(medoids + i * Samples + m) = temp;
                    }                    
                }
            }
        }
    }


  /*  // Print the final medoids array from each process takes place here
    printf("Process %d: Medoids = \n", rank);
    for (i = start_idx; i < end_idx; i++) {
        printf("Medoid %d: ", i);
        for (j = 0; j < Samples; j++) {
            printf("%lf ", *(medoids + i * Samples + j));
        }
        printf("\n");
    }*/
/*
//Gene data
    printf("Process %d: gene_data array:\n", rank);
for (i = 0; i < Genes; i++) {
    for (j = 0; j < Samples; j++) {
        printf("%f ", *(gene_data + i * Samples + j));
    }
    printf("\n");
}*/
}







int main(int argc, char *argv[]){

    int myrank, numranks, result;
    int i,j,k;
    int process_job, each_chunk_pos;
    double starttime, endtime;

    MPI_Init(NULL, NULL);
    MPI_Status status;
    MPI_Request r_request, s_request;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int rnd_num;
    int n =Genes, fs = Samples, k1=K, *label,lab;
    MPI_Request request1, request2, request3, request4; 
    double num1;
    MPI_File fh,fh1;
    FILE *fp;
    int si,ei;

    //Determination of each chunk of the input file each process is going to read
    process_job=Genes/world_size;
    if(world_rank==world_size-1) {
        process_job=process_job+Genes%world_size;
    }

    gene_data = (double*) calloc(Genes * Samples, sizeof(double));
    medoids = (double*) calloc(K * Samples, sizeof(double));
    cluster_idx = (int*) calloc(Genes, sizeof(int));

    // This is taken from slide number 11 of MPI file I/O that uses MPI I/O to read from the input file.
    //We also experimented with bin format of the file. However, txt is also faster. 



/*
    result=MPI_File_open(MPI_COMM_WORLD, "input.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(result != MPI_SUCCESS) {printf("Error in opening the file\n"); exit(-1);}
    result=MPI_File_read_at(fh, world_rank*process_job*Samples*sizeof(double),gene_data,process_job*Samples, MPI_INT, &status);
    if(result != MPI_SUCCESS) {printf("Error in reading the file\n"); exit(-1);}
    MPI_File_close(&fh);
   

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&fh);*/


//The gene data is not proper
if (world_rank == 0) {


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


//We are verifying this to see if the I/O has worked fine or not. Uncomment this to see the gene data
 /*   printf("Gene Data:\n");
    for (i = 0; i <Genes; i++) {
        printf("Gene %d: ", i + 1);
        for (j = 0; j < Samples; j++) {
            printf("%lf ", *(gene_data + i * Samples + j));
        }
        printf("\n");
    }
*/

   

}








  	if (world_rank==0)
	{	
    
            starttime = clock_now();
              int lower =0;
              int upper =process_job-1;
              srand(time(0));

  
  
    
//Random initialize of medoids takes place here
		for (int i = 0; i < K; i++) {

			int rnd_num = (rand()%(upper-lower + 1)) + lower;
	
		
			for (j=0;j<Samples;j++){ 
        		*(medoids+i*Samples+j) = *(gene_data+rnd_num*Samples+j); 
        	} 
      
    }
  
  }

   
 //This is for the findclosestmedoids function

    si = (world_rank)*(Genes/world_size);

    if(world_rank==(world_size-1))
    {
        ei = Genes;
    }
    else
    {
        ei = (world_rank+1)*(Genes/world_size);
    }

    
 /*We run the findclosestmedoids for 10 iterations. 
 The broadcast to all the processes and also the computation of final result using MPI_Allreduce is done in
 the function itself. Therefore we do not do it here. But we meed to do all gather to gather the results from all processes after
 the computeMed function */



	//MPI_Bcast(medoids, K*Samples, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for (i=0;i<10;i++){

        findclosestmedoids((double *)gene_data, (double *)medoids, &cluster_idx[0],
        world_rank,process_job,world_size,si,ei);

        computeMedoids((double *)gene_data, &cluster_idx[0], (double *)medoids,world_rank, world_size);
        MPI_Barrier(MPI_COMM_WORLD);


        //All gather must be done here//

         MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, medoids, K * Samples /world_size, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int a=0; a<K;a++){
        
          for (int b=0;b<Samples;b++){

            *(medoids+a*Samples+b)=*(medoids+a*Samples+b)/world_size;
          }
        }
      

      }
      
    // Print the medoids computed by the root process
    if (world_rank == 0) {

        endtime = clock_now();

        printf("Total time taken to run K-Medoids on %d genes with %d clusters is %e seconds.\n", Genes, K, (endtime-starttime)/512000000.0f);
        printf("Total clockcycles taken to run K-Medoids on %d genes with %d clusters  is %f cycles.\n", Genes, K, (endtime-starttime));
        printf("Medoids:\n");
        for (i = 0; i < K; i++) {
            printf("Medoid %d: ", i + 1);
            for (j = 0; j < Samples; j++) {
                printf("%lf ", *(medoids + i * Samples + j));
            }
            printf("\n");
        }


   /*       
    for (i = 0; i < Genes; i++) {
    
        printf("%d\n",cluster_idx[i]);

    }*/
       
    }

    // Free allocated memory
    if (world_rank == 0) {
        free(gene_data);
        free(cluster_idx);
        free(medoids);
    }

    MPI_Finalize();
    

   
    return 0;
}


