#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include"clockcycle.h"

#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif

//Define a macro to compute the minimum
#ifndef MIN1
#define MIN1(x, y) ((x < y) ? 1: 0)
#endif


//Macro
#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

// Calculates Euclidean distance between two points. Using this function throughout
double euclidean_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i ++ ) {
        distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(distance);
}


double manhattan_distance(double *p1, double *p2, int dim) {
      double distance = 0.0;
      for (int i = 0; i < dim; i++) {
      distance += fabs(p1[i] - p2[i]);
      }
      return distance;
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
    if(count==0)
    {
        printf(" %d\n",count);
    }
    return distance_sum / count;
}


// Calculate the average nearest cluster distance for the data point with index 'idx' in 'data'
//  to its nearest cluster
// argument 'k' is the number of clusters
double average_nearest_cluster_distance(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double cluster_distance[k];

    //Initialize their values to zero first
    for (int i = 0; i < k; i++) {
    cluster_distance[i] = 0.0;
}
    int count[k];
       for (int i = 0; i < k; i++) {
    count[i] = 0.0;
}
    for (int i = 0; i < numdata; i ++ ) {
        cluster_distance[labels[i]] += manhattan_distance(data[idx], data[i], dim);
        count[labels[i]] ++ ;
    }
    double res = INFINITY;
    for (int i = 0; i < k; i ++ ) {
        if (i == labels[idx]) continue;
        if (count[i]>0)
            res = MIN(res, cluster_distance[i] / count[i]);
    }

    return res;
}




// Calculate silhouette coefficient for the data point with index 'idx' in 'data'
double silhouette_coefficient_single_point(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double a = average_cluster_distance(data, labels, numdata, dim, idx);
    double b = average_nearest_cluster_distance(data, labels, numdata, k, dim, idx);
    // printf("%lf %lf\n",a,b);
    return (b - a) / MAX(a, b);
}

// Calculate average silhouette coefficient
double silhouette_coefficient(double **data, int *labels, int si, int ei, int k, int dim,int numdata) {
    double totalsc = 0.0;
    for (int i = si; i<ei; i ++ ) {
        totalsc += silhouette_coefficient_single_point(data, labels, numdata, k, dim, i);
    }

    return totalsc/(ei-si+1);
}

int main() { 

//the start index and the ending index for the elements that each process has to consider
    int si,ei;
//mpi initializations, rank determinor
    MPI_Init(NULL, NULL);
    MPI_Status status;
    MPI_Request r_request, s_request;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    double start_time, end_time;
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
//n is the number of elements, fs is the number of features, k is the number of clusters
    int n =7129, fs = 34, k=10, *label,lab;
    double **data,temp;
    
//initialize the spaces for the data (stored as 2d array) and labels
    FILE *fp;
    data = (double**)malloc(n*sizeof(double*));

    for(int i=0;i<n;i++)
    {
        data[i]=(double*)malloc(fs*sizeof(double));
    }
    label = (int*)malloc(n*sizeof(int));


    start_time = clock_now();
// if the rank is 0, i am initializing data and broadcasting it to all the processes
    if(world_rank==0)
 
    {    

      //Proper file handling is necessary to avoid issues of segmentation faults
      
      char filename[100];   
      sprintf(filename, "output_gene_data_%d.txt", k);
    
        //this is the initialization of data with the normalized results
        fp = fopen(filename, "r");
        if(fp!=NULL)
        {
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<fs;j++)
            {
                fscanf(fp,"%lf",&temp);
                data[i][j]=temp;
            }
        }
        }
        else{
           printf("Error opening file\n");
        }

        fclose(fp);



        //this is the initialization of the labels, put the file having cluster index for each data point here
          char filename1[100];   
         sprintf(filename1, "gene_med_ind_%d_clusters.txt", k);
        fp = fopen(filename1, "r");
        if(fp!=NULL)
        {
        for(int i=0;i<n;i++)
        {

            fscanf(fp,"%d",&lab);
            label[i]=lab;
        }
        }
        else {
    printf("Error opening file\n");
       }


        fclose(fp);



        //broadcasting data and the labels. here each data element is sent as a seperate broadcast call
        for(int i=0;i<n;i++)
        {
            MPI_Bcast(data[i], fs, MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(label,n,MPI_INT,0,MPI_COMM_WORLD);
    }
    else
    {
        //recieve data here
        for(int i=0;i<n;i++)
        {
            MPI_Bcast(data[i], fs, MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(label,n,MPI_INT,0,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);


    //determining the range range of points for each process to consider. It is similar to a window

    si = (world_rank)*(n/world_size);

    if(world_rank==(world_size-1))
    {
        ei = n;
    }
    else
    {
        ei = (world_rank+1)*(n/world_size);
    }
    // here we find the sillhouete coefficient for each of the points assigned to the process and use reduce to sum them all
    double a;
    double result;
    a = silhouette_coefficient(data, label, si,ei,k,fs,n);

    printf("proc:%d %lf\n",world_rank,a);
    MPI_Reduce(&a, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    end_time = clock_now();
    
    //here we are finding the average of the results across the processes
    
    if(world_rank==0)
    {
    printf("MPI_Reduce_sum which is same as the silhouette coefficient %lf \n",result/world_size);
    printf("Total time taken to compute the Silhoutte coefficient on %d genes with %d clusters is %e seconds.\n", n, k, (end_time-start_time)/512000000.0f);
    printf("Total clockcycles taken to compute the Silhoutte coefficient on %d genes with %d clusters  is %f cycles.\n",n, k, (end_time-start_time));
    }
    
     MPI_Finalize();
}