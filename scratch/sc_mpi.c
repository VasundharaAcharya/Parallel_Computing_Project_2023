#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
// #ifndef CLOCKCYCLE_H
// #define CLOCKCYCLE_H


// double start_time,end_time;
// uint64_t clock_now()
// {
//   unsigned int tbl, tbu0, tbu1;

//   do {
//     __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
//     __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
//     __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
//   } while (tbu0 != tbu1);
//   return (((uint64_t)tbu0) << 32) | tbl;
// }
// #endif // CLOCKCYCLE_H



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
    if(count==0)
    {
        printf("fdsddddddddd %d\n",count);
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
        if (count[i]>0)
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


    int si,ei;

  MPI_Init(NULL, NULL);
    MPI_Status status;
     MPI_Request r_request, s_request;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int n =7129, fs = 34, k=15, *label,lab;

    double **data,temp;

    FILE *fp;
    data = (double**)malloc(n*sizeof(double*));

    for(int i=0;i<n;i++)
    {
        data[i]=(double*)malloc(fs*sizeof(double));
    }

    label = (int*)malloc(n*sizeof(int));


    //   start_time = clock_now();

    if(world_rank==0)
    {
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

        for(int i=0;i<n;i++)
        {
            MPI_Bcast(data[i], fs, MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(label,n,MPI_INT,0,MPI_COMM_WORLD);
    // printf("sillhouette coefficient:%lf\n",silhouette_coefficient(data,label,n,k,fs));
    }
    else
    {

        for(int i=0;i<n;i++)
        {
            MPI_Bcast(data[i], fs, MPI_DOUBLE,0,MPI_COMM_WORLD);
        }
        MPI_Bcast(label,n,MPI_INT,0,MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);




    si = (world_rank)*(n/world_size);

    if(world_rank==(world_size-1))
    {
        ei = n;
    }
    else
    {
        ei = (world_rank+1)*(n/world_size);
    }
    double a;
    double result;
    a = silhouette_coefficient(data, label, si,ei,k,fs,n);

    printf("proc:%d %lf\n",world_rank,a);
    MPI_Reduce(&a, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//   end_time = clock_now();

    
    if(world_rank==0)
    {
    printf("MPI_Reduce_sum %lf \n",result/world_size);
    }
    
     MPI_Finalize();




    
};
