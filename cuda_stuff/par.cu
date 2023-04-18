#include <stdio.h>
#include <stdlib.h>

#define K 10
#define n 7129
#define nf 34


#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif


//here i am doing a reduction to find an array that has the minimum element and minimum element index for each block
__global__ void find_min(double *m_arr, int *mind, double *mval)
{
    // shared data for each block
    __shared__ double sdata[1024];
    __shared__ int sind[1024];
    int tid = threadIdx.x;
    size_t idx = threadIdx.x + blockDim.x * blockIdx.x;
    // putting data into shared memory
    if(idx<n)
    {
        sdata[tid]=m_arr[idx];
        sind[tid]=idx;
    }
    else
    {
        sdata[tid]=INFINITY;
        sind[tid]=-1;
    }

    __syncthreads();


    // reduction part 
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        __syncthreads();
        if (tid < s) // parallel sweep reduction
            if(sdata[tid]>=sdata[tid + s])
            {
                sdata[tid] = MIN(sdata[tid],sdata[tid + s]);
                sind[tid] = sind[tid+s];
            }
    }

    //if the thread id is 0, I am writing it to the min value and min index array
    if (tid == 0){mval[blockIdx.x]=sdata[0]; mind[blockIdx.x]=sind[0];}
}



//here i am making an array for average distances where each index i correcsponds to the average distance of point i
//I have these conditions while filling the array
//case 1: if two points have same cluster index, i am adding the distance to the point in consideration's avg distance
//case 2: if 2 points have different cluster index, i am giving the value of infinity
__global__ void find_avg_dist(double *data, int k, int *labels, double *avgs)
{
    int ind = blockIdx.x*blockDim.x + threadIdx.x;
    double as=0;
    int c = 0;

    if (labels[ind] != k)
    {
        avgs[ind] = INFINITY;
    }

    else
    {

        for (int i = 0; i < n; i++)
        {
            if (i != ind and labels[i] == k)
            {
                c += 1;
                double ms = 0;

                for (int j = 0; j < nf; j++)
                {
                    ms += fabs(*(data + nf * ind + j) - *(data + nf * i + j));
                }
                as += ms;
            }
            else
            {
                avgs[ind] = INFINITY;
            }
        }
        if (c == 0)
        {
            avgs[ind] = INFINITY;
        }
        else
        {
            avgs[ind] = as / c;
        }
    }
}



__global__ void findclosestmedoids_ker(double *data, double *medoids, int *med_labels,int si,int ei)
{
    int max_ind;
    int i = si+( blockIdx.x*blockDim.x + threadIdx.x);

    if(i<=ei)
    {
        double md=INFINITY;

        for(int j=0;j<K;j++)
        {
            double s=0;
            for(int k=0;k<nf;k++)
                s+=abs(*(data+i*nf+k)-*(medoids + j*nf+k));
            
            if(s<md)
            {
                md=s;
                max_ind = j;
            }
        }

        med_labels[i]=max_ind;
    }
}


extern "C" void computeMedoids(double* data, int* labels, double* medoids, int rank, int size) 
{
    int blockSize = 1024;
    int nblocks = (n+blockSize-1)/blockSize;
    int cE,cudaDeviceCount;

    if( (cE = cudaGetDeviceCount( &cudaDeviceCount)) != cudaSuccess )
    {
    printf(" Unable to determine cuda device count, error is %d, count is %d\n",
    cE, cudaDeviceCount );
    exit(-1);
    }
    if( (cE = cudaSetDevice( rank % cudaDeviceCount )) != cudaSuccess )
    {
    printf(" Unable to have rank %d set to cuda device %d, error is %d \n",
    rank, (rank % cudaDeviceCount), cE);
    exit(-1);
    }

    // Divide the K clusters across the processes
    int medoids_per_proc = K / size;
    int remainder = K % size;
    int si = rank * medoids_per_proc;
    int ei = si + medoids_per_proc;
   
    if (rank == size - 1) {
        ei += remainder;
    }

    double *avgs = (double*)malloc(n*sizeof(double));
    int *mind = (int*)malloc(nblocks*sizeof(int));
    double *mval =  (double*)malloc(nblocks*sizeof(double));


    for(int i=si;i<ei;i++)
    {
        // here for each cluster I am first computing the avg distance array and then I am finding the min index
        find_avg_dist<<<nblocks,blockSize>>>(data, i, labels,avgs);
        cudaDeviceSynchronize();
        find_min<<<nblocks,blockSize>>>(avgs, mind,mval);
        cudaDeviceSynchronize();

        // finding min index from the blocksize array 
        double mv= INFINITY;
        int mi;
        for(int j=0;j<nblocks;j++)
        {
            if(mval[j]<mv)
            {
                mv=mval[j];
                mi=mind[j];
            }
        }

        //assigning the medoid for the particular cluster index
        for(int j=0;j<nf;j++)
            *(medoids + nf*i+j)=*(data + nf*mi+j);
        
    }


        printf("done\n");






}

extern "C"  void findclosestmedoids(double *data, double *medoids, int *idx, int rank, int size,int process_job, int si,int ei)
{
    int blockSize = 1024;
    int size1 = (ei-si) + 1;
    int nblocks = (size1+blockSize-1)/blockSize;
    int cE,cudaDeviceCount;

    printf("%d %d\n",si,ei);
    if( (cE = cudaGetDeviceCount( &cudaDeviceCount)) != cudaSuccess )
    {
    printf(" Unable to determine cuda device count, error is %d, count is %d\n",
    cE, cudaDeviceCount );
    exit(-1);
    }
    if( (cE = cudaSetDevice( rank % cudaDeviceCount )) != cudaSuccess )
    {
    printf(" Unable to have rank %d set to cuda device %d, error is %d \n",
    rank, (rank % cudaDeviceCount), cE);
    exit(-1);
    }

    findclosestmedoids_ker<<<nblocks,blockSize>>>(data, medoids, idx,si,ei);

    cudaDeviceSynchronize();





}


