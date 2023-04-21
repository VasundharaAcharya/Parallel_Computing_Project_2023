#include <stdio.h>
#include <stdlib.h>

extern "C" {
    void computeMedoids(double* data, int* labels, double* medoids, int rank, int size);
    void findclosestmedoids(double *data, double *medoids, int *idx , int rank, int process_job,int size, int si,int ei);
}

#define K 10
#define n 7129
#define nf 34
#define REDUNDANT_SIZE 100


#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif


//here i am doing a reduction to find an array that has the minimum element and minimum element index for each block
__global__ void find_min(double *m_arr, int *mind, double *mval)
{
    // shared data for each block
    __shared__ double sdata[1024 + 10];
    __shared__ int sind[1024 + 10];
    int tid = threadIdx.x;
    size_t idx = threadIdx.x + blockDim.x * blockIdx.x;
    // putting data into shared memory
    if(idx<n)
    {
        sdata[tid]=m_arr[idx];
        if(m_arr[idx]==INFINITY)
            sind[tid]=-1;
        else
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
    int ind = blockIdx.x * blockDim.x + threadIdx.x;
    double as = 0;
    int c = 0;

    if (ind < n)
    {

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
    else
    {
        avgs[ind]=INFINITY;
    }
}

__global__ void findclosestmedoids_kernel(double *data, double *medoids, int *ids,int si,int ei)
{
    int i = si+( blockIdx.x*blockDim.x + threadIdx.x);
    if(i<ei+( blockIdx.x*blockDim.x + threadIdx.x))
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
                *(ids+i)=j;
            }
        }
    }
}


void computeMedoids(double* data, int* labels, double* medoids, int rank, int size) 
{
    //here i am initializing the 
    int blockSize = 1024;
    int nblocks = (n+blockSize-1)/blockSize;
    int cudaDeviceCount;
    cudaError_t cE;


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
    int medoids_per_proc = (K + size - 1) / size; // 2 = 15 / 6
    int remainder = K % size; // 4 = 10 % 6
    int si = rank * medoids_per_proc;
    int ei = si + medoids_per_proc;
    if (rank >= remainder) 
    {
        int diff = rank - remainder;
        si -= diff;
        ei = si + K / size;
    }
   

    #ifdef DEBUG_CUDA
    printf("  rank %d: computeMedoids(): cudaDeviceCount = %d,  medoids_per_proc = %d,  si = %d,  ei = %d\n", 
        rank, cudaDeviceCount, medoids_per_proc, si, ei);
    #endif

    double *avgs, *mval;
    int *mind;

    cudaMallocManaged(&avgs, nblocks*blockSize*sizeof(double) + REDUNDANT_SIZE*sizeof(double));
    cudaMallocManaged(&mind, nblocks*sizeof(int) + REDUNDANT_SIZE*sizeof(int));
    cudaMallocManaged(&mval, nblocks*sizeof(double) + REDUNDANT_SIZE*sizeof(double));

        // double *avgs = (double*)malloc(nblocks*blockSize*sizeof(double));
        // int *mind = (int*)malloc(nblocks*sizeof(int));
        // double *mval =  (double*)malloc(nblocks*sizeof(double));

    for(int i=si;i<ei;i++)
    {
        // here for each cluster I am first computing the avg distance array and then I am finding the min index

        find_avg_dist<<<nblocks,blockSize>>>(data, i, labels,avgs);
    #ifdef DEBUG_CUDA
    printf("    rank %d: computeMedoids(): i = %d: finish find_avg_dist()\n", 
        rank, i);
    #endif
        cudaDeviceSynchronize();



        find_min<<<nblocks,blockSize>>>(avgs, mind, mval);
    #ifdef DEBUG_CUDA
    printf("    rank %d: computeMedoids(): i = %d: finish find_min()\n", 
        rank, i);
    #endif
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

        
    #ifdef DEBUG_CUDA
    printf("    rank %d: computeMedoids(): i = %d: mv = %.03lf, mi = %d\n", 
        rank, i, mv, mi);
    #endif

        // printf("index of medoid:%d cluster index:%d\n",mi,i);

        // printf("mval:\n");
        // for(int h=0;h<nblocks;h++)
        // {
        //     printf("%lf ",mval[h]);
        // }
        // printf("\n");


        // printf("mind:\n");
        // for(int h=0;h<nblocks;h++)
        // {
        //     printf("%d ",mind[h]);
        // }
        // printf("\n");

        //assigning the medoid for the particular cluster index
        if(mi!=-1)
        {
        for(int j=0;j<nf;j++)
            *(medoids + nf*i+j)=*(data + nf*mi+j);
        }

        

        
    }

        cudaFree(avgs);
        cudaFree(mind);
        cudaFree(mval);

}

void findclosestmedoids(double *data, double *medoids, int *idx , int rank, int process_job,int size, int si,int ei)
{
    int blockSize = 1024;
    // printf("%d %d\n",si,ei);
    int size1 = (ei-si) + 1;
    int nblocks = (process_job+blockSize-1)/blockSize;
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


    findclosestmedoids_kernel<<<nblocks,blockSize>>>(data, medoids, idx,si,ei);
    cudaDeviceSynchronize();


    // printf("done\n");







}