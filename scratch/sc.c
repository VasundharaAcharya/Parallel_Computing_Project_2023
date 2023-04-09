#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>



double manhattan_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i++) {
        distance += abs(p1[i] - p2[i]);
    }
    return distance;
}


double average_cluster_distance(double **data, int i, int *labels, int n, int dim, int cluster) {
    double distance_sum = 0.0;
    int count = 0;
            for (int j = 0; j < n; j++) {
                if (labels[j] == cluster && i != j) {
                    distance_sum += manhattan_distance(data[i], data[j], dim);
                    count++;
                }
            }
        
    return distance_sum / count;
}





double average_nearest_cluster_distance(double **data, int i, int *labels, int n, int dim) {
    double md = INFINITY;
    int count = 0;
    for (int j = 0; j < 15; j++) {
        if (labels[i] != j) {

            double distance = average_cluster_distance(data,i, labels, n, dim, j);
            printf("%d\n",j);
            if(distance<md)
                md=distance;
        }
    }
    return md;
}



double silhouette_coefficient(double **data, int *labels, int n, int dim, int point) {
    int cluster = labels[point];
    double a = average_cluster_distance(data, point, labels, n, dim, cluster);
    printf("hi\n");
    double b = average_nearest_cluster_distance(data, point, labels, n, dim);
    printf("hi\n");
    double s;
    if(a>=b)
        s = (b - a) / a;
    else
        s = (b - a) / b;
    return s;
}



int main()
{
    int X = 7129;
    int Y = 34;
    int **labels = (int**)malloc(8*sizeof(int*));
    int lab;
    double ***data,temp;

    data = (double***)malloc(8*sizeof(double**));
    //initializing the labels
    for(int i=0;i<8;i++)
    {
        data[i] = (double**)malloc(1000*sizeof(double*));
    }
    for(int i=0;i<8;i++)
    {
        labels[i] = (int*)malloc(1000*sizeof(int));
    }




    for(int i=0;i<8;i++)
    {

    for (int j=0;j<1000;j++)
    {


        data[i][j]=(double*)malloc(Y*sizeof(double));
    }
    }



    FILE *fp;
    fp = fopen("gene_med_ind.txt","r");
    for (int i=0;i<8;i++)
    {
    for (int j=0;j<1000;j++){
        fscanf(fp,"%d",&lab);
        labels[i][j]=lab;
    }
    }

    fp=fopen("training.txt","r");

    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 1000; j++)
        {
            for (int k = 0; k < Y; k++)
            {
                if (fgetc(fp) != EOF)
                {
                    fscanf(fp, "%lf", &temp);
                    data[i][j][k] = temp;
                }
            }
        }
    }

    //find avg sc
    double su=0;
    int co=0;
    for(int i=0;i<120;i++)
    {
        su += silhouette_coefficient(data[0],labels[0],1000,34,i);
        co+=1;
    }
    printf("%lf",su/co);
}