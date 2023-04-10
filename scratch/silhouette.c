/**
 * compile from terminal/shell as follows:
 * 
 *   bash$ gcc -Wall -Werror silhouette.c
 * 
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

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
    double res = 1e9;
    for (int i = 0; i < k; i ++ ) {
        if (i == labels[idx]) continue;
        res = MIN(res, cluster_distance[i] / count[i]);
    }
    return res;
}


#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif

// Calculate silhouette coefficient for the data point with index 'idx' in 'data'
double silhouette_coefficient(double **data, int *labels, int numdata, int k, int dim, int idx) {
    double a = average_cluster_distance(data, labels, numdata, dim, idx);
    double b = average_nearest_cluster_distance(data, labels, numdata, k, dim, idx);
    return (b - a) / MAX(a, b);
}

int main() { return EXIT_SUCCESS; };