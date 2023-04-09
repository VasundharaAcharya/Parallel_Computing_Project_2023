double euclidean_distance(double *p1, double *p2, int dim) {
    double distance = 0.0;
    for (int i = 0; i < dim; i++) {
        distance += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return sqrt(distance);
}


double average_cluster_distance(double **data, int *labels, int n, int dim, int cluster) {
    double distance_sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (labels[i] != cluster) continue;
        for (int j = 0; j < n; j++) {
            if (labels[j] != cluster || i == j) continue;
            distance_sum += euclidean_distance(data[i], data[j], dim);
            count++;
        }
    }
    return distance_sum / count;
}





double average_nearest_cluster_distance(double **data, int *labels, int n, int dim, int cluster) {
    double distance_sum = 0.0;
    int count = 0;
    for (int i = 0; i < n; i++) {
        if (labels[i] == cluster) continue;
        double distance = average_cluster_distance(data, labels, n, dim, labels[i]);
        distance_sum += distance;
        count++;
    }
    return distance_sum / count;
}


#ifndef MAX
#define MAX(x, y) ((x < y) ? y : x)
#endif

double silhouette_coefficient(double **data, int *labels, int n, int dim, int point) {
    int cluster = labels[point];
    double a = average_cluster_distance(data, labels, n, dim, cluster);
    double b = average_nearest_cluster_distance(data, labels, n, dim, cluster);
    double s;
    s = (b - a) / MAX(a, b);
    return s;
}
