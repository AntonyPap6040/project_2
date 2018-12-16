//
// Created by antonis on 30/11/2018.
//

#ifndef ERGASIA2_ASSIGNMENT_H
#define ERGASIA2_ASSIGNMENT_H

#include "initialization.h"

typedef struct cluster_struct{
    vector<double> centroid;
    vector<double>* cluster_data;
    int counter;
    int id;
    double si;
} cluster_struct;

class assignment {

public:

    void out(argClass args,string algo,string opath);

    void m1(argClass args,  initialization init, assignment assign);
    void m2(argClass args,  initialization init, assignment assign);
    void m3(argClass args,  initialization init, assignment assign);
    void m4(argClass args,  initialization init, assignment assign);
    void m5(argClass args,  initialization init, assignment assign);
    void m6(argClass args,  initialization init, assignment assign);
    void m7(argClass args,  initialization init, assignment assign);
    void m8(argClass args,  initialization init, assignment assign);
    void m9(argClass args,  initialization init, assignment assign);
    void m10(argClass args,  initialization init, assignment assign);
    void m11(argClass args,  initialization init, assignment assign);
    void m12(argClass args,  initialization init, assignment assign);


    assignment(int k, int num_of_points ,string m, vector<double>* centroids);
    void clean(int k, int num_of_points, string m,vector<double>* centroids) ;


    void update_centers_kmeans();
    void update_centers_kmedoids();

    void cl_kmeans(argClass args,  initialization init, assignment assign);

    item* Loyds_assignment(vector<double>* points,item * old_clustering);
    item* LSH( argClass args,  initialization init, assignment assign, item * old_clustering);
    item* cube( argClass args,  initialization init, assignment assign,item * old_clustering);


    double silhouette();

    cluster_struct * get_clusters(){
        return clusters;
    }
    void free_assign(){
        free(clusters);
    }


private:

        double sil;
        bool done;
        string metric;
        int centroid_num;
        int num_points;
        cluster_struct* clusters;
        double clustering_time;

};


#endif //ERGASIA2_ASSIGNMENT_H
