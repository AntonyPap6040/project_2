//
// Created by antonis on 28/11/2018.
//

#ifndef ERGASIA2_INITIALIZATION_H
#define ERGASIA2_INITIALIZATION_H

#include <random>
#include <chrono>
#include <fstream>
#include "item.h"
#include "euclidean.h"
#include "cosine.h"
#include <sstream>      // std::stringstream
#include <algorithm>
#include <unordered_set>
#include "initialization.h"
#include "argClass.h"


class initialization {
public:

    initialization(string m,int a,int k){
        metric=m;
        data_num = a;
        cluster_centroid_points = new vector<double>[k];

    }

    //----INITIALIZATION FUNCTIONS BELOW----

    //---SIMPLY RANDOM---
    void random_k_points(int k,int point_length,vector<double>* data_points);

    //---K-means++ ---
    void k_means_pp(int k,int point_length,vector<double>* data_points);

    vector<double>* get_centroids(){
        return cluster_centroid_points;
    }

    void calculate_min_max(int k);


    double get_min(){
        return min_dist_centers;
    }

    double get_max(){
        return max_dist_centers;
    }

private:
    vector<double>* cluster_centroid_points;
    string metric;
    double min_dist_centers; //minimum distance between centers
    double max_dist_centers; //minimum distance between centers

    int data_num;
};

bool double_equals(double a, double b, double epsilon = 0.000000000000001);
bool double_less(double a, double b, double epsilon = 0.000000000000001);
bool double_greater(double a, double b, double epsilon = 0.000000000000001);

#endif //ERGASIA2_INITIALIZATION_H
