//
// Created by antonis on 29/10/2018.
//

#ifndef ERGASIA1_HYPER_H
#define ERGASIA1_HYPER_H
#include <iostream>
#include "item.h"
#include <fstream>
#include "item.h"
#include "euclidean.h"
#include "cosine.h"
#include <chrono>
#include <random>
#include <sstream>      // std::stringstream
#include <algorithm>
#include <array>

class hypercube {

public:

    hypercube(){
        count_points=0;
    };

    void insert_vertex(vector<double> vec,item &point){

        f.push_back(vec);
        p.push_back(point);
        count_points++;

    };

    vector<item> get_item(){
        return p;
    };
    vector<vector<double>> get_f(){
        return f;
    };
    void set_size(int si){
        size=si;
    };
    void alloc(int al){
        //f=new vector<double>[al];
    };
    int has_points(){
        return count_points;
    };
private:
    int size;
    int count_points;
    vector<item> p;
    vector<vector<double>> f; //the f value
};

double ham_dist(vector<double> v1, vector<double> v2);

void calculate_f(vector<double> &FV,item* item1,int d2) ;

void calculate_c(vector<double> &FV,item* item1,int d2) ;


int get_vertex(vector<double> f,int d2);


#endif //ERGASIA1_HYPER_H
