//
// Created by antonis on 23/10/2018.
//
#ifndef ERGASIA1_ITEM_H
#define ERGASIA1_ITEM_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>



const int w=450; //window
const double M1 = pow(2.0,32.0)-5;

using namespace std;


class item {

public:
    item(int i=0,int j=4,int o=5);
    ~item();

    int k,L;
    string get_type();
    void set(string i);
    void set(int i);

    int get_idi();
    string get_ids();
    void push(int i);
    void print();
    int get_size();
    vector<double>  get_point() const;
    int TableSize;
    void freeE(); //freeing memory for euclidean method
    void freeC(); //freeing memory for cosine method

    void set_V(item key);
    void set_t();
    void set_V1(vector<float>* V1);
    void set_t1(float* t1);
    vector<float>* getV();
    float* gett();

    void set_r(item key);
    void set_r1(vector<float>* r1);
    vector<float>* getr();

    void the_rand_nums();
    void freeRand();

    void set_coord(vector<double> v){

        coordinates = v;

    }



    long long hie(const item *key,int i) const;

    string g() const {

        string d1;
        int d;
        for (int i = 0; i < k; i++) {
            d += hie(this,i);
            d1.append(to_string(d)); //concatenating all the hi for i in 0,k to get the g func

        }
        return d1;

    };

    int hic(const item *key,int i) const;

    string G() const {

        string d;
        for (int i = 0; i < k; i++) {
            d.append(to_string(hic(this,i))); //concatenating all the hi for i in 0,k

        }
        return d;

    };


    string get_g(){
        string d1;
        int d;
        for (int i = 0; i < k; i++) {
            d += hie(this,i);
            d1.append(to_string(d)); //concatenating all the hi for i in 0,k to get the g func

        }
        g1=d1;
        return g1;
    };

    double memory(){
        double mem=0;
        mem+=g1.size()+2*sizeof(bool)+sizeof(int)+ids.size()+coordinates.size()*sizeof(double)+2*sizeof(float)*5*k+sizeof(float)*k+k*sizeof(long long);
        return mem;
    }

    void cluster_assign(int cl,double d){
        belonging_to_cluster = cl;
        clust_dist=d;
        assigned = true;
    }
    bool is_assigned(){
        return assigned;
    }
    int get_cluster_ball(){
        return belonging_to_cluster;
    }
    double get_clust_dist(){
        return clust_dist;
    }

    void change(){
        changed = true;
    }

    bool is_changed(){
        return changed;
    }

    bool is_equal(item p){

        for(int i=0; i<p.get_size();i++){
            if(p.get_point()[i]!=coordinates[i]) return false;
        }

        return true;

    }
private:
    string g1;
    bool have_memory;
    bool isString;
    int idi;
    string ids;
    vector<double> coordinates;
    vector<float> * V;
    float *t;
    vector<float> * r;
    long long* random_num;
    bool assigned;
    int belonging_to_cluster;
    double clust_dist;
    bool changed;
};


vector<float> v_generator(item key) ;
float r_generator() ;
float t_generator();


#endif
