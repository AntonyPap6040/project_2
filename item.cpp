//
// Created by antonis on 23/10/2018.
//
#include "item.h"
#include <random>
#include <unordered_map>
#include <chrono>



float t_generator(){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distributionU(0.0, w); //uniformly distributed
    float t;

    t = distributionU(generator);
    return t;
}



vector<float> v_generator(item key) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::normal_distribution<float> distributionN(0.0, 1.0); //normal distribution
    vector<float> V;

    vector<double> p = key.get_point();
    float a;
    for (int i = 0; i < p.size(); i++) {

        a = distributionN(generator);
        V.push_back(a);

    }

    return V;

}

void item::set_V1(vector<float>* V1){
    V=V1;
}

void item::set_t1(float *t1){
    t=t1;
}

void item::set_r1(vector<float>* r1){
   r=r1;
}

void item::set_V(item key){

    V = new vector<float>[k];
    have_memory=true;
    for(int i=0;i<k;i++) {
        V[i]= v_generator(key);
    }

}

void item::the_rand_nums(){
    random_num = new long long[k];
    for(int a=0;a<k;a++){

        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        std::default_random_engine generator(seed);

        random_num[a] = generator();
    }
}

void item::set_t() {
    t = new float[k];
    have_memory=true;
    for(int i=0;i<k;i++) {
        t[i]=t_generator();
    }

}

void item::set_r(item key) {
    r = new vector<float>[k];
    have_memory=true;
    for(int i=0;i<k;i++) {
        r[i]= v_generator(key);
    }

}

vector<float>* item:: getV(){
    return V;
}

float* item::gett(){
    return t;

}

vector<float>* item::getr(){
    return r;

}

int item::hic(const item *key,int i) const { //h functions for cosina method

    //edo upologizoume ta hi gia i apo 1 mexri k
    vector<double> p = key->get_point();

    double m=0.0;
    if(r[i].size()<=p.size()) {
        for (int u = 0; u < r[i].size(); u++) {
            m += r[i][u] * p[u];
        }
    }
    else {
        for (int u = 0;u < p.size(); u++) {

            m += r[i][u] * p[u];
        }
    }

    if(m>=0.0) return 1;
    else return 0;

}

long long item::hie(const item *key,int i) const{   //h functions for euclidean method

    //edo upologizoume ta hi gia i apo 1 mexri k
    vector<double> p = key->get_point();


    double m=0.0;
    if(V[i].size()<=p.size()) {
        for (int u = 0; u < V[i].size(); u++) {
            m += V[i][u] * p[u];
        }
    }
    else {
        for (int u = 0;u < p.size(); u++) {

            m += V[i][u] * p[u];
        }
    }


    long long h = floor((m+t[i])/w);


    return h*random_num[i];
}



item::item(int i,int j,int o){
    isString=false;
    idi=i;
    k=j;
    L=o;
    have_memory=false;
    assigned=false;
    belonging_to_cluster=-1;
    changed=false;
    ids="";


}

item::~item(){


}

void item::freeE(){
    if(have_memory){
        delete[] V;
        delete[] t;
    }
}

void item::freeRand(){
    delete[] random_num;

}


void item::freeC() {
    if(have_memory){
        delete[] r;
    }
}

void item::set(string i){
    isString=true;
    idi=-1;
    ids=i;
}
void item::set(int i){
    isString=false;
    idi=i;
    ids="";
}

string item:: get_type(){

    if(isString) return "string";
    else return "int";
}

int item:: get_idi(){

    return idi;

}

string item:: get_ids(){

    return ids;

}

int item:: get_size(){
    return coordinates.size();
}

vector<double> item:: get_point() const {
    return coordinates;
}


void item::push(int i){

    coordinates.push_back(i);

}

void item::print() {

    for (int j = 0; j < coordinates.size(); j++) {
        cout << coordinates[j];
        cout << " ";
    }
}

//long long item::f()