//
// Created by antonis on 28/11/2018.
//

#ifndef ERGASIA2_ARGCLASS_H
#define ERGASIA2_ARGCLASS_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

#include "item.h"

using namespace std;

class argClass {
public:

    argClass(int a=5, int b=5, int c=5 ){
        clusters=a;
        L=b;
        hash_func_num=c;
    };



    //arg manipulation
    int give_args(int argc, char* argv[]){

        if(argc<9){
            cout<<"wrong number of args please try again";
            return -1;
        }

        //examining call (flags)
        for (int i=0 ; i<argc ; i++){
            string arg = argv[i];
            if (arg== "-i"){
                inputFile= argv[i+1];
                init_data(inputFile);
                set_data_values(inputFile);
            }
            else if (arg== "-c"){
                conf= argv[i+1];
                clusters = init_cluster_num(conf);
                hash_func_num = init_hashFunc_num(conf);
                L = init_hashTable_num(conf);
            }
            else if (arg== "-d"){
                metric= argv[i+1];

            }
            else if (arg== "-o"){
                output= argv[i+1];
            }

        }
        return 1;
    }
    string get_inputFile(){
        return inputFile;
    }
    string get_confFile(){
        return conf;
    }
    string get_metric(){
        return metric;
    }
    string get_outFile(){
        return output;
    }
    int get_clusters(){
        return clusters;
    }
    int get_L(){
        return L;
    }
    int get_hash_func_num(){
        return hash_func_num;
    }
    //initializes memory for input data
    void init_data(string input);
    //get all the values of the data points
    void set_data_values(string input);
    //gets cluster number from configuration file
    int init_cluster_num(string conf);
    //gets hashfunc number from configuration file
    int init_hashFunc_num(string conf);
    //gets hashtable number from configuration file
    int init_hashTable_num(string conf);

    vector<double>* get_data(){
        return data_points;
    };

    vector<item*> get_data_items(){
        return myItems;
    };


    int get_num_items(){
        return num_of_data_points;
    };

    int get_point_length(){
        return point_length;
    };
private:
    string inputFile;
    string conf;
    string output;
    string metric;
    int clusters;
    int L;
    int hash_func_num;
    int num_of_data_points;
    int point_length;

    vector<double>* data_points;
    vector<item*> myItems; //vector storing all the items of the dataset (vectors) used for lsh and hypercube


};


#endif //ERGASIA2_ARGCLASS_H
