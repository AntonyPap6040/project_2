//
// Created by antonis on 28/11/2018.
//

#include <fstream>
#include "item.h"
#include "euclidean.h"
#include "cosine.h"
#include <sstream>      // std::stringstream
#include <algorithm>

#include "argClass.h"
#include "initialization.h"
#include "assignment.h"


//c is not given so it is assumed to have the value of 1



int main(int argc, char* argv[]) {
    std::cout.precision(17);

    argClass args;

    if( args.give_args(argc,argv)==-1) return -1;

    cout<<args.get_L()<<endl;
    vector<double>* v = args.get_data();
    vector<item*> v2 = args.get_data_items();




    //1.
    initialization init(args.get_metric(), args.get_num_items(),args.get_clusters());
    init.random_k_points(args.get_clusters(),args.get_point_length(),args.get_data());
    assignment assign(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());


     cout<<assign.get_clusters()[0].counter<<endl;
     cout<<init.get_centroids()[0][0]<<endl;

    //2-1-1...
    assign.m1(args,init,assign);
    assign.m2(args,init,assign);
    assign.m3(args,init,assign);
    assign.m4(args,init,assign);
    assign.m5(args,init,assign);
    assign.m6(args,init,assign);


    //2.
    initialization init1(args.get_metric(), args.get_num_items(),args.get_clusters());
    init1.k_means_pp(args.get_clusters(),args.get_point_length(),args.get_data());
    assignment assign1(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());

     cout<<assign1.get_clusters()[0].counter<<endl;
     cout<<init1.get_centroids()[0][0]<<endl;

    assign1.m7(args,init1,assign1);
    assign1.m8(args,init1,assign1);
    assign1.m9(args,init1,assign1);
    assign1.m10(args,init1,assign1);
    assign1.m11(args,init1,assign1);
    assign1.m12(args,init1,assign1);


    cout<<"the end"<<endl;



}