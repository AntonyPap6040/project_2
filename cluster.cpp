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




    //2-1-1...
  /*  assign.m1(args,init,assign);
    cout<<"m1 ok"<<endl;
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    fflush(stdout);*/

    assign.m2(args,init,assign);
    cout<<"m2 ok"<<endl;
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    fflush(stdout);
 

   /*    assign.m3(args,init,assign);
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    cout<<"m3 ok"<<endl;
    fflush(stdout);

   /* assign.m4(args,init,assign);
    cout<<"m4 ok"<<endl;
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    fflush(stdout);

    assign.m5(args,init,assign);
    cout<<"m5 ok"<<endl;
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    fflush(stdout);

    assign.m6(args,init,assign);
    cout<<"m6 ok"<<endl;
    assign.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init.get_centroids());
    fflush(stdout);*/

    //2.
   /* initialization init1(args.get_metric(), args.get_num_items(),args.get_clusters());
    init1.k_means_pp(args.get_clusters(),args.get_point_length(),args.get_data());
    assignment assign1(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());



    assign1.m7(args,init1,assign1);
    cout<<"m7 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout);   
    assign1.m8(args,init1,assign1);
    cout<<"m8 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout); */

/*
    assign1.m9(args,init1,assign1);
    cout<<"m9 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout);

    assign1.m10(args,init1,assign1);
    cout<<"m10 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout);

    assign1.m11(args,init1,assign1);
    cout<<"m11 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout);

   assign1.m12(args,init1,assign1);
    cout<<"m12 ok"<<endl;
    assign1.clean(args.get_clusters(),args.get_num_items(),args.get_metric(),init1.get_centroids());
    fflush(stdout);*/


    cout<<"the end"<<endl;



}