//
// Created by antonis on 28/11/2018.
//

#include "initialization.h"

bool double_equals(double a, double b, double epsilon )
{
    return std::abs(a - b) < epsilon;
}

bool double_less(double a, double b, double epsilon )
{
    return std::abs(a - b) > epsilon;
}

bool double_greater(double a, double b, double epsilon )
{
    return std::abs(a - b) < epsilon;
}




void initialization::calculate_min_max(int k){

    double min=99999.0;
    double max=-99999.0;

    for(int i =1 ; i<k; i++){
        item p,q;

        q.set_coord(cluster_centroid_points[i-1]);
        p.set_coord(cluster_centroid_points[i]);
        double temp_dist;

        if(metric == "euclidean") temp_dist = eucl_dist(p,q);
        else if(metric == "cosine") temp_dist= cos_dist(p,q);

        if(temp_dist<min && temp_dist!=0){
            min = temp_dist;
        }

        if(temp_dist>max && temp_dist!=0){
            max = temp_dist;
        }
    }

    min_dist_centers = min;
    max_dist_centers = max;

}



void initialization::random_k_points(int k, int point_length,vector<double>* data_points){


    int mino=0;
    int maxo=data_num;
    int t=0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distribution(mino ,maxo); //uniformly distributed

    int r = floor(distribution(generator));
    int rr[k];

    //step 1
    for(t=0; t<k ; t++){
        rr[t]=r;
        int flag = 0;
        while(flag==0) {
            r =  floor(distribution(generator));
            int found=0;
            for (int i = 0; i <= t; i++) {
                if (r == rr[i]) {
                    found = 1;
                }
            }
            if(found==0) {
                flag=1;
            }
        }
        if(double_equals(data_points[r][0],0.0) && double_equals(data_points[r][1],0.0)){
            cluster_centroid_points[t] = data_points[r-1];
        }
        else{
            cluster_centroid_points[t] = data_points[r];

        }



    }
    calculate_min_max(k);
};




void initialization::k_means_pp(int k, int point_length,vector<double>* data_points){
    int mino=0;
    int maxo=data_num;
    int t=0;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    std::default_random_engine generator(seed);
    std::uniform_real_distribution<float> distribution(mino ,maxo); //uniformly distributed

    int r = floor(distribution(generator));
    //k-means++ algortihm for initialization

    //step 1
    if(t==0){   //step 1
        cluster_centroid_points[t] = data_points[r];
    }

    int n = data_num;
    double P[data_num+1];
    int rr[k];

    rr[t]=r;
   // cout<<rr[0]<<endl;
    t++;

    //steps 2-4
    while(t<k){
        double min[data_num] ;
        double max[data_num];

        for(int i=0; i<n-t; i++){
            item p,q;
            int cont=0;
            for (int l = 0; l < t; l++) {
                if (i == rr[l]) {
                    cont=1;
                    break;
                }
            }
            if(cont) continue;


            p.set_coord(data_points[i]);
            max[i]=-9999;
            min[i]=9999;

            double temp_distance=0.0;
            int min_centroid=0;
            double dist_sum=0.0; //D(i)^2

            for(int j=0; j<t; j++) {
                q.set_coord(cluster_centroid_points[j]);

                if (metric == "euclidean") {
                    temp_distance = eucl_dist(p,q);
                }
                else if(metric=="cosine"){
                    temp_distance = cos_dist(p,q);

                }

                if(temp_distance==0){
                    continue;
                }

                if(i>0) {
                    if (temp_distance < min[i] ) {
                        min[i] = temp_distance;
                        min_centroid = j;

                    }
                }
                else{
                    if (temp_distance < min[i] && min[i] != min[i - 1] ) {
                        min[i] = temp_distance;
                        min_centroid = j;

                    }
                }
                if(temp_distance>max[i]){
                    max[i] = temp_distance;
                }
            }

            dist_sum+=min[i]*min[i];
            P[i] = dist_sum;

        }


        for(int i=0; i<n-t; i++){
            P[i] = P[i] / max[i];
        }

        std::uniform_real_distribution<double> distributionX(0.0, P[n - t]); //uniformly distributed


        double x = distributionX(generator);

        for (int i = 0; i <= n - t; i++) {

            int cont=0;
            for (int l = 0; l < t; l++) {
                if (i == rr[l]) {
                    cont=1;
                    break;
                }
            }
            if(cont) continue;


            if(i==0){
                if (P[n - t] < x && x <= P[i]) {
                    r = i;
                }
            }
            else {
                if (P[i - 1] < x && x <= P[i]) {
                    r = i;
                }
            }

        }


        rr[t]=r;
        if(double_equals(data_points[r][0],0.0) && double_equals(data_points[r][1],0.0)){
            cluster_centroid_points[t] = data_points[r-1];
        }
        else{
            cluster_centroid_points[t] = data_points[r];

        }
        t++;



    }

    calculate_min_max(k);

}
