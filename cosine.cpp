//
// Created by antonis on 26/10/2018.
//

#include "cosine.h"

//norm function inspired from http://www.cplusplus.com/forum/beginner/227667/
double norm(vector<double> x)
{
    double  total = 0;
    for (int i = 0; i < x.size(); ++i)
    {
        total += x[i]*x[i];
    }

    return sqrt(total);
}


double cos_dist(item p,item q){
    vector<double> v1 = p.get_point();
    vector<double> v2 = q.get_point();
    double cospq=0.0;
    double S;

    /*//some processing if we have different dimensions
    if(v1.size()>v2.size()){
        for(int i=0;i<v1.size()-v2.size();i++){
            v2.push_back(0);
        }
    }
    else if (v1.size()<v2.size()){
        for(int i=0;i<v2.size()-v1.size();i++){
            v1.push_back(0);
        }
    }*/

    //calculating
    double m1 = norm(v1)*norm (v2);
    double m;

    for(int i=0; i<v1.size(); i++){
        m = v1[i]*v2[i];
        cospq = m / m1;
        S = 1 - cospq;
    }
    return S;

}