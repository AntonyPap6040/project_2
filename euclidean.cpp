//
// Created by antonis on 24/10/2018.
//


#include "euclidean.h"

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


double eucl_dist(item p,item q){
    vector<double> v1 = p.get_point();
    vector<double> v2 = q.get_point();
    double S=0.0;

    if(v1.size()>v2.size()){
        for(int i=0;i<v1.size()-v2.size();i++){
            v2.push_back(0);
        }
    }
    else if (v1.size()<v2.size()){
        for(int i=0;i<v2.size()-v1.size();i++){
            v1.push_back(0);
        }
    }

    for(int i=0; i<v1.size(); i++){
        double d = v1[i]-v2[i];
        S+=d*d;
    }

    return sqrt(S);

}