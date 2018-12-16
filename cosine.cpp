//
// Created by antonis on 26/10/2018.
//

#include "cosine.h"

//the source below
//https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
float sqrt3(const double x)  
{
  union
  {
    int i;
    double x;
  } u;

  u.x = x;
  u.i = (1<<29) + (u.i >> 1) - (1<<22); 
  return u.x;
} 



double cos_dist(item p,item q){
    vector<double> v1 = p.get_point();
    vector<double> v2 = q.get_point();
    double S=0.0;

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

    double m=0.0;


    double  total1 = 0.0;
    double total2=0.0;
    for (int i = 0; i < v1.size(); ++i)
    {
        total1 += v1[i]*v1[i];
        total2 += v2[i]*v2[i];
        m += v1[i]*v2[i];

    }


    double normv1= sqrt3(total1);
    double normv2= sqrt3(total2);
    double m1 = normv1*normv2;

    double cospq = m / m1;
    S = 1-cospq;
       

    return S;

}