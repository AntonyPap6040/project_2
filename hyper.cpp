//
// Created by antonis on 29/10/2018.
//

#include "hyper.h"
#include <bitset>

void calculate_f(vector<double> &FV,item *item1,int d2) {
    //f()
    //calculating the f:[f1(h1(p)...fd(hd(p)]

    long long hi = item1->hie(item1, 0);

    for (int i = 0; i < d2; i++) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        //unsigned long long seed = item1->hie(item1, 0);
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distributionU(0.0, 1.0); //uniformly distributed
        double fi;


        fi = distributionU(generator);
        long long check = floor(10*fi + hi);

        //paradoxi oti an exo to check perito tote tha valo timi 1 gia to fi(hi)
        if(check % 2 ==1){
            fi = 1;

        }
        else{
            fi = 0;

        }

        FV.push_back(fi);

    }
}

void calculate_c(vector<double> &FV,item *item1,int d2) {

    vector<double> d;
    for (int i = 0; i < d2; i++) {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        //unsigned long long seed = item1->hie(item1, 0);
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distributionU(0.0, 1.0); //uniformly distributed
        double fi;


        fi = distributionU(generator);
        FV.push_back(fi*item1->hic(item1,0));

    }
}


int get_vertex(vector<double> f,int d2){

    string binary;
    int vertex=0;
    for(int i=d2-1; i>=0; --i){
        binary.append(to_string(f[i]));
        if(f[i]==1) {
            double pow=1.0;

            for (int j = 0; j < i - d2 + 1; j++) {
                pow *= pow * f[i] * 2;
            }
            vertex+=pow;

        }


    }
    //cout<<vertex<<"aa";
    return floor(vertex);


}



double ham_dist(vector<double> v1, vector<double> v2) {

    int dist=0;
    if(v1.size()<v2.size()) {
        for (int j = 1; j < v1.size(); j++) {
            dist+= bitset<64>((long long) v1[j] ^ (long long)v2[j]).count();
        }
        dist += (v2.size() - v1.size())*4;
        //cout<<dist<<endl;
    }
    else{
        for (int j = 1; j < v2.size(); j++) {
            dist+= bitset<64>((long long) v1[j] ^ (long long)v2[j]).count();
        }
        dist += (v1.size() - v2.size())*4;
    }
    // cout<<dist<<endl;
    return dist;
}




//pigi parakato sinartisis https://stackoverflow.com/questions/31223985/c-double-to-binary-representation-reinterpret-cast
/*
string toBinary( char a)
{
    string s1;
    const size_t size = sizeof(a) * 8;
    for (int i = size - 1; i >= 0; --i){
        bool b = a & (1UL << i);
        s1.append(to_string(b));
    }
    return s1;
}

string toBinary( double d)
{
    string s;
    const size_t size = sizeof(d);
    for (int i = 0; i < size; ++i){
        char* c = reinterpret_cast<char*>(&d) + i;
        s.append(toBinary(*c));
    }
    return  s;
}


double ham_dist(vector<double> v1, vector<double> v2){

    int dist=0;
    if(v1.size()<v2.size()) {
        for (int j = 1; j < v1.size(); j++) {
            string s1,s2;

            s1 = toBinary(v1[j]);
            s2 = toBinary(v2[j]);

            if(s1.size()<s2.size()) {
                for (int i = 0; i < s1.size(); i++) {
                    if (s1[i] != s2[i]) {
                        dist++;
                    }
                }
                dist+=s2.size()-s1.size();

            }
            else{
                for (int i = 0; i <s2.size(); i++) {
                    if (s1[i] != s2[i]) {
                        dist++;
                    }
                }
                dist+=s1.size()-s2.size();
            }
        }
    }
    else{

        for (int j = 1; j < v2.size(); j++) {
            string s1,s2;

            s1 = toBinary(v1[j]);
            s2 = toBinary(v2[j]);

            if(s1.size()<s2.size()) {
                for (int i = 0; i < s1.size(); i++) {
                    if (s1[i] != s2[i]) {
                        dist++;
                    }
                }
                dist+=s2.size()-s1.size();

            }
            else{
                for (int i = 0; i <s2.size(); i++) {
                    if (s1[i] != s2[i]) {
                        dist++;
                    }
                }
                dist+=s1.size()-s2.size();
            }
        }
    }

    return dist;


}*/
