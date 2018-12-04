//
// Created by antonis on 26/10/2018.
//

#ifndef ERGASIA1_COSINE_H
#define ERGASIA1_COSINE_H

#include <unordered_map>
#include <iostream>
#include "item.h"
#include <random>

/*FUNCTION OVERLOADING BELOW FOR UNORDERED SET*/
// Custom comparator

double cos_dist(item p,item q);


struct cosComparator
{
    bool operator()(const item & obj1, const item & obj2) const
    {
        if (obj1.get_point().size()==obj2.get_point().size())
            return true;
        return false;
    }
};

//custom hasher (cosine)
struct cosHasher{
public:

    std::size_t operator ()(const item & obj) const{
        return  hash<string>()(obj.G());
    }

};

#endif //ERGASIA1_COSINE_H
