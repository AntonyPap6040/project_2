//
// Created by antonis on 24/10/2018.
//
#ifndef ERGASIA1_EUCLIDEAN_H
#define ERGASIA1_EUCLIDEAN_H

#include <unordered_map>
#include <iostream>
#include "item.h"
#include <random>


bool is_number(const std::string& s);

double eucl_dist(item p,item q);

/*FUNCTION OVERLOADING BELOW FOR UNORDERED SET*/
// Custom comparator
struct euclComparator
{
    bool operator()(const item & obj1, const item & obj2) const
    {
        if (obj1.get_point().size()==obj2.get_point().size())
            return true;
        return false;
    }
};

//custom hasher
struct euclHasher{
public:

    std::size_t operator ()(const item & obj) const{
            return  hash<string>()(obj.g());
    }

};

#endif
