//
// Created by antonis on 28/11/2018.
//

#include "argClass.h"
#include <algorithm>
#include <iterator>
#include <sstream>      // std::stringstream
#include "euclidean.h"


//configuration file

int argClass::init_cluster_num(string conf){
    ifstream ifs; //dataset stream
    string line;

    ifs.open(conf,std::ifstream::in);
    if(ifs.is_open()) {

        int count_lines=0;
        while ( getline (ifs,line) ){
            if(count_lines==0){
                string delimiter = ":";

                size_t pos = 0;
                std::string token;
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(pos+1, line.length());
                    if (token.length()>0)return stoi(token);
                    else return 15;
                }
            }
            count_lines++;
        }

    }
    else return 15;
    ifs.close();
}


int argClass::init_hashFunc_num(string conf){
    ifstream ifs; //dataset stream
    string line;

    ifs.open(conf,std::ifstream::in);
    if(ifs.is_open()) {

        int count_lines=0;
        while ( getline (ifs,line) ){
            if(count_lines==1){
                string delimiter = ":";

                size_t pos = 0;
                std::string token;
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(pos+1,line.length());
                    if(token.length()>0)return stoi(token);
                    else return 4;
                }
            }
            count_lines++;

        }

    }
    else return 4;
    ifs.close();

}

int argClass::init_hashTable_num(string conf){
    ifstream ifs; //dataset stream
    string line;
    ifs.open(conf,std::ifstream::in);
    if(ifs.is_open()) {

        int count_lines=0;
        while ( getline (ifs,line) ){
            if(count_lines==2){
                string delimiter = ":";

                size_t pos = 0;
                std::string token;
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(pos+1, line.length());
                    if (token.length()>0) return stoi(token);
                    else return 5;
                }
            }
            count_lines++;

        }

    }
    else return 5;
    ifs.close();
}

//--

void argClass::init_data(string input){
    std::ifstream myfile(input);

    // new lines will be skipped unless we stop it from happening:
    myfile.unsetf(std::ios_base::skipws);

    // count the newlines with an algorithm specialized for counting:
    unsigned line_count = std::count(
            std::istream_iterator<char>(myfile),
            std::istream_iterator<char>(),
            '\n');

    std::cout << "Lines: " << line_count << "\n";

    data_points = new vector<double>[line_count];


}

void argClass::set_data_values(string input){
    ifstream ifs; //dataset stream
    string line;
    ifs.open(input,std::ifstream::in);
    if(ifs.is_open()) {

        int count_lines=0;
        while ( getline (ifs,line) ){
            stringstream ss(line);

            std::string token;
            int count_stream=0;
            item *itemId = new item;

            while(std::getline(ss, token, ',')) {

                if (count_stream==0) {

                    bool is_num = is_number(token);
                    if (is_num) {
                        itemId->set(stoi(token));
                    } else {
                        itemId->set(token);

                    }
                }
                else {

                    double i;
                    i = stod(token);

                    data_points[count_lines].push_back(i);

                    itemId->set_coord(data_points[count_lines]);

                    myItems.push_back(itemId);

                }
                count_stream++;

            }
            point_length=count_stream;
            count_lines++;
        }
        num_of_data_points = count_lines;


    }

}
