//
// Created by antonis on 29/10/2018.
//
#include "hyper.h"
#include "assignment.h"

#include <unordered_set>

int e=1;

using namespace std;
//c is not given so it is assumed to have the value of 1


item* assignment::cube( argClass args,  initialization init, assignment assign,item* old_clustering) {

    if(old_clustering!=NULL) done=true;

    ifstream ifs; //dataset stream
    ifstream query_file; //query stream
    ofstream output; //output stream

    int M=150; //max number of points to be checked
    int probes=1; //max number of vertices to be checked

    //get parameters

    int n = args.get_num_items();
    string metric = args.get_metric();
    int L = args.get_L();
    int k = args.get_hash_func_num();
    int c=1;
    double R=init.get_min()/2; //Radius given

    double d2 =3;



    vector<item*> myItems; //vector storing all the items of the dataset (vectors)
    for (int i = 0; i < args.get_clusters(); i++) {
        item* p = new item(i,k,L);
        p->set_coord(init.get_centroids()[i]);
        myItems.push_back(p);
    }

    vector<item*> myQuery; //vector storing all the items of the dataset (vectors)
    for (int q = 0; q < args.get_num_items(); q++) {
        item* p = new item(q,k,L);
        p->set_coord(args.get_data()[q]);
        myQuery.push_back(p);
    }

    item* m_range = new item[args.get_num_items()*L];
    int count_range=0;

    if(metric=="euclidean") {  //we are going to use the euclidean hypercube

        int y = 0;

        hypercube hyper_eucl[myItems.size()];

        //vector<double>* f=new vector<double>[myItems.size()];

        for (int j = 0; j < myItems.size(); j++) {
            //  hyper_eucl[j].set_size(myItems.size());
            //hyper_eucl[j].alloc(myItems.size());
            if (y == 0) {
                myItems[j]->set_V(*myItems[j]);
                myItems[j]->set_t();
                myItems[j]->the_rand_nums();
                if (y < myQuery.size()) {

                    //apothikeuoyme ta V kai t kai gia ta queries
                    myQuery[y]->set_V1(myItems[j]->getV());
                    myQuery[y]->set_t1(myItems[j]->gett());
                    myQuery[y]->the_rand_nums();

                }

            } else {
                myItems[j]->set_V1(myItems[j - 1]->getV());
                myItems[j]->set_t1(myItems[j - 1]->gett());
                myItems[j]->the_rand_nums();

            }

            myItems[j]->TableSize = n / 2;

            y++;
        }
        //gia ola ta queries
        myQuery[0]->TableSize = n / 2;
        for (y = 1; y < myQuery.size(); y++) {
            myQuery[y]->set_V1(myQuery[y - 1]->getV());
            myQuery[y]->set_t1(myQuery[y - 1]->gett());
            myQuery[y]->TableSize = n / 2;
            myQuery[y]->the_rand_nums();

        }

        // if (std::equal(v1.begin(), v1.begin() + n, v2.begin())
        vector<double> *FV = new vector<double>[myItems.size()];

        for (int j = 0; j < myItems.size(); j++) {

            calculate_f(FV[j], myItems[j], d2);
            int vertex = get_vertex(FV[j], d2);
            hyper_eucl[vertex].insert_vertex(FV[j], *myItems[j]);
        }

        delete[] FV;



        bool assign_finished = false;
        double num=R/args.get_clusters(); //range/ball

        if(num>args.get_num_items()) {
            cout<<"error 1";
            return m_range;
        }

        while(assign_finished==false) {

            int keep_last = count_range;

            int points_assigned_per_iteration = 0;

            for (int q = 0; q < myQuery.size(); q++) {

                vector<double> FV;
                calculate_f(FV, myQuery[q], d2);
                int vertex = get_vertex(FV, d2);
                int count_vertices = 0;
                int points_checked = 0;
                double dist = 99999;




                while (count_vertices < probes && points_checked < M) {

                    double dist1 = 0;
                    for (int i = 0; i < hyper_eucl[vertex].get_f().size(); i++) {

                        item v = hyper_eucl[vertex].get_item()[i];
                        dist1 = eucl_dist(v, *myQuery[q]);
                        //dist1 = ham_dist(v.get_point(), myQuery[q]->get_point());

                        if (double_equals(R,0.0)==false) {
                            if ((dist1 < e * R)) {

                                item z = *myQuery[q];
                                if (myQuery[q]->is_assigned() == false) {

                                    item z1 = v;

                                    points_assigned_per_iteration++;
                                    myQuery[q]->cluster_assign(z1.get_idi(), dist1);


                                    m_range[count_range] = *myQuery[q];

                                    count_range++;
                                }


                            }
                        }


                        points_checked++;
                        if (points_checked > M) {
                            break;
                        }
                    }

                    int min_dist_ham = 99999;

                    for (int i = 0; i < myItems.size(); i++) {

                        if (hyper_eucl[i].has_points() > 0) {
                            item v = hyper_eucl[i].get_item()[0];
                            dist1 = ham_dist(v.get_point(), myQuery[q]->get_point());

                            if (dist1 < min_dist_ham){
                                min_dist_ham = dist1;
                                vertex = i;
                            }
                        }
                    }
                    count_vertices++;


                }


            }

            for(int i=keep_last;i<count_range;i++ ) {
                double min1 = m_range[i].get_clust_dist();
                int pos=m_range[i].get_cluster_ball();

                if(m_range[i].is_changed()==false) {
                    for (int j = i + 1; j < count_range; j++) {

                        if (m_range[j].is_equal(m_range[i])) {
                            if (m_range[j].get_clust_dist() < min1) {
                                min1 = m_range[j].get_clust_dist();
                                pos = m_range[j].get_cluster_ball();
                            }
                        }
                    }
                    for (int j = i + 1; j < count_range; j++) {

                        if (m_range[j].is_equal(m_range[i]) && m_range[j].is_changed()==false) {
                            m_range[j].change();
                            m_range[i].change();

                            m_range[j].cluster_assign(pos, min1);
                            m_range[i].cluster_assign(pos, min1);
                        }
                    }


                }


            }

            if(count_range>=args.get_num_items() ) assign_finished=true;
            else if((points_assigned_per_iteration/num<0.85) && (count_range>0.7*args.get_num_items())) assign_finished=true;
            else{
                if((2 > init.get_max())) R = R* 2/init.get_max();
                else R = R*2;

                num=R/args.get_clusters();
                if(num>args.get_num_items() || (num>init.get_max()))break;
            }



            /*

            for(int i=0; i<myQuery.size();i++){
                myQuery[i]->freeE();
                myQuery[i]->freeRand();
            }
            for(int i=0; i<myItems.size();i++){
                myItems[i]->freeE();
            }*/

        }


        for(int i=0; i<myQuery.size();i++){
            double min=99999;
            if(myQuery[i]->is_assigned() == false){
                double temp_dist;
                int num_center;
                for(int j =0 ; j<args.get_clusters(); j++){
                    temp_dist=eucl_dist(*myQuery[i],*myItems[j]);
                    if((temp_dist<min)){
                        min=temp_dist;
                        num_center = j;
                    }
                }
                myQuery[i]->cluster_assign(num_center,min);

                m_range[count_range] = *myQuery[i];

                count_range++;
            }
        }

    }


    if(metric=="cosine") {  //we are going to use the cosine hypercube
        int y = 0;

        hypercube hyper_cos[myItems.size()];



        y = 0;
        for (int j = 0; j < myItems.size(); j++) {

            if (y == 0) {
                myItems[j]->set_r(*myItems[j]);
                if (y < myQuery.size()) {

                    myQuery[y]->set_r1(myItems[j]->getr());
                }

            } else {
                myItems[j]->set_r1(myItems[j - 1]->getr());

            }

            myItems[j]->TableSize = n / 2;

            y++;
        }



        //gia ola ta queries
        myQuery[0]->TableSize = n / 2;
        for (y = 1; y < myQuery.size(); y++) {
            myQuery[y]->set_r1(myQuery[y - 1]->getr());
            myQuery[y]->TableSize = n / 2;

        }

        vector<double> *FV = new vector<double>[myItems.size()];

        for (int j = 0; j < myItems.size(); j++) {

            calculate_c(FV[j], myItems[j], d2);
            int vertex = get_vertex(FV[j], d2);
            hyper_cos[vertex].insert_vertex(FV[j], *myItems[j]);
        }

        delete[] FV;

        /*END OF PREPROCESSING*/

        /*RANGE SEARCH*/


        //for range search


        bool assign_finished = false;
        double num=R/args.get_clusters(); //range/ball

        if(num>args.get_num_items()) {
            cout<<"error 1";
            return m_range;
        }

        while(assign_finished==false) {

            int keep_last = count_range;

            int points_assigned_per_iteration = 0;

            for (int q = 0; q < myQuery.size(); q++) {

                vector<double> FV;
                calculate_c(FV, myQuery[q], d2);
                int vertex = get_vertex(FV, d2);
                int count_vertices = 0;
                int points_checked = 0;
                double dist = 99999;
                int keep_vertices[probes];
                int count_my_kept=0;
                double dist1 = 0;

                int min_dist_ham = 99999;

                if(probes>1){
                    for (int i = 0; i < myItems.size(); i++) {
                        if (hyper_cos[i].has_points() > 0) {

                            item v = hyper_cos[i].get_item()[0];
                            dist1 = ham_dist(v.get_point(), myQuery[q]->get_point());

                            if ((dist1 < min_dist_ham)) {
                                min_dist_ham = dist1;
                                keep_vertices[count_my_kept]=i;
                                count_my_kept++;
                                if(count_my_kept>=probes) count_my_kept=0;   //care
                            }
                        }
                    }
                }
                




                while (count_vertices < probes && points_checked < M) {

       
                    for (int i = 0; i < hyper_cos[vertex].get_f().size(); i++) {

                        item v = hyper_cos[vertex].get_item()[i];
                 
                        dist1 = cos_dist(v, *myQuery[q]);
                 

                        if (double_equals(R,0.0)==false) {
                            if ((dist1 <e * R)) {

                                item z = *myQuery[q];
                                if (myQuery[q]->is_assigned() == false) {

                                    item z1 = v;

                                    points_assigned_per_iteration++;
                                    myQuery[q]->cluster_assign(z1.get_idi(), dist1);


                                    m_range[count_range] = *myQuery[q];

                                    count_range++;
                                }


                            }
                        }


                        points_checked++;
                        if (points_checked > M || count_vertices>=probes) {
                            break;
                        }
                    }



                    vertex = keep_vertices[count_vertices];
                    count_vertices++;


                }


            }
            for(int i=keep_last;i<count_range;i++ ) {
                double min1 = m_range[i].get_clust_dist();
                int pos=m_range[i].get_cluster_ball();
                //cout<<pos<<"endd"<<endl;
                if(m_range[i].is_changed()==false) {
                    for (int j = i + 1; j < count_range; j++) {

                        if (m_range[j].is_equal(m_range[i])) {
                            if ((m_range[j].get_clust_dist() < min1)){
                                min1 = m_range[j].get_clust_dist();
                                pos = m_range[j].get_cluster_ball();
                            }
                        }
                   
                        if (m_range[j].is_equal(m_range[i]) ) {
                            m_range[j].change();
                            m_range[i].change();

                            m_range[j].cluster_assign(pos, min1);
                            m_range[i].cluster_assign(pos, min1);
                        }
                    }


                }

                // cout<<pos<<" k "<<endl;

            }

            if(count_range>=args.get_num_items() ) assign_finished=true;
            else if((points_assigned_per_iteration/num<0.85) && (count_range>0.7*args.get_num_items())) assign_finished=true;
            else{

                if((2 > init.get_max())) R = R* 2/init.get_max();
                else R = R*2;

                num=R/args.get_clusters();
                if(num>args.get_num_items() || (num>init.get_max()))break;
            }

            /*
            for(int i=0; i<myQuery.size();i++){
                myQuery[i]->freeC();
            }
            for(int i=0; i<myItems.size();i++){
                myItems[i]->freeC();
            }
    */
        }

        for(int i=0; i<myQuery.size();i++){
            double min=99999;
            if(myQuery[i]->is_assigned() == false){
                double temp_dist;
                int num_center;
                for(int j =0 ; j<args.get_clusters(); j++){
                    temp_dist=eucl_dist(*myQuery[i],*myItems[j]);
                    if((temp_dist<min)){
                        min=temp_dist;
                        num_center = j;
                    }
                }
                myQuery[i]->cluster_assign(num_center,min);

                m_range[count_range] = *myQuery[i];

                count_range++;
            }
        }

    }

/*
    for(int i=0; i<myItems.size();i++){
        delete myItems[i];
    }


    for(int i=0; i<myQuery.size();i++){
        delete myQuery[i];
    }
*/

    for (int i=0; i<count_range;i++){
        if(m_range[i].get_cluster_ball()<args.get_clusters()) {
            clusters[m_range[i].get_cluster_ball()].cluster_data[clusters[m_range[i].get_cluster_ball()].counter] = m_range[i].get_point();
            clusters[m_range[i].get_cluster_ball()].counter++;
        }


    }
    if(old_clustering!=NULL) {
        int count_d=0;
        for (int i = 0; i < count_range; i++) {

            for (int j=0; j<count_range;j++) {


                if(m_range[j].is_equal(old_clustering[i]) && m_range[j].get_cluster_ball()!=old_clustering[i].get_cluster_ball()){
                    count_d++;

                    break;

                }

            }

            if(count_d>190) {
                done = false;
                break;
            }
        }
       // cout<<count_d<<endl;

    }

    return m_range;
}

