#include <stdio.h>
#include <iostream>
#include <list>
#include "muscle2pop.h"
#include <tuple>
#include <stdexcept>
#include <gsl/gsl_interp.h>

void parameter::set(std::list<double> timep,std::list<double> rate){
    if (timep.size() != rate.size()){
        throw std::invalid_argument(  "time and rate not same length" );   
        return;
    }
    unsigned int size = rate.size();
    tp = new double[size];
    r = new double[size];
    unsigned int pos=0;
    for (auto elm : timep){
        tp[pos] = elm;
        pos++;
    }
    pos=0;
    for (auto elm : rate){
        r[pos] = elm;
        pos++;
    }
    interp = gsl_interp_alloc(gsl_interp_linear, size);
    gsl_interp_init(interp, tp,r,size);
}
double parameter::get(double now){
    double res;
    int e = gsl_interp_eval_e(interp, tp,r, now, NULL,&res);
    return res;
}

void T1::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    bool de = p.r_dm1();
    
    if (p.r_dm1()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
    }
}


void T2::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    bool de = p.r_dm2();
    
    if (p.r_dm2()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
    }
}


std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>>  test( parameters &p, double delta_t,
    double delta_t_out, double maxT, unsigned int seed, unsigned int num_t1, unsigned int num_t2){
    std::list<cell_base*> current_cells;
    std::list<cell_base*> old_cells;
    
    //std::random_device rd;
    std::mt19937 gen(seed);
    
    // inizsalizing
    p.set_random(&gen,delta_t);
    p.max_age = maxT;
 
    std::vector<std::vector<int>> t1_time_age;
    std::vector<std::vector<int>> t2_time_age; 
    //#define DT  0.1
    //#define OODT  1/DT
    double DT = delta_t_out;
    double OODT = 1/delta_t_out;
    unsigned int t_int = maxT*OODT+1.5;
    for (unsigned int i=0;i<t_int;i++){
        t1_time_age.push_back(std::vector<int> (t_int,0));
        t2_time_age.push_back(std::vector<int> (t_int,0));
    }
    //init cells
    for (unsigned int i=0;i<num_t1;i++){
        current_cells.push_back(new T1(p,0));
    }
     for (unsigned int i=0;i<num_t2;i++){
       current_cells.push_back(new T2(p,0));
    }


    for (double t=0;t<=maxT;t+=delta_t){
        p.get(t);
        auto it = current_cells.begin();
        unsigned int elem=0;
        while (it != current_cells.end()){
            auto it_c = it;
            auto element = *it;
            it++;
            elem++;
            element->events(it_c,current_cells,old_cells);        
        }
        if (p.r_fm1())
            current_cells.push_front(new T1(p,p.c_time));
        if (p.r_fm2())
            current_cells.push_front(new T2(p,p.c_time));
        
        if (elem>1e9) break;
        //clear old data and write to output array
        for (auto& oc : old_cells) {
            check_cell(DT, OODT, oc, t1_time_age, t2_time_age);
            delete oc;
        }
        old_cells.clear();
    }
     
    for (auto& cc : current_cells) {
        check_cell(DT, OODT, cc,t1_time_age,t2_time_age);
        delete cc;
    }
    std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> results(t1_time_age,t2_time_age);
    return results;
}

void check_cell(double DT, double OODT, cell_base *cell, std::vector<std::vector<int>>& t1_time_age,std::vector<std::vector<int>>& t2_time_age){
    unsigned int i_max = cell->t_death*OODT+0.5; // 1.5 for overlapping cell live
    if (i_max > t1_time_age.size())
        i_max = t1_time_age.size();
    if (cell->type == t_1){
        for (unsigned int i=cell->t_birth*OODT+0.5;i<i_max;i++){
            unsigned int age = (i*DT - cell->t_birth)*OODT+0.5;
            t1_time_age[i][age] += 1; 
        }
    }else{
        for (unsigned int i=cell->t_birth*OODT+0.5;i<i_max;i++){
            unsigned int age = (i*DT - cell->t_birth)*OODT+0.5;
            t2_time_age[i][age] += 1; 
        }
    }
}
