#include <stdio.h>
#include <list>
#include <gsl/gsl_interp.h>
#include <limits>
#include <random>

enum cell_type { t_1, t_2 };
enum event { die, live };
class parameter{
   public:
  ~parameter(){
      gsl_interp_free(interp);
      delete[] tp;
      delete[] r;
  }
  void set(std::list<double> timep, std::list<double> rate);
  double get(double now);
  private:
    gsl_interp *interp;
    double *tp,*r;
};


class parameters {
 public:
  parameters():
    b_fm1(NULL),
    b_fm2(NULL),
    b_dm1(NULL),
    b_dm2(NULL) {}
  ~parameters(){
    delete b_dm1;
    delete b_dm2;
    delete b_fm1;
    delete b_fm2;
  }
  void set_random(std::mt19937 *r_gen,double dt){
    delta_t = dt;
    gen = r_gen;
    delete b_dm1;
    b_dm1 = new std::bernoulli_distribution(dm1*dt);
    delete b_dm2;
    b_dm2 = new std::bernoulli_distribution(dm2*dt);
  }
  void get(double now){
    c_time = now;
    fm1 = fm1_t.get(now);
    fm2 = fm2_t.get(now);
    delete b_fm1;
    b_fm1 = new std::bernoulli_distribution(fm1*delta_t);
    delete b_fm2;
    b_fm2 = new std::bernoulli_distribution(fm2*delta_t);
    }
  bool r_fm1(){return (*b_fm1)(*gen);}
  bool r_fm2(){return (*b_fm2)(*gen);}
  bool r_dm1() {return (*b_dm1) (*gen);}
  bool r_dm2() {return (*b_dm2) (*gen);}
  parameter fm1_t;
  parameter fm2_t;
  double max_age;
  double delta_t;
  double c_time;
  double fm1;
  double fm2;
  double dm1;
  double dm2;
  std::mt19937 *gen;
  std::bernoulli_distribution *b_fm1;
  std::bernoulli_distribution *b_fm2;
  std::bernoulli_distribution *b_dm1;
  std::bernoulli_distribution *b_dm2;
  private:   
};

class cell_base {
  public:
    cell_base(parameters &A, double btime, cell_type c_t) : 
        p(A),
        type(c_t),
        t_birth(btime),
        t_death(A.max_age + 100) {} 
    cell_base(parameters &&A, double btime, cell_type c_t) = delete; 
    virtual ~cell_base(){}

    virtual void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){std::cout<<"virt;";return;}
    parameters &p;
    cell_type type;
    double t_birth;
    double t_death;
  private:
};


class T1: public cell_base {
  public:
    T1(parameters &A, double btime) : cell_base(A,btime,t_1) {}
    ~T1() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

class T2: public cell_base {
  public:
    T2(parameters &A, double btime) : cell_base(A,btime,t_2) {}
    ~T2() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> test(parameters &p, double delta_t, double delta_t_out, double maxT, unsigned int seed, unsigned int num_2n, unsigned int num_4n);
void check_cell(double DT, double OODT, cell_base *cell, std::vector<std::vector<int>>& t1_time_age,std::vector<std::vector<int>>& t2_time_age);
