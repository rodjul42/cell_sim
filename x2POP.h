#include <stdio.h>
#include <list>
#include <gsl/gsl_interp.h>
#include <limits>
#include <random>

enum cell_type { t_q, t_p };
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
    b_lambdaq(NULL),
    b_lambdap(NULL),
    b_deltaq(NULL),
    b_deltap(NULL) {}
  ~parameters(){
    delete b_lambdaq;
    delete b_lambdap;
    delete b_deltaq;
    delete b_deltap;
  }
  void set_random(std::mt19937 *r_gen,double dt){
    delta_t = dt;
    gen = r_gen;
    delete b_lambdaq;
    b_lambdaq = new std::bernoulli_distribution(lambdaq*dt);
    delete b_lambdap;
    b_lambdap = new std::bernoulli_distribution(lambdap*dt); 
    delete b_deltaq;
    b_deltaq = new std::bernoulli_distribution(deltaq*dt);
    delete b_deltap;
    b_deltap = new std::bernoulli_distribution(deltap*dt);
  }
  void get(double now){
    c_time = now;
    }
  bool r_lambdaq(){return (*b_lambdaq)(*gen);}
  bool r_lambdap(){return (*b_lambdap)(*gen);}
  bool r_deltap() {return (*b_deltap) (*gen);}
  bool r_deltaq() {return (*b_deltaq) (*gen);}
  double max_age;
  double delta_t;
  double c_time;
  double lambdaq;
  double lambdap;
  double deltaq;
  double deltap;
  std::mt19937 *gen;
  std::bernoulli_distribution *b_lambdaq;
  std::bernoulli_distribution *b_lambdap;
  std::bernoulli_distribution *b_deltaq;
  std::bernoulli_distribution *b_deltap;
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


class CQ: public cell_base {
  public:
    CQ(parameters &A, double btime) : cell_base(A,btime,t_q) {}
    ~CQ() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

class CP: public cell_base {
  public:
    CP(parameters &A, double btime) : cell_base(A,btime,t_p) {}
    ~CP() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> test(parameters &p, double delta_t, double delta_t_out, double maxT, unsigned int seed, unsigned int num_q, unsigned int num_p);
void check_cell(double DT, double OODT, cell_base *cell, std::vector<std::vector<int>>& q_time_age,std::vector<std::vector<int>>& p_time_age);
