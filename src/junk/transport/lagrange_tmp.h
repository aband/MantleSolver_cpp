#ifndef LAGRANGE_TEMP_H_
#define LAGRANGE_TEMP_H_

#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

static inline int combinatorial_monomial_map(int i, int j) { return (i+j)*(i+j+1)/2 + j; }

class Combinatorials {
  private:
    static int my_num;
    static std::vector<double> my_combinatorials;

  public:
    // Constructors
    void init(int new_num);
    void extend(int new_num) { init(new_num); }
    
    Combinatorials() {};
    Combinatorials(int new_num) { init(new_num); };

    void clear();

    // Access functions
    int num() const { return my_num; }
    int monomial_map(int n, int k) const { return combinatorial_monomial_map(n,k); }

    double eval(int n, int k);
    double comb(int n, int k) { return eval(n,k); }
    double combinations(int n, int k) { return eval(n,k); }
    double choose(int n, int k) { return eval(n,k); }
    double operator() (int n, int k) { return eval(n,k); }

    // Output for testing
    void write_raw(std::ofstream& fout) const;
    int write_raw(std::string& filename) const;
  };

// ================================================================================

static inline int lagrangeBasisDeriv_map(int n, int i) { return n*(n+1)/2 + i; }
static inline int lagrangeBasisDeriv_map_odd(int n, int i) { return (n+1)*(n-1)/4 + i; }

class LagrangeBasisDeriv {
private:
  static int max_degree;
  static std::vector<double> my_lagrangeBasisDeriv_middle;
  static std::vector<double> my_lagrangeBasisDeriv_endHi;

public:
  // Constructors
  void init(int new_max_degree);
  void extend(int new_max_degree) { init(new_max_degree); }

  LagrangeBasisDeriv() {};
  LagrangeBasisDeriv(int new_max_degree) { init(new_max_degree); };

  void clear();

  // Access functions
  int maxDegree() const { return max_degree; }
  int map(int n, int k) const { return lagrangeBasisDeriv_map(n,k); }
  int map_odd(int n, int k) const { return lagrangeBasisDeriv_map_odd(n,k); }

  double middle(int n, int k);
  double endLo(int n, int k);
  double endHi(int n, int k);

  void write_raw(std::ofstream& fout) const;
  int write_raw(std::string& filename) const;
};

#endif
