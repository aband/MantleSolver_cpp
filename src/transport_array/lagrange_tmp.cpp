// Lagrange helper function from Todd

#include "lagrange_tmp.h"

int Combinatorials::my_num = -1;
std::vector<double> Combinatorials::my_combinatorials;

void Combinatorials::init(int new_num) {
  if(my_combinatorials.size() == 0) my_num = -1;
  
  if(new_num <= my_num) return;
  my_num = (new_num > 8) ? new_num + 2 : 10;

  // Compute combinatorials (Binomial coefficients) from Pascal's triangle

  my_combinatorials.clear();
  my_combinatorials.resize(monomial_map(my_num,my_num)+1);

  for(int n = 0; n <= my_num; n++) {
    my_combinatorials[monomial_map(n,0)] = 1;
    my_combinatorials[monomial_map(n,n)] = 1;
    
    for(int k = 1; k < n; k++) {
      my_combinatorials[monomial_map(n,k)]
	= my_combinatorials[monomial_map(n-1,k-1)] + my_combinatorials[monomial_map(n-1,k)];
    }
  }
}

void Combinatorials::clear() {
  my_combinatorials.clear();
  my_num = -1;
}

double Combinatorials::eval(int n, int k) {
  if(k < 0 || k > n) return 0;
  if(n > my_num) extend(n);

  return my_combinatorials[monomial_map(n,k)];
}

void Combinatorials::write_raw(std::ofstream& fout) const {
  fout << "my_num = " << my_num << endl;
  fout << "my_size = " << my_combinatorials.size();

  if(my_combinatorials.size() > 0) {
    fout << "Combinatorials:";
    for(int n = 0; n <= my_num; ++n) {
      fout << "\n n = " << n << ":";
      for(int k = 0; k <= n; ++k) {
	fout << " " << my_combinatorials[monomial_map(n,k)];
      }
    }
  }
}

int Combinatorials::write_raw(string& filename) const {
  std::ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}

// =====================================================================

int LagrangeBasisDeriv::max_degree = -1;
std::vector<double> LagrangeBasisDeriv::my_lagrangeBasisDeriv_middle;
std::vector<double> LagrangeBasisDeriv::my_lagrangeBasisDeriv_endHi;

void LagrangeBasisDeriv::init(int new_max_degree) {
  if(new_max_degree <= max_degree) return;
  max_degree = (new_max_degree > 5) ? new_max_degree + 1 + (new_max_degree % 2) : 7; // odd!

  Combinatorials comb;

  // Compute my_lagrangeBasisDeriv_endHi

  my_lagrangeBasisDeriv_endHi.clear();
  my_lagrangeBasisDeriv_endHi.resize(map(max_degree+1,max_degree+1));

  my_lagrangeBasisDeriv_endHi[0] = 0;
  
  double sumJInv = 0;
  for(int n = 1; n <= max_degree; n++) {
    sumJInv += 1.0 / (double)n;

    my_lagrangeBasisDeriv_endHi[map(n,n)] = sumJInv;

    double sign = (n % 2) ? -1 : 1;
    for(int i=0; i<n; i++) {
      my_lagrangeBasisDeriv_endHi[map(n,i)] = sign * comb(n,i) / (double)(n-i);
      sign = -sign;
    }
  }
  
  // Compute my_lagrangeBasisDeriv_middle (n odd!)

  my_lagrangeBasisDeriv_middle.clear();
  my_lagrangeBasisDeriv_middle.resize(map_odd(max_degree+2,max_degree+2));

  for(int n = 1; n <= max_degree; n+=2) {
    
    double constantFactor = pow(2,1-n);
    for(double k=3; k<=n; k+=2) {
      constantFactor *= k/(k-1);
    }
    if((n+1)/2 % 2) constantFactor = -constantFactor;

    for(int i = 0; i <= n; i++) {
      my_lagrangeBasisDeriv_middle[map_odd(n,i)] = constantFactor * comb(n,i) / pow((n - 2*i),2);
      constantFactor = -constantFactor;
    }
  }
}

void LagrangeBasisDeriv::clear() {
  my_lagrangeBasisDeriv_middle.clear();
  my_lagrangeBasisDeriv_endHi.clear();
  max_degree = -1;
}

double LagrangeBasisDeriv::endHi(int n, int i) {
  if(i < 0 || i > n) return 0;
  if(n > max_degree) extend(n);
  return my_lagrangeBasisDeriv_endHi[map(n,i)];
}

double LagrangeBasisDeriv::endLo(int n, int i) {
  if(i < 0 || i > n) return 0;
  if(n > max_degree) extend(n);
  return -my_lagrangeBasisDeriv_endHi[map(n,n-i)];
}

double LagrangeBasisDeriv::middle(int n, int i) {
  if(i < 0 || i > n) return 0;
  if(n > max_degree) extend(n);
  return my_lagrangeBasisDeriv_middle[map_odd(n,i)];
}

void LagrangeBasisDeriv::write_raw(ofstream& fout) const {
  fout << "max_degree = " << max_degree << "\n";

  fout << "my_lagrangeBasisDeriv_middle\n";
  for(int n=1; n<max_degree; n+=2) {
    fout << "  n=" << n << ":";
    for(int i=0; i<=n; i++) {
      fout << "  " << my_lagrangeBasisDeriv_middle[map_odd(n,i)];
    }
    fout << "\n";
  }
  
  fout << "my_lagrangeBasisDeriv_endHi\n";
  for(int n=0; n<max_degree; n++) {
    fout << "  n=" << n << ":";
    for(int i=0; i<=n; i++) {
      fout << "  " << my_lagrangeBasisDeriv_endHi[map(n,i)];
    }
    fout << "\n";
  }
}

int LagrangeBasisDeriv::write_raw(string& filename) const {
  ofstream fout(filename);
  if( !fout ) return 1;
  write_raw(fout);
  return 0;
}
