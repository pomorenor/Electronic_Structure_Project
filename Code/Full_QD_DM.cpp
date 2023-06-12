#include <iostream>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>


typedef std::complex<double> complex_number;
typedef complex_number state_type[4];


//FQDDM stands for full
struct FQDDM
{
  double m_gamma;
  double m_Omega;
  double m_omega_0;
  double m_omega_l;

  FQDDM(double gamma = 1.0, double Omega = 1.0, double omega_0 = 1.0): m_gamma(gamma), m_Omega(Omega), m_omega_0(omega_0) {}

  void operator()(const state_type &rho, state_type &drhodt, double t)
  {
    const complex_number I(0.0,1.0);
    drhodt[0] = 2.0*m_gamma*rho[3] -I*m_Omega*();
    drhodt[0] = 0.0;
    drhodt[0] = 0.0;
    drhodt[0] = 0.0;

  }

};





int main(void)
{
  /*
  state_type numeritos ={{1, 2},{3, 4},{1, 2},{3, 4}};

  for(int ii = 0; ii < 4; ++ii){
    std::cout << numeritos[ii] << std::endl;
  }
  */


  return 0;
}
