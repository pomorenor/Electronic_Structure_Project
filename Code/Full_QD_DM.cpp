#include <iostream>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <vector>
#include <math.h>


typedef std::complex<double> complex_number;
typedef complex_number state_type[4];

double Omega(std::vector<double> dipole_molecule, std::vector<double> E_Field, double alpha, double epsilon);
double gamma(std::vector<double> dipole_molecule, double omega_0);

//FQDDM stands for full
struct FQDDM
{
  double m_gamma;
  double m_Omega;
  double m_omega_0;
  double m_omega_l;

  FQDDM(double gamma = 1.0, double Omega = 1.0, double omega_0 = 1.0, double omega_l = 1.0): m_gamma(gamma), m_Omega(Omega), m_omega_0(omega_0), m_omega_l(omega_l) {}

  void operator()(const state_type &rho, state_type &drhodt, double t)
  {
    const complex_number I(0.0,1.0);
    drhodt[0] = 2.0*m_gamma*rho[3] - I*m_Omega*(rho[2]*std::exp(I*m_omega_l*t) - rho[1]*std::exp(-I*m_omega_l*t));
    drhodt[1] = -m_gamma*rho[1] + I*m_omega_0*rho[1] - I*m_Omega*(rho[3] - rho[0])*std::exp(I*m_omega_l*t) ;
    drhodt[2] = -m_gamma*rho[1] - I*m_omega_0*rho[1] + I*m_Omega*(rho[3] - rho[0])*std::exp(I*m_omega_l*t);
    drhodt[3] = -2.0*m_gamma*rho[3] + I*m_Omega*(rho[2]*std::exp(I*m_omega_l*t) - rho[1]*std::exp(-I*m_omega_l*t));
  }
};


struct observer
{
  std::ostream& m_out;

  observer(std::ostream &out): m_out(out){}
  template <class State>
  void operator()(const State &x, double t) const
  {
    //m_out << t;
    //m_out << "\t" << x[0].real() << "\t" << x[0].imag() << "\t" << x[1].real() << "\t" << x[1].imag() << "\t" << x[2].real() << "\t" << x[2].imag() << "\t" << x[3].real() << "\t" << x[3].imag();
    //m_out << "\t" << std::norm(x[0]) << "\t" << std::norm(x[1]) << "\t" << std::norm(x[2]) << "\t" << std::norm(x[3]);
    //m_out << "\t" << x[0].imag() + x[3].imag();
    m_out << "\t" << x[0].real() << "\t" << x[3].real() << "\t" <<  x[0].real() + x[3].real() ;
    m_out << "\n";
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

  double gamma = 0.62;
  double Omega = 0.01;
  double omega_0 = 0.4;
  double omega_l = 0.2;

  state_type rho = {{1.0/2.0*1.0,0.0},{0.0,0.0},{0.0,0.0},{1.0/2.0*1.0,0.0}};
  const double dt = 0.1;
  typedef boost::numeric::odeint::runge_kutta4< state_type > stepper_type;
  boost::numeric::odeint::integrate_const(stepper_type(), FQDDM(gamma,Omega,omega_0, omega_l), rho, 0.0, 100.0, dt, observer(std::cout));


  return 0;
}


double Omega(std::vector<double> dipole_molecule, std::vector<double> E_Field, double alpha, double epsilon )
{
  double Omega_hat = std::inner_product(dipole_molecule.begin(), dipole_molecule.end(), E_Field.begin(), 0);
  return epsilon*alpha*Omega_hat;
}

double gamma(std::vector<double> dipole_molecule, double omega_0)
{
  double dipole_norm = 0.0;
  double Gamma = 0.0;
  dipole_norm = std::inner_product(dipole_molecule.begin(), dipole_molecule.end(), dipole_molecule.begin(), 0);
  Gamma = dipole_norm*(std::pow(omega_0,3)/(3*M_PI));
  return Gamma;
}
