#include <iostream>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>


typedef std::complex<double> state_type;

//QDDM stands for Quantum Dynamics of Density Matrix
// This Functor implements de TD-ODE
// This functor is for the rho_ee component TD
struct QDDM_ee
{
  double m_gamma;

  QDDM_ee(double gamma = 1.0): m_gamma(gamma) {}

  // This first differential equation is the one associated with rho_ee
  void operator()(const state_type &rho, state_type &drhodt, double t)
  {
    const std::complex<double> I(0.0,1.0);
    drhodt = -2.0*m_gamma*rho;
  }

};

// This functor is for the rho_eg component TD
struct QDDM_eg
{
  double m_gamma;
  double m_omega_0;

  QDDM_eg(double gamma = 1.0, double omega_0 = 100): m_gamma(gamma), m_omega_0(omega_0) {}

  void operator()(const state_type &rho, state_type &drhodt, double t)
  {
    const std::complex<double> I(0.0, 1.0);
    drhodt = -m_gamma*rho - I*m_omega_0*rho;
  }


};





// We now need something to print
struct streaming_observer
{
    std::ostream& m_out;

    streaming_observer( std::ostream &out ) : m_out( out ) { }

    template< class State >
    void operator()( const State &x , double t ) const
    {
        m_out << t;
        m_out << "\t" << x.real() << "\t" << x.imag() ;
        m_out << "\n";
    }
  };


int main(void)
{
  state_type x = (0.0,1.0);

  const double dt = 0.1;

  typedef boost::numeric::odeint::runge_kutta4< state_type > stepper_type;

  boost::numeric::odeint::integrate_const(stepper_type(), QDDM_eg(1.0), x, 0.0, 10.0, dt ,streaming_observer(std::cout));

  return 0;
}
