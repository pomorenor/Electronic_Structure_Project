#include <iostream>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>


typedef std::complex<double> state_type;

//QDDM stands for Quantum Dynamics of Density Matrix
// This Functor implements de TD-ODE
struct QDDM
{
  double m_gamma;

  /*
  void operator()()
  {
    const std::complex<double> I(0.0, 1.0);
    std::cout << I << std::endl;
    std::cout << imag(I) <<std::endl;
  }
*/
  QDDM(double gamma = 1.0): m_gamma(gamma) {}

  // This first differential equation is the one associated with rho_ee
  void operator()(const state_type &rho, state_type &drhodt, double t)
  {
    const std::complex<double> I(0.0,1.0);
    drhodt = -2.0*m_gamma*rho;
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

  boost::numeric::odeint::integrate_const(stepper_type(), QDDM(1.0), x, 0.0, 10.0, dt ,streaming_observer(std::cout));

  return 0;
}
