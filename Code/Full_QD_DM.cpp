#include <iostream>
#include <complex>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

typedef std::complex<double> complex_number;
typedef complex_number state_type[4];

int main(void)
{

  state_type numeritos ={{1, 2},{3, 4},{1, 2},{3, 4}};

  for(int ii = 0; ii < 4; ++ii){
    std::cout << numeritos[ii] << std::endl;
  }



  return 0;
}
