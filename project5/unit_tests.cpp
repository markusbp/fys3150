#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <armadillo>
#include "quantum_dots.hpp"

TEST_CASE( "Testing QD model, for wavefuncs 1 and 2", "[Test QD]" )
{
  double tol = 0.2; // tolerance, acceptance rate
  int n = 1000000;
  double alpha1 = 0.88;
  double alpha = 1.00;
  double beta = 0.270;
  double omega = 1.00;
  double accept1 = 0;
  double accept2 = 0;

  SECTION("Testing psi1 and psi2 acceptance rate")
  {
    Psi1 trial1(n, alpha1, omega, 0);
    Psi2 trial2(n, alpha, beta, omega, 0);

    for(double omega = 0.01; omega < 1.0; omega += 0.05)
    {
      std::cout << "Omega = " << omega << std::endl;;
      trial1.set_omega(omega);
      trial2.set_omega(omega);
      trial2.metropolis();
      trial1.metropolis();
      accept1 = (double) trial1.accepted_moves/n;
      accept2 = (double) trial2.accepted_moves/n;
      std::cout << "Accept1: " << accept1 << std::endl;
      std::cout << "Accept2: " << accept2 << "\n" <<  std::endl;
      REQUIRE(fabs(accept1 - 0.5) < tol);
      REQUIRE(fabs(accept2 - 0.5) < tol);

    }
  }

}
