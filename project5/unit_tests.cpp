#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <armadillo>
#include "quantum_dots.hpp"

TEST_CASE( "Testing QD model, for wavefuncs 1 and 2", "[Test QD]" )
{
  double tol = 0.15; // tolerance, acceptance rate
  int n = 1000000;
  double alpha = 1.00;
  double beta = 0.2;
  double omega = 1.00;

  SECTION("Testing psi1 acceptance rate")
  {
    Psi1 trial1(n, alpha, omega, 0);
    trial1.metropolis();
    REQUIRE(fabs((double) trial1.accepted_moves/n - 0.5) < tol);
  }

  SECTION("Testing psi1 acceptance rate")
  {
    Psi2 trial2(n, alpha, beta, omega, 0);
    trial2.metropolis();
    REQUIRE(fabs((double) trial2.accepted_moves/n - 0.5) < tol);
  }
}
