// Ideas: Test ALL configs of 2x2, test init of matrix
// Test energy transition
// expected value should be well-behaved over long mc-cycles.
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <armadillo>
#include "ising_model.hpp"

TEST_CASE( "Testing 2x2 Ising Model init", "[Test Ising init]" )
{
  double tol = 1e-16; // tolerance

  int size = 2; // 2x2 lattice
  int mc_cycles = 1; // toy number
  double temp = 1;  // toy temperature
  IsingModel test_model(size, mc_cycles, temp);
  SECTION("Default spin config is ground state with energy -8")
  {
    REQUIRE(test_model.energy == -8);
  }
  SECTION("One flipped spin should carry energy zero")
  {
    // This also tests that quantities are reset properly :)
    test_model.spins(0,0) = -test_model.spins(0,0); // flip one spin
    test_model.set_ordered_spin_config("custom"); // re-initialize due to custom flip
    REQUIRE(fabs(test_model.energy) < tol);
  }
}

TEST_CASE("Testing 2x2 Ising Model sanity", "[Test Ising run]")
{
  double tol = 1e-2; // tolearnce
  int size = 2;
  int mc_cycles = 1e6;
  double temp = 2.718281828;
  IsingModel test_run(size, mc_cycles, temp);

  SECTION("Mean magnetization should go to zero in 2x2 case")
  {
    test_run.metropolis(0);
    std::cout <<"Mean mag.: " << test_run.mean_mag(mc_cycles) << std::endl;
    REQUIRE(fabs(test_run.mean_mag(mc_cycles-1)) < tol);
  }
}

// TODO: Create more unit tests! Impl. heat-cap, magsep.,
