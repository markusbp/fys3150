#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "test_file.hpp"

// Unit tests to be passed.

TEST_CASE( "Testing", "[Test all]" )
{
    SECTION("Gauss Legendre should integrate P1 with only 2 points exactly from -1 to 1")
    {
      int n = 2;
      double tolerance = 1e-14;
      REQUIRE(test_gauss_legendre(n) < tolerance);
    }

    SECTION("Gauss Laguerre should integrate e^(-x) with only 2 points exactly from 0 to infty")
    {
      int n = 2;
      double tolerance = 1e-14;
      REQUIRE(test_gauss_laguerre(n) < tolerance);
    }

    SECTION("Monte Carlo should integrate a function decently for large n")
    {
      int n = 1e7;
      double tolerance = 1e-3;
      REQUIRE(test_monte_carlo(n) < tolerance);
    }
}
