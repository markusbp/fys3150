#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <armadillo>
#include "tools.hpp"
#include "jacobi_rotation.hpp"

TEST_CASE( "Testing", "[Test all]" )
{
    double tolerance = 1e-10;
    int n = 30;   // Test for a 30x30 matrix case
    SECTION("Testing set up of tridiagonal matrix")
    {
    REQUIRE(test_tridiagonal(n) < tolerance);
    }

    SECTION("Check that max_nondiagonal can find maximum matrix element vs. arma method")
    {
      REQUIRE(test_max_nondiagonal(n) < tolerance);
    }
    SECTION("Arma eigenvalues should match Jacobi's method Eigenvalues.")
    {
      REQUIRE(test_eigenvals(tolerance, n) < tolerance);
    }
    SECTION("Frobenius norm of matrix should be preserved when subject to Jacobi rotation")
    {
      REQUIRE(test_norm(tolerance, n) < tolerance);
    }
}
