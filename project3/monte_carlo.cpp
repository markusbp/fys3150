#include <iostream>
#include <cmath>
#include <chrono>
#include <armadillo>

/*
  Script calculating the expected value of particle separation in Helium atom,
  an integral, using Monte Carlo integration.

  usage:
  ./monte_carlo integral_type n/num_steps upper/stepsize
  - integral_type is a string determining the type of calculation to be done
  can be "brute_force", "importance", "profile_brute_force" or "profile_importance"
  if brute_force or importance
    - calculates the integral using either brute force approach, or importance sampling
    - n is the number of integration points
    - upper is the upper limit of the radial variable random sampling, used only by brute_force
  if profile..
    - Logs time, and integral for various values of n, given by
    - num_steps, the number of different n_values to consider
    - stepsize, the stepsize.
    logs results in csv file

    example:
    ./monte_carlo brute_force 10000 3
    prints: Integral is approximately 0.184826

*/

double random_float(double lower, double upper)
{
  double temp =  (double) rand()/(RAND_MAX);
  return (lower + (upper - lower)*temp);
}

double random_radius()
{
  return (double) rand();
}

double random_angle()
{
  double temp = (double) rand()/(RAND_MAX);
  return 2*M_PI*temp;
}

double brute_force(int n, double upper, double &variance)
{
  double mean = 0;
  double r1, r2, r_12, theta1, theta2, phi1, phi2;
  r1 = 0; r2 = 0; r_12 = 0; theta1 = 0; theta2 = 0; phi1 = 0; phi2 = 0;

  double r1r1 = 0;
  double r2r2 = 0;
  double sinsin = 0;
  double alpha = 0;
  double mean_squared = 0;
  double temp = 0;
  double jacobi_determinant = 4*upper*upper*pow(M_PI, 4);

  for(int i = 0; i<n; i++)
  {
    r1 = random_float(0, upper);
    r2 = random_float(0, upper);
    r1r1 = r1*r1;
    r2r2 = r2*r2;
    theta1 = (double) random_angle()/2.0;
    theta2 = (double) random_angle()/2.0;
    sinsin = sin(theta1)*sin(theta2);
    phi1 = random_angle();
    phi2 = random_angle();

    alpha = cos(theta1)*cos(theta2) + sinsin*cos(phi1 - phi2);
    r_12 = sqrt(r1r1 + r2r2 - 2*r1*r2*alpha);
    temp = (double) exp(-4*(r1 + r2))/r_12*r1r1*r2r2*sinsin;
    mean  += temp;
    mean_squared += temp*temp;
  }

  mean = (double) mean/n*jacobi_determinant;
  variance = (double) mean_squared/n*jacobi_determinant*jacobi_determinant - mean*mean;
  variance = variance/n;
  return mean;
}

double importance(int n, double &variance)
{
  double r1, r2, r_12, theta1, theta2, phi1, phi2, mean;
  r1 = 0; r2 = 0; r_12 = 0; theta1 = 0; theta2 = 0; phi1 = 0; phi2 = 0, mean = 0;

  double r1r1 = 0;
  double r2r2 = 0;
  double sinsin = 0;
  double alpha = 0;
  double temp = 0;
  double mean_squared = 0;
  double jacobi_determinant = 4*pow(M_PI, 4);

  for(int i = 0; i<n; i++)
  {
    r1 = - log( 1 - (double) rand()/RAND_MAX );
    r2 = - log( 1 - (double) rand()/RAND_MAX );
    r1r1 = r1*r1;
    r2r2 = r2*r2;
    theta1 = (double) random_angle()/2.0;
    theta2 = (double) random_angle()/2.0;
    sinsin = sin(theta1)*sin(theta2);
    phi1 = random_angle();
    phi2 = random_angle();

    alpha = cos(theta1)*cos(theta2) + sinsin*cos(phi1 - phi2);
    r_12 = sqrt(r1r1 + r2r2 - 2*r1*r2*alpha);
    temp = (double) 1.0/r_12*r1r1*r2r2*sinsin*exp(-3*(r1 + r2));
    mean += temp;
    mean_squared += temp*temp;
  }
  mean = (double) mean/n*jacobi_determinant;
  variance = (double) mean_squared/n*jacobi_determinant*jacobi_determinant - mean*mean;
  variance = variance/n;
  return mean;
}

int main(int argc, char *argv[])
{
  std::string integral_type = argv[1];
  double integral = 0;
  srand(time(NULL)); // Initialize RNG

  if(integral_type == "brute_force")
  {
    int n = atoi(argv[2]);
    double var = 0;
    double upper = atof(argv[3]);
    integral = brute_force(n, upper, var);
    std::cout << "Integral is approximately " << integral << std::endl;
    std::cout << "Stddev is " << sqrt(var) << std::endl;
  }
  else if(integral_type == "importance")
  {
    int n = atoi(argv[2]);
    double var = 0;
    integral = importance(n, var);
    std::cout << "Integral is approximately " << integral << std::endl;
    std::cout << "Stddev is " << sqrt(var) << std::endl;

  }
  else if(integral_type == "profile_importance")
  {
    int num_steps = atoi(argv[2]);
    int stepsize = atoi(argv[3]);
    double var = 0;

    arma::mat integral = arma::zeros(num_steps, 5);
    arma::vec execution_time = arma::zeros(num_steps);
    arma::vec n_values = arma::zeros(num_steps);
    arma::vec stddev = arma::zeros(num_steps);
    std::cout << "Profiling..." << std::endl;

    double average_time = 0;

    for(int i = 0; i < num_steps; i++)
    {
      for(int j = 0; j < 5; j++)
      {
        auto start = std::chrono::high_resolution_clock::now(); // time execution
        integral(i,j) = importance((i+1)*stepsize, var);
        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = stop - start;
        average_time += elapsed.count();
      }
      execution_time(i) = (double) average_time/5.0;
      n_values(i) = (i+1)*stepsize;
      stddev(i) = sqrt(var);
      std::cout << "n = " << n_values(i) << std::endl;
      average_time = 0;
    }

    integral.save("results/monte_carlo_importance_integral", arma::csv_ascii);
    execution_time.save("results/monte_carlo_importance_timing", arma::csv_ascii);
    n_values.save("results/monte_carlo_n_values", arma::csv_ascii);
    stddev.save("results/monte_carlo_importance_std", arma::csv_ascii);

  }
  else if(integral_type == "profile_brute_force")
  {

    int num_steps = atoi(argv[2]);
    int stepsize = atoi(argv[3]);
    double var = 0;
    double upper = atof(argv[4]);

    arma::mat integral = arma::zeros(num_steps, 5);
    arma::vec execution_time = arma::zeros(num_steps);
    arma::vec n_values = arma::zeros(num_steps);
    arma::vec stddev = arma::zeros(num_steps);
    std::cout << "Profiling..." << std::endl;

    double average_time = 0;

    for(int i = 0; i < num_steps; i++)
    {
      for(int j = 0; j < 5; j++)
      {
        auto start = std::chrono::high_resolution_clock::now(); // time execution
        integral(i,j) = brute_force((i+1)*stepsize, upper, var);
        auto stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = stop - start;
        average_time += elapsed.count();
      }
      execution_time(i) = (double) average_time/5.0;
      n_values(i) = (i+1)*stepsize;
      stddev(i) = sqrt(var);
      std::cout << "n = " << n_values(i) << std::endl;
      average_time = 0;
    }
    integral.save("results/monte_carlo_brute_force_integral", arma::csv_ascii);
    execution_time.save("results/monte_carlo_brute_force_timing", arma::csv_ascii);
    n_values.save("results/monte_carlo_brute_force_n_values", arma::csv_ascii);
    stddev.save("results/monte_carlo_brute_force_std", arma::csv_ascii);
  }
  else
  {
    std::cout << "Error" << std::endl;
  }
  return 0;
}
