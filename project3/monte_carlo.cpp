#include <iostream>
#include <cmath>

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

double brute_force(int n, double upper)
{
  double mean = 0;
  double r1, r2, r_12, theta1, theta2, phi1, phi2;
  r1 = 0; r2 = 0; r_12 = 0; theta1 = 0; theta2 = 0; phi1 = 0; phi2 = 0;

  double r1r1 = 0;
  double r2r2 = 0;
  double sinsin = 0;
  double alpha = 0;
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
    mean  += (double) exp(-4*(r1 + r2))/r_12*r1r1*r2r2*sinsin;
  }

  mean = (double) mean/n*jacobi_determinant;
}

double importance(int n)
{
  double r1, r2, r_12, theta1, theta2, phi1, phi2, mean;
  r1 = 0; r2 = 0; r_12 = 0; theta1 = 0; theta2 = 0; phi1 = 0; phi2 = 0, mean = 0;

  double r1r1 = 0;
  double r2r2 = 0;
  double sinsin = 0;
  double alpha = 0;
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
    mean  += (double) 1.0/r_12*r1r1*r2r2*sinsin*exp(-3*(r1 + r2));
  }
  mean = (double) mean/n*jacobi_determinant;
}

int main(int argc, char *argv[])
{
  std::string integral_type = argv[1];
  int n = atoi(argv[2]);
  double integral = 0;
  srand(time(NULL)); // Initialize RNG

  if(integral_type == "brute_force")
  {
    double upper = atof(argv[3]);
    integral = brute_force(n, upper);
    std::cout << "Integral is approximately " << integral << std::endl;
  }
  if(integral_type == "importance")
  {
    integral = importance(n);
    std::cout << "Integral is approximately " << integral << std::endl;
  }
  return 0;
}
