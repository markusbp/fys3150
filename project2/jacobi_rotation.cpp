#include <armadillo>
#include <math.h>

void max_nondiagonal(arma::mat &matrix, int (&max_index)[2], double &max_val ,int n)
{
  // Find maximum non-diagonal element of matrix matrix
  // Returns indices of max element by setting values of array max_index
  double current = 0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (i == j)
      {
        continue; // skip diagonal elements
      }
      current = fabs(matrix(i, j)); // absolute value of current element
      if (fabs(max_val) < current) // compare to see if it is greater than max
      {
        max_val = matrix(i,j);
        max_index[0] = i;
        max_index[1] = j;
      }
    }
  }
}

void jacobi_rotation(arma::mat &sym_matrix, double tolerance,  int n)
{
  double max_val = tolerance*10; // Set to some value which has to be overwritten
  int max_index[2] = {-1, -1}; // Initial indices, garbage value
  double c, cc, s, ss, t, tau; // cos, sin, tan, and tau parameter
  double aik, ail, akk, all; // Placeholders for elements to be updated
  aik = 0; ail = 0; akk = 0; all = 0;
  int k, l, count;
  k = 0; l = 0; count = 0; c = 0; cc = 0; s = 0; ss = 0; t = 0;
  //std::cout << "Diagonalization started..." << std::endl;
  while(fabs(max_val) > tolerance)
  {
    // Initialize to zero to guarantee finding larger element, but still execute loop
    max_val = 0;
    max_nondiagonal(sym_matrix, max_index, max_val, n); // find current largest element
    k = max_index[0]; // row, column of largest element
    l = max_index[1];
    tau = (sym_matrix(l, l) - sym_matrix(k, k))/(2.0*sym_matrix(k, l));

    if(tau > 0)
    {
      t=1/(tau+sqrt(1+tau*tau));
    }
    else
    {
      t=-1/(-tau+sqrt(1+tau*tau));
    }
    c = 1.0/sqrt(1 + t*t);
    s = t*c;
    cc = c*c;
    ss = s*s;
    for (int i = 0; i < n; i++)
    {
      if ((i==k) || (i==l)){continue;}
      aik = sym_matrix(i, k); ail = sym_matrix(i, l);
      akk = sym_matrix(k, k); all = sym_matrix(l, l);
      sym_matrix(i, k) = c*aik - s*ail;
      sym_matrix(i, l) = c*ail + s*aik;
      sym_matrix(k, i) = sym_matrix(i, k); sym_matrix(l, i) = sym_matrix(i, l);
    }
    sym_matrix(k, k) = cc*akk - 2*sym_matrix(k, l)*c*s + ss*all;
    sym_matrix(l, l) = cc*all + 2*sym_matrix(k, l)*c*s + ss*akk;
    sym_matrix(k,l)  = 0; sym_matrix(l, k) = sym_matrix(k, l);
    count++; // Keep track of the number of iterations
  }
  //std::cout << "Diagonalization finished in " << count << " iterations." << std::endl;
}
