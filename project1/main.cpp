#include <iostream>   // input and output
#include <cmath>      // math functions
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include "../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/lib.h"
std::ofstream ofile;

// Implementation of Gaussian elimination of tridiagonal matrix
// Runtime arguments : output_filename.txt n algorithm_type
// output_filename is name of output file for storing results
// n is number of calculation steps
// algorithm_type is the chosen algorithm, can be either "general" or "simple"
// Runtime example for code:
/*
Execution time for simple algorithm: 5.401e-06 s for n=100
Maximum relative error (log10): -3.08803683155941
Execution time, LU decomposition: 0.001802405
Average absolute difference between LU and Gaussian elimination methods: 4.58690377347359e-15
*/


void general_elimination(double *a, double *b, double *c, double *bprime, double *fprime, double *f, double *u, int n)
{
    // Function to compute forward and backward substitution for general tridiagonal matrix
    double a_bprime; // Extra variable, to avoid calculating a[i]/bprime[i-1] twice in each iteration
    // Forward substitution
    for (int i = 1; i <= n-1; ++i)
    {
        a_bprime = a[i]/bprime[i-1];
        fprime[i] = f[i] - fprime[i-1]*a_bprime;
        bprime[i] = b[i] - c[i-1]*a_bprime;
    }
    // Backward substitution
    u[n-1] = fprime[n-1]/bprime[n-1]; // Due to boundary conditions
    for (int i=n-2; i >= 0; i--)
    {
        u[i] = (fprime[i] - c[i]*u[i+1])/bprime[i];
    }
}

void simple_elimination(double *bprime, double *fprime, double *f, double *u, int n)
{
    // Forward substitution
    for (int i = 2; i <= n; ++i)
    {
        fprime[i-1] = f[i-1] + fprime[i-2]/bprime[i-2];
        bprime[i-1] = (i+1.0)/i;
    }
    u[n-1] = fprime[n-1]/bprime[n-1]; // Boundary condition
    // Backward substitution
    for (int i=n-2; i >= 0; i--)
    {
        u[i] = (fprime[i] + u[i+1])/bprime[i];
    }
}

void initialize_variables(double *a,double *b, double *c,double *u, double *f,double *fprime, double *bprime, double *solution,  int n)
{
    /*
     Initialize all variables related to simple and general algorithms.
     Does not initialize variables associated with LU decomp, to save memory.
     Takes in pointers to all required variables, initializes them in place.
    */
    double x, h;
    h = 1.0/(n+1.0); // stepsize
    for (int i = 0; i < n ; ++i)
    {
        a[i] = -1.0;    // elements below diagonal
        b[i] = 2.0;     // diagonal elements
        c[i] = -1.0;    // elements above diagonal
        u[i] = 0.0;     // approximation/numerical solution
        x = h*(i+1);    // x-coordinate
        fprime[i] = 0.0;    // temporary values during substitution
        bprime[i] = 0.0;    // variable values during substitution
        f[i] = 100.0*exp(-10.0*x)*h*h;  // source term
        solution[i] = 1.0 - (1.0 - exp(-10.0))*x  - exp(-10.0*x); // analytical solution
    }
}

void initialize_lu_variables(double **A, int *rowswaps, int n)
{
    // Initialize matrix to be lower upper decomposed
    for (int i= 0; i < n; i++)
    {
        rowswaps[i] = 0;
        for (int i = 0; i < n; i++)
            A[i] = new double[n];
    }
    for (int i = 1 ; i<n-1 ; i++)
    {
        for (int j = 0; j <n-1; j++)
        {
            A[i][j] = 0; // Make sure all elements are initialized to zero
        }
        A[0][i] = 0; // Initialize tridiagonal
        A[i][i] = 2;
        A[i][i-1] = -1;
        A[i][i+1] = -1;
    }
    A[0][0] = 2;
    A[0][1] = -1;
    A[n-1][n-1] = 2;
    A[n-1][n-2] = -1;
}

void delete_lu_variables(double **A, int *rowswaps, int n)
{
    // Clear variables associated with LU decomposition
    delete [] rowswaps;
    for (int i = 0; i < n; i++)
        delete[] A[i];
    delete[] A;
}

void delete_variables(double *a, double *b, double *c,double *u, double *f,double *fprime, double *bprime, double *solution)
{
    // Clear variables
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] bprime;
    delete [] fprime;
    delete [] solution;
    delete [] f;
    delete [] u;
}

int main (int argc, char *argv[])
{
    std::string output_file = argv[1];
    std::string algorithm_type = argv[3];
    int n = std::atoi(argv[2]);
    ofile.open(output_file);// std::ios_base::app);
    // Initialize all arrays/pointers
    double *a, *b, *c, *bprime, *fprime, *f, *u, *solution, *relative_error, max_relative_error, error_diff;
    // Allocate memory
    a = new double[n];
    b = new double[n];
    c = new double[n];
    f = new double[n];
    u = new double[n];
    bprime = new double[n];
    fprime = new double[n];
    solution = new double[n];
    relative_error = new double[n];

    max_relative_error = 0;
    error_diff = 0;

    // Initialize constants and variables with values
    initialize_variables(a,b,c,u,f,fprime, bprime, solution, n);

    // Forward substitution
    fprime[0] = f[0];
    bprime[0] = b[0];
    auto start = std::chrono::high_resolution_clock::now(); // time execution
    // Perform calculation, either using simple or general algorithm, chosen from command line
    if (algorithm_type == "general")
    {
        general_elimination(a, b, c, bprime, fprime, f, u, n);
    }
    else if (algorithm_type == "simple")
    {
        simple_elimination(bprime, fprime, f, u, n);
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout <<  "Execution time for " << algorithm_type << " algorithm: "
              << elapsed.count() << " s for n=" << n << std::endl;

    std::string sol_val;
    for (int i= 0;i < n;i++)
    {
        relative_error[i] = log10(fabs((u[i] - solution[i])/solution[i]));
        max_relative_error = (i == 0) ? relative_error[i] : std::max(max_relative_error, relative_error[i]);
        sol_val = std::to_string(u[i]);
        ofile << sol_val << std::endl;
    }
    std::cout<< "Maximum relative error (log10): " << std::setprecision(15) <<max_relative_error << std::endl;

    // Test LU decomposition from course library
    double **A; // Initialize full tridiagonal matrix
    int *rowswaps;
    double evenodd;
    rowswaps = new int[n]; // Keeps track of operations in LU decomp
    A = new double*[n];

    initialize_lu_variables(A, rowswaps, n); // setup all variables
    auto startlu = std::chrono::high_resolution_clock::now(); // time it
    ludcmp(A, n, rowswaps, &evenodd);
    lubksb(A, n, rowswaps, f); // changes f in place
    auto finishlu = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedlu = finishlu - startlu;
    std::cout <<  "Execution time, LU decomposition: " << elapsedlu.count() << std::endl;
    delete_lu_variables(A, rowswaps, n);
    // Compute average absolute difference between solutions
    for (int i = 0;i < n; i++)
    {
        error_diff += fabs(f[i] - u[i]);
    }
    error_diff = error_diff/double(n); // average value
    std::cout << "Average absolute difference between LU and Gaussian elimination methods: " << error_diff << std:: endl;
    // Clear memory
    delete_variables(a,b,c,u, f,fprime, bprime, solution);
    ofile.close();
    return 0;
}


