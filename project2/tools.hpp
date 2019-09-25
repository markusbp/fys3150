void write_to_file(arma::mat &matrix, std::string filename, int n);
double test_eigenval_decomp(arma::mat &sym_matrix, arma::vec &eigenvals,
                          arma::mat &eigenvecs, double a, double d, int n);
void test_max_nondiagonal(int n);
double test_norm(double tolerance, int n);
