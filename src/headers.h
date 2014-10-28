using namespace arma;
using namespace Rcpp;
using namespace std;

#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif



struct PExLen{
  mat pexLen;
  mat p;
};


vec do_approx(const mat& x, const vec& y, const vec xo);
vec fitter(const mat& A,  const vec& b);

void nnls_fit(double* a,  int mda,  int m,  int n, double* b,
	      double* x, double* rnorm, double* w, double* zz, int* index,
	      int mode, int nsetp);

vec cpois_fitter(const mat& X, const vec& y, int maxit, double error_limit, int with_se, 
		 int with_res, double Q1, vec& std_err, vec& resd);

void getStart(const List tx, const NumericVector exLen,  const unsigned int iter, mat& p, const unsigned int df);

void fitCr(const mat& X, const mat& Y,  const NumericVector exLen_r, const mat& iTheta, const vec& Len, const int& maxit, 
	   double error_limit, int std_err, const int resd, const double& Q1, const List tx,  const unsigned int iter, 
	   mat& p, const unsigned int df);

mat Zi2V(const mat& Zi, const vec& wt, const uvec& r);

vec ZiVy(const mat& Zi, const vec& y, const vec& wt, const uvec r);

vec Zi_b(const mat& Xf, const mat& B, const  std::vector<std::string> index, 
	 const std::vector<std::string> tx_names, int ncolY);

vec Zb(const mat& Xf, const mat& B, int ncolY);

List Yfun(mat Y, const mat& B, const double Q, const mat& Xf);

List getR(const List tx, int deriv);

vec smspline(vec x, vec y, int df); // => smooth.spline

mat getXf(const mat& X, const mat& iTheta, const vec& len_in, const vec& exLen);

void getZi(cube& Zi_in, const mat& Xf);

mat getLambda(const vec& lambda, const cube& Zi_r, const mat& Xf, const mat& B, const mat& Y, const double Q, 
	      const mat& iTheta, List tx_r, mat X, List R_r, 
	      const std::vector< std::string >& tx_names);

void getsi(const cube& Zi_r, List R_r, const mat& X, const vec& wt, vec& si, int d, const mat& p);

PExLen getPexLen(NumericVector exLen_r, const mat& p, List tx_r, const mat& X);

double integrateNorm(vec args, double lower, double upper);
