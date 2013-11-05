/* ** Functions borrowed by R **

Functions used:
1. approx => R_approx (we have a version based on GSL but it's more convenient to use the R one)
2. smooth.spline => qsbart

*/

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "headers.h"

typedef union {void *p; DL_FUNC fn;} fn_ptr;

DL_FUNC R_ExternalPtrAddrFn(SEXP s){
     fn_ptr tmp;
     tmp.p =  EXTPTR_PTR(s);
     return tmp.fn;
};


////////////////////////////////////////////////////////////////////////////////


/////////////////////
// @@ 1. approx @@ //
/////////////////////

vec approx(const vec x, const vec y, const vec xo)
{

  int method = 1;
  double f = 0;
  unsigned int rule = 2;
  double yleft,yright = NA_REAL;
  
  int nout = xo.n_elem;
  int nxy = x.n_elem;
  
  // x and y must be sorted (ascending) according to x.
  uvec index = sort_index(x);

  vec sx(x);
  vec sy(y);

  sx = sort(sx);
  sy = sy.elem(index);

  if(rule == 2)
    {
      yleft = sy[0];
      yright = sy[(sy.n_elem-1)];
    }

  vec yo(xo);
  
  static void (*R_approx)(double *sx, double *sy, Sint *nxy, double *xout, Sint *nout,
			  Sint *method, double *yleft, double *yright, double *f) = NULL;
  
  if (!R_approx) {
    Function getNativeSymbolInfo("getNativeSymbolInfo");
    List nativeSymbolInfo = getNativeSymbolInfo("R_approx");
    SEXP xp = nativeSymbolInfo["address"];
    DL_FUNC R_approx_p = R_ExternalPtrAddrFn( xp ) ;
    
    R_approx = (void(*)(double *, double *, Sint *, double *, Sint *,
			Sint *, double *, double *, double *)) R_approx_p;
  }
  (*R_approx)(const_cast<double *>(sx.begin()),
	      const_cast<double *>(sy.begin()), &nxy, yo.begin(), &nout, &method, &yleft,  &yright, 
	      &f);
  
  return yo;

}

vec do_approx(const mat& x, const vec& y, const vec xo)
{
  vec out = zeros<vec>(x.n_rows);

  for(unsigned int nr = 0; nr < x.n_rows; nr ++)
    {
      vec iX = conv_to<vec>::from(x.row(nr));
      out(nr) = as_scalar(approx(iX, y, xo)); // xo is actually a scalar
    }

  return out;
}




////////////////////////////
// @@ 2. smooth_spline @@ //
////////////////////////////

// // This is not really the same story as in R, where we have signif(x,6L) before sort & unique


// Source: http://stackoverflow.com/questions/7535299/how-to-format-numbers-to-significant-digits-using-stl
double toPrecision(double num, int n) {

    if(num == 0) {
      return 0;
    }

    double d = std::ceil(std::log10(num < 0 ? -num : num));
    int power = n - (int)d;
    double magnitude = std::pow(10., power);
    double shifted = ::round(num*magnitude);

    return shifted/magnitude;
}

vec unisort(const vec x_in, const int  digits)
{

  vec sx_in(x_in);
  sx_in = sort(x_in);
  

  std::vector<double> input = conv_to<std::vector<double> >::from(sx_in);
  
  std::vector<double> output(input.size(),0); // initialize at zeros

  std::vector<int> digvec(input.size(),digits);

  std::transform(input.begin(),input.end(),digvec.begin(), output.begin(), toPrecision);

  std::vector<double>::iterator it;

  it = unique(output.begin(),output.end());

  output.resize(it - output.begin());

  vec out = conv_to<vec>::from(output);

  return out;

}

rowvec fun1(const uvec udx, const vec y, const vec w)
{
  rowvec out(3);
  vec tmp;
  out(0) = accu(w.elem(udx));
  tmp = w.elem(udx) % y.elem(udx);
  out(1) = accu(tmp);
  tmp = w.elem(udx) % y.elem(udx) % y.elem(udx);
  out(2) = accu(tmp);

  return out;
}

int n_knots(int n)
{
  if(n < 50)
    return(n);

  double a1 = log2(50);
  double a2 = log2(100);
  double a3 = log2(140);
  double a4 = log2(200);

  double ans;
  
  if(n < 200){
    ans = pow(2, a1 + (a2 - a1) * (n - 50)/150);
    return ans;}
  
  if(n < 800){
    ans = pow(2,a2 + (a3 - a2) * (n - 200)/600);
    return ans;}
  
  if(n < 3200){
    ans = pow(2,a3 + (a4 - a3) * (n - 800)/2400);
    return ans;}

  double bb = n - 3200;
  ans = 200 + pow(bb,0.2);
  
  return (int)ans;

}


vec smspline(const vec x, const vec y, const int df)
{
    int ldnk = 1;
    int ld4 = 4;
    vec contr_sp(5);
    contr_sp(0) = -1.5; contr_sp(1) =  1.5; contr_sp(2) = 0.0001; contr_sp(3) = 0.00000002; contr_sp(4) = 500;
  
    
    int n = x.n_elem;
    vec w = ones<vec>(n); // we don't use weights

    vec ux = unisort(x,6);
    int nx = ux.n_elem;


    uvec ox;
    mat tmp;

    if(nx == n)
      {
	// ox = linspace(0,n-1,n); // not needed at the moment
	tmp.insert_cols(0,w);
	tmp.insert_cols(1,w % y);
	tmp.insert_cols(2,w % y % y);

      } else
      {
	uvec tox(n);
	
	for(unsigned int i=0; i < ux.n_elem; i++)
	  {
	    uvec idx = find(x == ux(i));
	    uvec vals(idx.n_elem);
	    vals.fill(i);
	    tox.elem(idx) = vals;
	    tmp.insert_rows(i,fun1(idx,y,w));
	  }
	
	ox = tox;
      }


    vec wbar = tmp.col(0);
    vec ybar = tmp.col(1);

    ybar.elem(find(wbar > 0)) = ybar.elem(find(wbar > 0))/wbar.elem(find(wbar > 0));

    vec tvec = tmp.col(2) - (wbar % ybar % ybar);
    double yssw = accu(tvec);
    
    double rux = ux(ux.n_elem - 1) - ux(0);
    vec xbar = (ux - ux[0])/rux;

    int nknots = n_knots(nx);
    int nk = nknots + 2;
    
    colvec knot(3);
    knot.fill(xbar(0));

    uvec idx = conv_to<uvec>::from(linspace(0,nx-1,nknots));
    knot = join_cols(knot,xbar.elem(idx));
    
    vec tvec1(3);
    tvec1.fill(xbar(xbar.n_elem - 1));
    knot = join_cols(knot,tvec1);

    double spar = 0;

    double dofoff = 0;
    int icrit = 1;

    if(df > 1 && df <= nx)
      {
	icrit = 3;
	dofoff = (double) df;
      } else
      {
	cout << "you must supply 1 < df <= n,  n = #{unique x} = " << nx << endl;
	cout << "df set to " << nx << endl;
	dofoff = (double) nx;
      }

    ivec iparms(3);
    iparms(0) = icrit;
    iparms(1) = 0;
    iparms(2) = (int) contr_sp(4);

    double penalty = 1;

    uvec indices;
    indices << 0 << 1 << 2 << 3;
    vec parms = contr_sp.elem(indices);
    
    vec scratch(17 * nk + 1);
    vec lev(nx);
    vec coef(nk);
    vec ty(nx);
    ivec ier(1);
    vec crit(1);

    static void (*rbart)(double *penalty, double *dofoff, 
    			  double *xbar, double *ybar,double *wbar, 
    			  double *yssw, Sint *nx, double *knot, 
    			  Sint *nk, double *coef, double *ty, double *lev, 
    			  double *crit, Sint *iparms, double *spar, double *parms , 
    			  double *scratch, Sint *ld4, Sint *ldnk, Sint *ier) = NULL;
    
    if (!rbart) {
      Function getNativeSymbolInfo("getNativeSymbolInfo");
      List nativeSymbolInfo = getNativeSymbolInfo("rbart");
      SEXP xp = nativeSymbolInfo["address"];
      DL_FUNC rbart_p = R_ExternalPtrAddrFn(xp) ;
      
      rbart = (void(*) (double *, double *, 
    			 double *, double *, double *, 
    			 double *, Sint *, double *, 
    			 Sint *, double *, double *, double *, 
    			 double *, Sint *, double *, double *, 
    			 double *, Sint *, Sint *, Sint *)) rbart_p;
    }
    (*rbart)(&penalty, &dofoff,
    	      const_cast<double *>(xbar.begin()),  const_cast<double *>(ybar.begin()), 
	      const_cast<double *>(wbar.begin()), 
    	      &yssw, &nx, const_cast<double *>(knot.begin()), &nk, coef.begin(),  
    	      ty.begin(), lev.begin(), const_cast<double *>(crit.begin()), const_cast<int *>(iparms.begin()), 
	      &spar,  const_cast<double *>(parms.begin()),
    	      scratch.begin(), &ld4, &ldnk, ier.begin());
   
    vec fitted;

    if(nx == n)
	fitted = ty;
    else
      fitted = ty.elem(ox);
   
    return fitted;

}


///////////////////
// @@ 3. nnls @@ //
///////////////////

void nnls_fit(double* a,  int mda,  int m,  int n, double* b,
	      double* x, double* rnorm, double* w, double* zz, int* index,
	      int mode, int nsetp)
{
  
  
  static void (*nnls)(double* a,  Sint* mda,  Sint* m,  Sint* n, double* b,
		      double* x, double* rnorm, double* w, double* zz, Sint* index,
		      Sint* mode, Sint* nsetp) = NULL;
  
  if (!nnls) {
    Function getNativeSymbolInfo("getNativeSymbolInfo");
    List nativeSymbolInfo = getNativeSymbolInfo("nnls");
    SEXP xp = nativeSymbolInfo["address"];
    DL_FUNC nnls_p = R_ExternalPtrAddrFn( xp ) ;
    
    nnls = (void(*)(double *,  Sint *,  Sint *,  Sint *, double *,
		    double *, double *, double *, double *, Sint *,
		    Sint *, Sint *)) nnls_p;
  }

  (*nnls)(a, &mda,  &m,  &n, b, x, rnorm, w, zz, index, &mode, &nsetp);
  
}

// // [[Rcpp::export]]
// SEXP smsplineR(SEXP in_x, SEXP in_y, SEXP in_df)
// {
//   NumericVector rx(in_x);
//   NumericVector ry(in_y);
//   NumericVector rdf(in_df);

//   vec x = as<vec>(rx);
//   vec y = as<vec>(ry);

//   int df = as<int>(rdf);
  
//   vec ans;

//   ans = smspline(x,y,df);

//   return(wrap(ans));
// }


// INTEGRATION //

typedef void integr_fn(double *x, int n, void *ex);

void plen(double *x, int n, void *ex)
{
  double *v;
  v = (double*)ex;

  double y1 = v[0]; // not used by plen1
  double y2 = v[1];
  double rlen = v[2];
  double mulen = v[3];
  double sdlen = v[4];
  int wplen = v[5]; // 1= plen1, 2=plen2

    if(wplen == 1)
      {
	for(int i=0; i<n; i++) {
	  x[i] = R::pnorm(y2-x[i]+rlen,mulen,sdlen,1,0)- R::pnorm(2*rlen,mulen,sdlen,1,0);
	}
      }
    else
      {
	for(int i=0; i<n; i++) {
	  x[i] = R::pnorm(y2-x[i]+rlen,mulen,sdlen,1,0) - R::pnorm(y1-x[i]+rlen,mulen,sdlen,1,0);
	}
      }
}


double integrateNorm(vec args, double lower, double upper)
{


  static void (*Rdqags)(integr_fn f, void *ex, double *a, double *b,
			double *epsabs, double *epsrel,
			double *result, double *abserr, int *neval, int *ier,
			int *limit, int *lenw, int *last, int *iwork, double *work) = NULL;
  
  Function getNativeSymbolInfo("getNativeSymbolInfo");
  List nativeSymbolInfo = getNativeSymbolInfo("Rdqags");
  SEXP xp = nativeSymbolInfo["address"];
  DL_FUNC Rdqags_p = R_ExternalPtrAddrFn( xp ) ;
  
  Rdqags = (void(*)(integr_fn *, void *, double *, double *,
		    double *, double *,double *, double *, int *, int *,
		    int *, int *, int *, int *, double *)) Rdqags_p;
  
  

  if((args(5) != 1) & (args(5) != 2))
    {
      throw range_error("Last value must be in either 1 or 2");
      return 0;
    }
  
  if(args.n_elem != 6)
    {
      throw range_error("first argument must be vector of length 6");
      return 0;
    }
  
  if(upper < lower)
    {
      throw range_error("Check the limits of integration (must be lower < upper)");
      return 0;
    }

  if((args(5) == 2) & (args(0) >= args(1)))
    {


      printf("ylim1=%5.0f/ylim2=%5.0f\n",args(0),args(1));

      throw range_error("Check the y limits (must be y1 < y2)");
      return 0;
    }
  
  if((args(5) == 2) & (args(0) < upper))
    {
      throw range_error("upper limit of integration must be smaller than y1");
      return 0;
    }
  

  
  void *ex;
  ex = args.begin();
  
  double epsabs = 0.0001;
  double epsrel = 0.0001;
  double result = 0.0;
  double abserr = 0.0;
  int neval = 0;
  int ier = 0;
  int limit = 100;
  int lenw = 4*limit;
  int last = 0;
  int iwork = 100;
  double work = 400.0;
  
  (*Rdqags)(plen, ex, &lower, &upper, &epsabs, &epsrel, &result, &abserr,
	    &neval, &ier, &limit, &lenw, &last, &iwork, &work);
  
  return(result);
}
