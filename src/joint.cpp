#include <RcppArmadillo.h>
#include <R_ext/Utils.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "ranker.h"
#include "headers.h"

#include <fstream> // to write to file: for debugging only
#include <stdio.h>


// REMOVE THIS
// #include <ctime>
// #include <time.h>
// std::clock_t start = std::clock();
// std::cout<< "5: Update-" << it <<": " << ( ( std::clock() - start ) / (double)CLOCKS_PER_SEC ) << endl;
// start = std::clock();



vec fitter(const mat& A,  const vec& b)
{
  
  int m = A.n_rows;
  int n = A.n_cols;

  double *a = new double[m*n];
  double *w = new double[n];
  int *indx = new int[n];
  double *x = new double[n];
  
  for (int i=0; i < n; i++)
    {
      for(int j=0;j<m;j++)
        {
	  a[i*m+j] = A(j,i);
        }

      w[i] = 0.0;
      indx[i] = 0;
      x[i] = 0.0;
    }
  
  double *rhs = new double[m];
  double *zz = new double[m];
  for (int j=0;j<m;j++)
    {
      rhs[j] = b(j);
      zz[j] = 0.0;
    }

  double rnorm = 0;
  
  int mode = 0;
  int nsetp = 0;
  
  nnls_fit(a, m, m, n, rhs, x, &rnorm, w, zz, indx, mode, nsetp);
  vec sol(n);
   
  for (int i = 0; i < n; i++) sol(i) = x[i];
  
  delete[] a;
  delete[] rhs;
  delete[] x;
  delete[] w;
  delete[] zz;
  delete[] indx;

  return sol;
}


////////////////////////////////////////////////////////////////////////////////
// NOTE: X must have exons on rows and tx on columns, so t(design)!!!
// Tested on 21/7/2011 => OK

vec cpois_fitter(const mat& X, const vec& y, int maxit, double error_limit, int with_se, 
		  int with_res, double Q1, vec& std_err, vec& resd)
{
  
  int nc = X.n_cols;

  vec betas = zeros<vec>(nc);
  vec old(betas);
  vec secureBetas(betas);

  betas = fitter(X,y);

  vec xb = (y), mu(y);
  vec res = zeros<vec>(y.n_elem);
  vec wt = ones<vec>(y.n_elem);

  vec ares(res);
  vec nwt(wt);
  vec iY(y);

  mat WT = ones<mat>(y.n_elem,X.n_cols);
  mat iX = zeros<mat>(X.n_rows,X.n_cols);

  vec Err = zeros<vec>(nc);
  vec aB = zeros<vec>(nc);

  double q3, mErr;

  for(int it=0; it < maxit; it++)
    {

      xb = X*betas;
      mu = xb;
     
      res = y - xb;

      for(unsigned int i=0; i < mu.n_elem; i++)
	{
	  if(mu(i) < 0.1) mu(i) = 0.1; // with sugar we could use ifelse
	}
      
      wt = 1/mu;

      ares = abs(res);
      // std::transform(res.begin(),res.end(),ares.begin(), absolute);


      // CHECK: 1) need to input a vec but we have a double here => cast
      // 2) What about NA?

      // q3 = quantile(const_cast<double *>(ares.begin()), ares.n_elem, Q1);
      q3 = quantile(conv_to< std::vector<double> >::from(ares),Q1);

      nwt = (q3/(ares + 0.01))/mu;

      for(unsigned int i=0; i < wt.n_elem; i++)
	{
	  if(ares(i) >= q3) {wt(i) = nwt(i);}
	  if(wt(i) < 0.001) {wt(i) = 0.001;}
	  if(wt(i) > 10) {wt(i) = 10;}
	}

      old = betas;
      
      for(unsigned int i = 0; i < X.n_cols; i++)
	{WT.col(i) = wt;}


      iX = trans(X % WT) * X;
      iY = trans(X % WT) * y;

      betas = fitter(iX,iY);


      secureBetas = betas;
      
      for(unsigned int i = 0; i < secureBetas.n_elem; i++) 
	{ 
	  if(secureBetas(i) < 0.001) {secureBetas(i) = 0.001;}
	}

      
      Err = abs(betas - old);
      aB = abs(secureBetas);
      
      mErr = as_scalar(max(Err/aB));


      if(mErr < error_limit)
	break;
      
    }



  if(with_se == 1)
    {
      mat Iinv;
      vec wt2(wt.n_elem);
      wt2 = wt % wt % res % res;
      
      mat WT2(wt.n_elem,X.n_cols);
      
      for(unsigned int i = 0; i < X.n_cols; i++)
	{WT2.col(i) = wt2;}
      
      mat J = trans(X % WT2) * X;
      mat iX = trans(X % WT) * X; // we should be able to use the previously computed one. Check it out.

      mat Imat(X.n_cols,X.n_cols);
      Imat.eye(); // identity matrix

      for(unsigned int i=0; i < Imat.n_cols; i++)
      	{
	  vec tmp = conv_to<vec>::from(Imat.row(i));
      	  rowvec otmp = conv_to< rowvec >::from(fitter(iX,tmp)); // problema con vettore estratto!
      	  Iinv.insert_rows(i,otmp);
      	}

      mat cov_matrix = Iinv * J * Iinv;
      std_err = sqrt(diagvec(cov_matrix));


    }

  if(with_res == 1)
    {
      resd = res;
    }


  return betas;
  
}




////////////////////////////////////////////////////////////////////////////////
// Tested on  Aug 2011 => OK

RcppExport SEXP fitJoint_R(SEXP iTheta_in, SEXP maxit_in, SEXP  error_limit_in, SEXP X_in, 
			   SEXP Y_in,SEXP len_in,SEXP std_err_in, SEXP resd_in, SEXP lambda_in,
			   SEXP d_in, SEXP Q1_in, SEXP DF_in, SEXP exLen_in, 
			   SEXP tx_in, SEXP use_conditional_in,SEXP ridge_lambda_in)
			   // , SEXP iGene)
{
  BEGIN_RCPP


  const NumericMatrix X_r(X_in);
  const mat X(X_r.begin(),X_r.nrow(),X_r.ncol(),true);


  const NumericMatrix Y_r(Y_in);
  mat Y(Y_r.begin(),Y_r.nrow(),Y_r.ncol(),true);

  const NumericVector exLen_r(clone(exLen_in));
  const vec exLen(exLen_r.begin(),exLen_r.size(),true);

  const NumericVector lambda_r(lambda_in);
  const vec lambda(lambda_r.begin(),lambda_r.size(),true);

  const NumericVector len_r(len_in);
  const vec len(len_r.begin(),len_r.size(),true);


  const List tx_r(clone(tx_in));

  const NumericMatrix iTheta_r(iTheta_in);
  List dimnames = iTheta_r.attr("dimnames");
  mat iTheta(iTheta_r.begin(),iTheta_r.nrow(),iTheta_r.ncol(),true); // make a copy

  vec DF_v = as<vec>(DF_in); // needed for approxR...to lazy to adapt to a scalar...

  const unsigned int ncY = Y.n_cols;
  const unsigned int nrY = Y.n_rows;
  const unsigned int nEx = X.n_cols;
  const unsigned int nTx = X.n_rows;
  
  mat iLen(nTx,nEx);

  for(unsigned i=0;i < iLen.n_rows; i++)
    {
      iLen.row(i) = conv_to<rowvec>::from(exLen);
    }

  // scalars

  const int d = as<int>(d_in);
  const int DF = as<int>(DF_in);
  const int maxit = as<int>(maxit_in);
  const double err = as<double>(error_limit_in);
  const unsigned int std_err = as<unsigned int>(std_err_in);
  const int maxit_fit = 100;       // cpois maxit set to 100

  const unsigned int resd = as<unsigned int>(resd_in);
  const double Q1 = as<double>(Q1_in);

  /////////////////////////////////////////////

  vec ym = sum(Y,1);
  vec obsSum = ym;

  std::vector< std::string > tx_names = tx_r.attr( "names" ) ;

  List R_r = getR(tx_r,d);

  // ncY = N samples
  // nEx = N exons
  // nTx = N tx


  bool fLambda = lambda.n_elem == 1; // Lambda is fixed, so no df (dfall) estimation

  vec se_out;
  vec res_out;

  if(resd == 0)
    res_out.fill(0);

  if(std_err == 0)
    se_out.fill(0);

  int iter = 0;

  for(iter=0; iter < maxit; iter++)
    {

      ////////////////////////////////////////////////////////////////////////////////////////////////////

      mat p = zeros<mat>(nTx,nEx);
      vec si= zeros<vec>(nTx);

      double checkError1;
      double checkError2;
      
      cube Zi_r(ncY, nrY, nTx); // <n samples>-x-<n exons>-x-<n tx>
      
      ////////////////////////////////////////////////////////////////////////////////////////////////////


      // estimate p
      fitCr(X, Y, exLen_r, iTheta, len, maxit_fit, err, std_err, resd, Q1,tx_r,iter, p , DF);

      // get Xf
      const mat Xf = getXf(X, iTheta, len, exLen);

      // get Zi matrix
      getZi(Zi_r, Xf);

      if(!fLambda)
	{
	  mat dfall = getLambda(lambda, Zi_r, Xf, p, Y, Q1, iTheta, tx_r, X, R_r,tx_names);
	  si = 1/do_approx(dfall,lambda,DF_v);

	}
      else
	{
	  si.fill(1/as_scalar(lambda));
	}

      
      vec val = zeros<vec>(iTheta.n_rows);
      uvec sel = sum(iTheta,1) == val; // there might be some precision issue here!!!

      if(sum(sel))
      	{
      	  vec tt(sum(sel));
      	  tt.fill(0.05);
      	  si.elem(find(sel > 0)) = tt;
      	}

      for(int it=0; it < 10; it++) // updates p matrix 10 times
	{


	  List fitY = Yfun(Y,p,Q1,Xf);
	  vec yc = as<vec>(fitY["Y"]);
	  vec wt = as<vec>(fitY["wt"]);

	  mat hold_p = p; // copy p to hold it during the loop along tx. p will be updated
	  
	  for(unsigned int idx = 0; idx < nTx; idx ++) // loops along tx
	    {

	      std::vector<std::string> index; // hold the idx transcript name. This could be simply a char
	      index.push_back(tx_names[idx]);

	      mat Zi = Zi_r.slice(idx);

	      const StringVector itx = tx_r[idx]; // store the exons names for the idx transcript
	      
	      bool allzero = sum(iTheta.row(idx)) == 0;

	      if(!(itx.size() == 1 || allzero))
		{
		  
		  vec iEx = conv_to< vec>::from(X.row(idx));
		  uvec isone = find(iEx != 0);
		  mat iR = R_r[idx];

		  vec iP = conv_to< vec>::from(hold_p.row(idx));
		  
		  if(Xf.n_cols > 1)
		    {
		      yc = as<vec>(fitY["Y"]) - Zi_b(Xf,hold_p,index,tx_names,ncY);
		    }


		  vec iZiVy = ZiVy(Zi,yc, wt, isone);


		  mat iZi2 = Zi2V(Zi,wt,isone) + (iR / as_scalar(si[idx]));


		  vec iFit = fitter(iZi2,iZiVy);

		  iP.elem(isone) = iFit;
		  p.row(idx) = conv_to<rowvec>::from(iP);
		}
	    } // end iteration along tx (idx)
	  


	  if(DF == 0)
	    getsi(Zi_r,R_r,X,wt,si,d,p);

	  
	} // end update p matrix -> OK

      PExLen pexEp = getPexLen(exLen_r,p,tx_r,X);
      const mat pexLen = pexEp.pexLen;
      p = pexEp.p;

      mat oldLen = iLen;
      iLen = pexLen;
      mat secureP = pexLen;

      uvec findzero = find(secureP==0);
      vec nozero = secureP.elem(findzero);
      nozero.fill(min(iLen.elem(find(iLen > 0))));
      secureP.elem(findzero) = nozero;

      findzero = find(secureP < 0.001);
      nozero = secureP.elem(findzero);
      nozero.fill(0.001);
      secureP.elem(findzero) = nozero;
      
      checkError1 = as_scalar(max(max(abs(iLen - oldLen)/secureP)));

      if(checkError1 < err)
      	break;

      mat oldTheta = iTheta;

      if(X.n_cols == 1) // only one exon, and hence 1 transcript!      
	iTheta = Y;
      else
	{
	  mat iX = trans(X % iLen);

	  for(unsigned int idsmp=0; idsmp < len.size(); idsmp++)
	    {
	      mat iiX = iX * len(idsmp) / 1e+09;
	      vec iY = Y.col(idsmp);
	      iTheta.col(idsmp) = cpois_fitter(iiX, iY, maxit_fit, err, 0, 0, Q1, se_out, res_out);
	    }
	}

      
      mat secureTheta = iTheta;

      uvec findzeroT = find(secureTheta < 0.001);
      vec nozeroT = secureTheta.elem(findzeroT);
      nozeroT.fill(0.001);
      secureTheta.elem(findzeroT) = nozeroT;

      checkError2 = as_scalar(max(max(abs(iTheta - oldTheta)/secureTheta)));
      
      if(checkError2 < err)
      	break;

    }

  int kiter = 1;

  if(iter == maxit)
    kiter = 0;

  NumericMatrix iTheta_out = wrap(iTheta);
  iTheta_out.attr("dimnames") = dimnames;
  iTheta_out.attr("iter") = wrap(iter+kiter); // to make it comparable to R version

  // logfile.close();
  // remove(ogene);
  return iTheta_out;
  
  END_RCPP
}
