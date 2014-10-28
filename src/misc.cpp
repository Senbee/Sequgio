// Time-stamp: <02-03-2013 17:00:55 on Goliath.med.unibs.it>

#include <RcppArmadillo.h>
#include <cmath>
#include <cstdlib>
#include "ranker.h"
#include "headers.h"

// #include <ctime>

 // If we use GSL: approx & smooth.spline
 // #include <gsl/gsl_interp.h>



 ////////////////////////////////////////////////////////////////////////////////


mat Zi2V(const mat& Zi, const vec& wt, const uvec& r)
 {

   mat newt(wt);
   newt.reshape(Zi.n_rows,Zi.n_cols); // from n sample*n exons vector to n sample-x-n exons matrix

   mat tZi = Zi % Zi % newt;
   rowvec asum = sum(tZi,0);

   mat ans = diagmat(asum.elem(r));

   return ans;
 }


vec ZiVy(const mat& Zi, const vec& y, const vec& wt, const uvec r)
{

  mat newt = wt % y;

  newt.reshape(Zi.n_rows,Zi.n_cols); // from n_sample*n_exons vector to n_sample-by-n_exons matrix

  mat tZi = Zi % newt;

  rowvec asum = sum(tZi,0);

  vec ans = conv_to< vec >::from(asum.elem(r));



  return ans;
}


vec Zi_b(const mat& Xf, const mat& B, const  std::vector<std::string> index, 
	 const std::vector<std::string> tx_names, int ncolY)
{
  // we assume order in B rows = Xf cols (they are tx names)
   mat nXf;
   mat nB;

   int k = 0;

   for(unsigned int i=0; i < Xf.n_cols; i++)
     {
       bool keep = true;

       for(unsigned int j = 0; j < index.size(); j++)
	 {
	   if(tx_names[i] == index[j])
	     {
	       keep = false;
	       break;
	     }
	 }

       if(keep)
	 {
	   nXf.insert_cols(k,Xf.col(i));
	   nB.insert_rows(k,B.row(i));
	   k++;
	 }
     }


   mat tnB = trans(nB);
   mat tmpB;

   for(unsigned int j =0; j < tnB.n_rows; j++)
     for(int i = 0; i < ncolY; i++)
       tmpB.insert_rows(i+j*ncolY,tnB.row(j));

   mat abb = nXf % tmpB;

   vec ans = sum(abb,1);

   return ans;
 }


vec Zb(const mat& Xf, const mat& B, int ncolY)
 {
   mat tB = trans(B);
   mat tmpB;

   for(unsigned int j =0; j < tB.n_rows; j++)
     for(int i = 0; i < ncolY; i++)
       tmpB.insert_rows(i+j*ncolY,tB.row(j));

   mat abb = Xf % tmpB;

   vec ans = sum(abb,1);

   return ans;
 }


List Yfun(mat Y, const mat& B, const double Q, const mat& Xf)
{

  static double (*absolute)(double) = &std::abs;
  
  int ncolY = Y.n_cols;
  
  vec eta = Zb(Xf,B,ncolY);
  
  vec mu(eta);
  Y = trans(Y);
  Y.reshape(Y.n_cols*Y.n_rows,1);
  
  vec res = Y - mu;
  
  for(unsigned int i=0; i < mu.n_elem; i++)
    {
      if(mu(i) < 0.1) {mu(i) = 0.1;}
    }
  
  vec wt = 1/mu;
  vec ares(res.n_elem);
  
  std::transform(res.begin(),res.end(),ares.begin(), absolute);

  double q3 = quantile(conv_to< std::vector<double> >::from(ares), Q);
  vec nwt = (q3/(ares + 0.01))/mu;
  
  for(unsigned int i=0; i < wt.n_elem; i++)
    {
      if(ares(i) >= q3) {wt(i) = nwt(i);}
      // if(wt(i) < 0.1) {wt(i) = 0.1;}
      if(wt(i) > 10) {wt(i) = 10;}
    }
  
  
  vec y_out = eta + res;
  
  return List::create(Named("wt")=wt,
			    Named("Y")=y_out);
}


////////////////////////////////////////////////////////////////////////////////
// Tested on 28 Jul 2011 => OK

// NumericMatrix fun2(const StringVector x, int d)
mat fun2(const StringVector x, int d)
{

  int len = x.size()-1;
  
  mat d1 = eye<mat>(len,len);

  mat izeri = zeros<mat>(len,1);
  mat d2 = join_rows(d1,izeri);
  mat d3 = join_rows(izeri,d1);
  d2 = -1*d2;
  mat delta = d2 + d3;

  mat out;

  if(d == 1)
    out = trans(delta) * delta;
  else
    {
      mat delta2 = delta.submat(span(1,len-1),span(1,len)); // remove first row & col
      mat delta3 = delta2 * delta;
      out = trans(delta3) * delta3;
    }

  return out;

 }

List getR(const List tx, int deriv)
{
  // tx is a named list with transcript exon names
  List output(tx.size()); // output list
  
  // second input in transform must be of the same length as the first. So we need to create a vector
  // of the same length as input holding a repetition on deriv value
  vec d(tx.size());
  d.fill(deriv);

  std::transform(tx.begin(), tx.end(), d.begin(),output.begin(), fun2);

  output.names() = tx.names();
  
  return output;
}



////////////////////////////////////////////////////////////////////////////////


void fitCr(const mat& X, const mat& Y,  const NumericVector exLen_r, const mat& iTheta, const vec& Len, const int& maxit, 
	   double error_limit, int std_err, const int resd, const double& Q1, const List tx,  const unsigned int iter, 
	   mat& p, const unsigned int df)
{
  static double mscale = 1e+09;

  // mat tmp_betas(X.n_rows,X.n_cols);
  // mat out_betas(X.n_rows,X.n_cols);
  vec beta(p.n_rows);

  rowvec rLen = conv_to<rowvec>::from(Len);

  mat ixMat(X.n_rows,iTheta.n_cols);
  mat iLen(iTheta.n_rows,rLen.n_elem);

  mat aMat;
  aMat.copy_size(ixMat);

  mat taMat(aMat.n_cols,aMat.n_rows);

  vec ixCol(X.n_rows);
  rowvec rY(Y.n_cols);
  vec iY(Y.n_cols);

  for(unsigned int i = 0; i < X.n_cols; i++)
    {
      ixCol = X.col(i);

      rY = Y.row(i);
      iY = conv_to<vec>::from(rY);

      for(unsigned int j=0; j < iTheta.n_cols; j++)
	ixMat.col(j) = ixCol;
      
      for(unsigned int k=0; k < iTheta.n_rows; k++)
	iLen.row(k) = rLen;

      double konst = exLen_r[i] / mscale;
      aMat = ixMat % iTheta % iLen * konst;

      // Not used yet /////
      vec se_out;
      vec res_out;
      
      if(resd == 0)
      	res_out.fill(0);
      
      if(std_err == 0)
      	se_out.fill(0);
      ////////////////////
      
      taMat = trans(aMat);

      beta = cpois_fitter(taMat, iY, maxit, error_limit, std_err, resd, Q1, se_out, res_out);

      p.col(i) = beta;

    }

  getStart(tx, exLen_r,  iter, p, df);

}


////////////////////////////////////////////////////////////////////////////////
// p_in(matrix)  and tx_in (list) must have the same names: names(tx_in) == rownames(p_in)
// and tx elements must have the same order as in p_in

void getStart(List tx, const NumericVector exLen,  const unsigned int iter, mat& p, const unsigned int df)
{
  
  vec ifit;

  std::vector<std::string> nomi = tx.names();

  static std::map<std::string,int> exLen_map;

  std::vector<std::string> ex_nomi = exLen.names();

  // create map holding names exon position. p & exLen must have the same ordered names
  for(unsigned int i=0; i < ex_nomi.size(); i++)
    exLen_map[ex_nomi[i]] = i;

  for(int i=0; i < tx.size(); i++)
    {

      std::vector<int> cLen;
      std::vector<double> cP;
      std::vector<std::string> itx = tx[i];

      rowvec hold_p = p.row((unsigned int) i);
      std::vector<unsigned int> hold_pos;


      for(std::vector<std::string>::iterator j=itx.begin(); j!=itx.end();++j)
	{
	  unsigned int pos = exLen_map.find(*j)->second;
	  cLen.push_back(exLen[pos]);
	  cP.push_back(p(i,pos));
	  hold_pos.push_back(pos);
	}


      vec pSum(cLen.size());
      std::partial_sum(cLen.begin(), cLen.end(),pSum.begin());
      
      uvec vhold_pos = conv_to< uvec >::from(hold_pos);

      vec x_in(pSum);
      vec y_in  = conv_to< colvec >::from(cP);

      if(cLen.size() > 3)
      	{
	  if(iter == 0) // first iteration
	    {
	      ifit = smspline(x_in,y_in,3);
	    }
	  else
	    {
	      ifit = smspline(x_in,y_in,df);
	    }
	  hold_p.elem(vhold_pos) = conv_to< rowvec >::from(ifit);
	}
      p.row((unsigned int) i) = hold_p;
    }

}




////////////////////////////////////////////////////////////////////////////////

mat getXf(const mat& X, const mat& iTheta, const vec& len_in, const vec& exLen)
{

  rowvec len = conv_to<rowvec>::from(len_in);

  mat LenMat;
  for(unsigned int j=0; j < iTheta.n_rows; j++)
    LenMat.insert_rows(j,len);

  // It's probably better to set output dimension and then insert values. Let's try sometime...
  mat output;
  
  for(unsigned int i=0; i < X.n_cols; i++)
    {
      mat tmpX;

      for(unsigned int j=0; j < iTheta.n_cols; j++)
	tmpX.insert_cols(j,X.col(i));

      mat cMat = trans((tmpX % iTheta % LenMat * exLen(i)) / 1e+09);
      
      output = join_cols(output,cMat);
    }

  return output;

}



////////////////////////////////////////////////////////////////////////////////

void getZi(cube& Zi_in, const mat& Xf)
{


  unsigned int i,j,st,en,nr,nc;
  
  mat inZi = zeros<mat>(Zi_in.n_rows,Zi_in.n_cols);


  nr = Zi_in.n_rows;
  nc = Zi_in.n_cols;


  for(i = 0; i < Xf.n_cols; i++)
    {
      st = 0;
      en = nr - 1;

      for(j = 0; j < nc; j++)
	{
	  inZi.col(j) = Xf(span(st,en),i);
	  st += nr;
	  en += nr;
	}
      Zi_in.slice(i) = inZi;
    }
}



////////////////////////////////////////////////////////////////////////////////

// The getLambda function loops along lambda values. But within loops also over tx. Must return a matrix
// n of tx -times- n of lambda values (5)

mat getLambda(const vec& lambda, const cube& Zi_r, const mat& Xf, const mat& B, const mat& Y, const double Q, 
	      const mat& iTheta, List tx_r, mat X, List R_r, 
	      const std::vector< std::string >& tx_names)
{


  int ncolY =  Y.n_cols;

  ////////////////////////////////////

  mat dfall = zeros<mat>(Xf.n_cols,lambda.n_elem);
  mat iB(B);
  mat b(B);
  vec yc = zeros<vec>(Y.n_cols * Y.n_rows);
  vec wt = ones<vec>(Y.n_cols * Y.n_rows);

  mat Zi(Zi_r.n_rows,Zi_r.n_cols);

  // see .getTestBeta in SC code
  for(unsigned int la = 0; la < lambda.n_elem; la++) // loops along lambda values
    {

      iB = B;

      for(int it=0; it < 5; it++) // updates p matrix 5 times
	{
	  List fitY = Yfun(Y,iB,Q,Xf);
	  yc = as<vec>(fitY["Y"]);
	  wt = as<vec>(fitY["wt"]);

	  b = iB;

	  for(unsigned int idx = 0; idx < Xf.n_cols; idx ++) // loops along tx
	    {
	      
	      // this must be checked. Why use a vector??
	      std::vector<std::string> index;
	      index.push_back(tx_names[idx]);

	      Zi = Zi_r.slice(idx);

	      StringVector itx = tx_r[idx];

	      bool allzero = sum(iTheta.row(idx)) == 0;

	      if(!(itx.size() == 1 || allzero))
		{

		  vec iEx = conv_to< vec>::from(X.row(idx));
		  uvec sel = find(iEx != 0);
		  mat iR = R_r[idx];

		  vec iP = conv_to< vec>::from(b.row(idx));

		  if(Xf.n_cols > 1)
		    {
		      yc = as<vec>(fitY["Y"]) - Zi_b(Xf,b,index,tx_names,ncolY);
		    }

		  vec iZiVy = ZiVy(Zi, yc, wt, sel);

		  mat Lambda(iR.n_rows,iR.n_cols);
		  Lambda.fill(lambda(la));

		  mat iZi2 = Zi2V(Zi, wt, sel) + (iR % Lambda);
		  vec iFit = fitter(iZi2,iZiVy);
		  // vec iFit = solve(iZi2,iZiVy);
		  // uvec isneg = find(iFit < 0);
		  // if(sum(isneg))
		  //   iFit.elem(isneg) = zeros<vec>(isneg.n_elem);

		  iP.elem(sel) = iFit;
		  iB.row(idx) = conv_to<rowvec>::from(iP);
		}
	    } // end iteration along tx (idx)
	} // end update p matrix


      for(unsigned int idx = 0; idx < Xf.n_cols; idx ++) // loops along tx
	{	

	  bool allzero = sum(iTheta.row(idx)) == 0;
	  
	  if(!allzero)
	    {
	      
	      Zi = Zi_r.slice(idx);
	      vec iEx = conv_to< vec>::from(X.row(idx));
	      uvec sel = find(iEx != 0);
	      // mat iR = as<mat>(R_r[tx_names[idx]]);
	      mat iR = R_r[idx];

	      mat Lambda(iR.n_rows,iR.n_cols);
	      Lambda.fill(lambda(la));
	      mat iZi2 = Zi2V(Zi, wt,sel);// + (iR % Lambda);

	      mat iMat = iZi2 + (iR % Lambda);
	      
	      uvec iszero = find(iMat == 0); // to avoid computational singularity
	      
	      if(sum(iszero))
		{
		  vec repzero(iszero.n_elem);
		  repzero.fill(0.001);
		  iMat.elem(iszero) = repzero;
		}
	      
	      mat oMat = inv(iMat);
	      mat aMat = oMat % iZi2;
	      dfall(idx,la) = sum(sum(aMat));
	    }
	}
    } // end iteration along lambda (la)

  return dfall;

}

////////////////////////////////////////////////////////////////////////////////

void getsi(const cube& Zi_r, List R_r, const mat& X, const vec& wt, vec& si, int d, const mat& p)
{
  

  std::vector< std::string > cnames = R_r.attr("names");

  mat Zi(Zi_r.n_rows,Zi_r.n_cols);

  for(unsigned int idx = 0; idx < X.n_rows; idx ++) // loops along tx
    {	

      vec iEx = conv_to< vec>::from(X.row(idx));
      uvec sel = find(iEx != 0);
      Zi = Zi_r.slice(idx);
      // mat iR = as<mat>(R_r[cnames[idx]]);
      mat iR = R_r[idx];

      vec ib = conv_to<vec>::from(p.row(idx));
      ib = ib.elem(sel);
      mat iZi2 = Zi2V(Zi,wt,sel);
      mat iZi2R = iZi2 + (iR/si[idx]);
      mat dMat = eye<mat>(iR.n_rows,iR.n_cols);
      mat sMat(iR.n_rows,iR.n_cols);
      for(unsigned int cId = 0; cId < dMat.n_cols; cId++)
	{
	  sMat.col(cId) = fitter(iZi2R,dMat.col(cId));
	}
      mat aMat = sMat % iR;
      double aSum = sum(sum(aMat));
      mat rMat(iR.n_rows,iR.n_cols);
      rMat.fill(aSum);

      mat bMat = iR + rMat;
      mat cMat = conv_to<rowvec>::from(ib) * bMat * ib;
      double isi =  as_scalar(cMat);

      si(idx) = isi / (double)(ib.n_elem - d) ;
    }
  
}



////////////////////////////////////////////////////////////////////////////////

PExLen getPexLen(NumericVector exLen_r, const mat& p, List tx_r, const mat& X)
{

  vec exLen = as< vec >(exLen_r);
  rowvec rexLen = conv_to< rowvec >::from(exLen);

  mat exLen_mat = zeros<mat>(p.n_rows,p.n_cols);
  
  for(unsigned int r=0; r < p.n_rows;r++)
    {
      exLen_mat.row(r) = rexLen;
    }

  mat pexLen_tmp = p % exLen_mat;

  // double min_ex = min(exLen);
  double min_ex = 1.0;

  mat ans(p.n_rows,p.n_cols);
  uvec sel, nosel;
  vec iEx(X.n_cols);
  
  double txLen=0;
  rowvec pex(pexLen_tmp.n_cols);

  for(int idx = 0; idx < tx_r.size(); idx ++) // loops along tx
    {
      iEx = conv_to< vec >::from(X.row(idx));
      sel = find(iEx != 0);
      nosel = find(iEx == 0);

      std::vector< std::string > exName = tx_r[idx];

      txLen = 0;
      for(unsigned int tx_i=0; tx_i < exName.size(); tx_i++)
	{
	  txLen += exLen_r[exName[tx_i]];
	}

      pex = pexLen_tmp.row(idx);
      uvec wmin = find( pex < min_ex);
      rowvec fmin(wmin.n_elem);
      fmin.fill(min_ex);
      pex.elem(wmin) = fmin;

      vec issel = pex.elem(sel);
      uvec isone = find(issel == 1);

      if(isone.n_elem == pex.n_elem)
	{
	  ans.row(idx) = exLen;
	}
      else
	{
	  pex.elem(nosel) = exLen.elem(nosel);
	  pex.elem(sel) = pex.elem(sel)/sum(pex.elem(sel)) * txLen;
      
	  ans.row(idx) = pex;
	}
	
    }

  mat newp = ans/exLen_mat;
  
  PExLen out;
  
  out.pexLen = ans;
  out.p = newp;

  return out;
}


