/* 
   Count reads for every 'exon/region' and save them directly into
   shared big.matrix object
*/

#include <Rcpp.h>
#include <tr1/unordered_map>

// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>


using namespace Rcpp;
using namespace std;

// http://gallery.rcpp.org/articles/using-bigmemory-with-rcpp/

// [[Rcpp::export]]
void countExBM(SEXP mat_reads, SEXP regid_, SEXP exnames_, SEXP colID_, 
			  XPtr<BigMatrix> pBigMat)
{

  int rowID;

  // Create the matrix accessor so we can get at the elements of the matrix.
  MatrixAccessor<int> ma(*pBigMat);


  //NOTE:for a MatrixAccessor ma, the i-th row and j-th column is accessed with ma[j][i]

  CharacterMatrix mat_rd(mat_reads);
  CharacterVector regid(regid_);

  const unsigned int colID = as<unsigned int>(colID_);


  IntegerVector exrow(exnames_); // NAMED integer vector with count matrix row # and ex names as rownames attribute
  vector<string> exname = as< vector<string> >(exrow.attr("names")); 
  
  tr1::unordered_map<string,int> exmap;

  for(size_t i=0;  i < exrow.size(); i++)
    {
      exmap[exname[i]] = exrow[i]; 
    }


  // int totReads = 0;
  
  for(size_t i=0; i < mat_rd.nrow(); i++)
    {
      string str1 = as<string>(mat_rd(i,0));
      string str2 = as<string>(mat_rd(i,1));
      string rid = as<string>(regid(i));
      
      string iKey =  str1 + "." + str2 + "__" + rid;
      
      if(exmap.count(iKey) > 0)
	{
	  rowID = exmap[iKey];
	  ma[colID][rowID] =  1 + ma[colID][rowID];
	  // totReads += 1;
	}
      
    }

  //  return(wrap(0));
}




// RcppExport SEXP countEx(SEXP mat_reads, SEXP regid_, SEXP exnames_)
// {

//   CharacterMatrix mat_rd(mat_reads);
//   CharacterVector regid(regid_);
//   CharacterVector exnames(exnames_);

//   // Build a set to store 'known' exon names. Only those present there would be counted
//   vector<string> exnames_v = as< vector<string> >(exnames);
//   tr1::unordered_set<string> exset(exnames_v.begin(),exnames_v.end());
  
//   tr1::unordered_map<string,int> cvec; // hash to hold counts
  
//   string iKey;
  
//   for(size_t i=0; i != mat_rd.nrow(); i++)
//     {
//       string str1 = as<string>(mat_rd(i,0));
//       string str2 = as<string>(mat_rd(i,1));
//       string rid = as<string>(regid(i));

//       iKey =  str1 + "." + str2 + "__" + rid;
      
//       if(exset.count(iKey) > 0)
// 	cvec[iKey] += 1;

//     }

//   return wrap(cvec);
// }
