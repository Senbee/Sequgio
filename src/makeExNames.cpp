#include <Rcpp.h>
#include <boost/unordered_set.hpp>

using namespace Rcpp;
using namespace std;


// This function assumes that txdb is already order according to "rank".

// [[Rcpp::export]]
SEXP makeExNames(SEXP txdb_, SEXP reg_id_)
{

  BEGIN_RCPP

  List txdb = List(txdb_);
  CharacterVector reg_id = CharacterVector(reg_id_);

  boost::unordered_set<string> cvec; // set to hold names

  string iKey;

  for(size_t itx=0; itx < txdb.size(); ++itx)
    {
      CharacterVector idf = txdb(itx);
      
      string reg = as<string>(reg_id(itx));
      vector<string> exnames = as< vector<string> >(idf);
      
      for(vector<string>::iterator iex1=exnames.begin(); iex1 < exnames.end(); ++iex1)
	for(vector<string>::iterator iex2=iex1; iex2 < exnames.end(); ++iex2) 
	  // it is possible to have reads within the same exon!!!
	  {
	    iKey = *iex1+"."+*iex2+"__"+reg;
	    cvec.insert(iKey);
	  }
    }
  

  vector<string> out_exnames(cvec.size());
  
  copy(cvec.begin(),cvec.end(),out_exnames.begin());

  return(wrap(out_exnames));

  END_RCPP  
}

