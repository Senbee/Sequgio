#include <Rcpp.h>
#include <tr1/unordered_map>
#include <map>

// #include <string>
// #include <iostream>
// #include <sstream>

using namespace Rcpp;
using namespace std;

// string Append(int _1, int _2)
// {
//   stringstream converter;
  
//   converter << _1 << "." << _2;
  
//   return converter.str();
// }


// Create a vector with exon.exon names and counts associated. Input: character matrix with 2 columns: left exon_name, right exon_name

RcppExport SEXP countEx(SEXP mat_reads, SEXP regid_)
{

  // IntegerMatrix mat_rd(mat_reads);
  CharacterMatrix mat_rd(mat_reads);
  CharacterVector regid(regid_);
  
  tr1::unordered_map<string,int> cvec; // hash to hold counts

  string iKey;

  for(int i=0; i != mat_rd.nrow(); i++)
    {
      // iKey =  Append(mat_rd(i,0),mat_rd(i,1));
      string str1 = as<string>(mat_rd(i,0));
      string str2 = as<string>(mat_rd(i,1));
      string rid = as<string>(regid(i));

      iKey =  str1 + "." + str2 + "__" + rid;
      
      cvec[iKey] += 1;

    }

  return wrap(cvec);
}


// This function assumes that txdb is already order according to "rank".

RcppExport SEXP makeExNames(SEXP txdb_, SEXP reg_id_, SEXP totLen_)
{

  BEGIN_RCPP

  List txdb = List(txdb_);
  CharacterVector reg_id = CharacterVector(reg_id_);

  int totLen = as<int>(totLen_);


  size_t k=0;
  vector<string> out_exnames(totLen);
  
  for(size_t itx=0; itx < txdb.size(); ++itx)
    {
      CharacterVector idf = txdb(itx);
      
      string reg = as<string>(reg_id(itx));
      vector<string> exnames = as< vector<string> >(idf);
      
      for(vector<string>::iterator iex1=exnames.begin(); iex1 < exnames.end(); ++iex1)
	for(vector<string>::iterator iex2=iex1; iex2 < exnames.end(); ++iex2) 
	  // it is possible to have reads within the same exon!!!
	  {
	    out_exnames[k] = *iex1+"."+*iex2+"__"+reg;
	    k++;
	  }
    }
  
  return(wrap(out_exnames));

  END_RCPP  
}
