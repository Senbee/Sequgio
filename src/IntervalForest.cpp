#include "IntervalForest.h"

IntervalForest _makeIntervalForest(SEXP r_listData)
{

  ////////////////////////////////////////////////////////////////////////////////

  Environment GA("package:GenomicAlignments");
  CharacterVector CIGAR_OPS = GA["CIGAR_OPS"];
  
  SEXP flag = R_NilValue;
  SEXP f = R_NilValue;
  LogicalVector with_ops(1), drop_empty_ranges(1), reduce_ranges(1);
  with_ops[0] = 0;
  drop_empty_ranges[0] = 1;
  reduce_ranges[0] = 0;
  IntegerVector space(1);
  space[0] = 3;

  
  ////////////////////////////////////////////////////////////////////////////////
  
  IntervalForest forest;

  S4 listData(r_listData);
  
  S4 partitioning = listData.slot("partitioning");
  DataFrame unlistData = listData.slot("unlistData");
  
  vector<string> seqNames = as< vector<string> >(partitioning.slot("NAMES"));
  
  // partitioning@end = stores the end of the 'partition'
  vector<int> partEnd = as< vector<int> >( partitioning.slot("end") );
  vIntIter partEnd_itr = partEnd.begin();
  partEnd.insert(partEnd_itr,0); // partEnd will store start at position i and end (not included) in position i+1 
  
  vector<string> strand = as< vector<string> >(unlistData["strand"]);

  
  IntegerVector pos = unlistData["start"];

  SEXP cigar = unlistData["cigar"];
  
  SEXP cigar_rng = cigar_ranges(cigar, flag, space, pos, f, CIGAR_OPS, drop_empty_ranges, 
				reduce_ranges,with_ops);
  
  
  CompressedIRangesList_holder cigar_hold = hold_CompressedIRangesList(cigar_rng);
  
  IntegerVector width = unlistData["width"];
  
  vector<string> exname = as< vector<string>  >(unlistData["exon_name"]);
  vector<string> txname = as< vector<string>  >(unlistData["tx_name"]);

  // loop along chromosomes

  for(pair<vStrIter, vIntIter> itr(seqNames.begin(),partEnd.begin()); itr.first != seqNames.end(); 
      ++itr.first, ++itr.second)
    {

      string chrName = *itr.first;
      int iStart = *itr.second;
      int iEnd = *(itr.second+1);
      
      // create an intervalVector
      intervalVector intervals;
      
      for(int j=iStart; j != iEnd; ++j)
	{

	  IRanges_holder elt = get_elt_from_CompressedIRangesList_holder(&cigar_hold,j);
	  int len_elt = get_length_from_IRanges_holder(&elt);
	  
	  vector<string> j_vec{txname[j], exname[j], strand[j]};

	  for(int k = 0; k != len_elt; ++k)
	    {
	      int j_start = get_start_elt_from_IRanges_holder(&elt,0);
	      int j_width = get_width_elt_from_IRanges_holder(&elt,0);
	      int j_stop = j_start + j_width - 1;
	      
	      intervals.push_back(interval( j_start, j_stop, j_vec ));

	    }

	}
      
      intervalTree tree = intervalTree(intervals);
      forest[chrName]  = tree;
    }


  return forest;

}




// [[Rcpp::export]]
SEXP makeForest(SEXP r_listData)
{

  XPtr< IntervalForest > forest_ptr(new IntervalForest,true);
  *forest_ptr =  _makeIntervalForest(r_listData);

  return forest_ptr;

}



/*
  makeTree

  INPUT: unlistData slot from a GRangeList
  OUTPUT: pointer to a IntervalTree object

 */


// [[Rcpp::export]]
// SEXP makeTree(SEXP r_unlistData)
// {
//   //  BEGIN_RCPP
//   intervalVector intervals;

//   // We need to partition within seqnames!!!
//   // for the test we just use chr4
//   S4 unlistData(r_unlistData);

//   // S4 seqRle = unlistData.slot("seqnames");
//   // IntegerVector seqVals = seqRle.slot("values");
//   // vector< string > seqLev = as< vector< string > >(seqVals.attr("levels"));
  
 
//   S4 strandRle = unlistData.slot("strand");
//   IntegerVector strandVal = strandRle.slot("values");
//   IntegerVector strandLen = strandRle.slot("lengths");

//   vector< int > strandValues;

//   for(pair<intIter, intIter> itr(strandLen.begin(),strandVal.begin()); itr.first != strandLen.end(); 
//       ++itr.first, ++itr.second)
//     for(int i = 0; i != *itr.first; ++i )
//       {
//       strandValues.push_back(*itr.second);
//     }

//   vector< string > strandLevs = as< vector< string > >(strandVal.attr("levels"));

//   DataFrame metadata = unlistData.slot("elementMetadata");
//   S4 ranges = unlistData.slot("ranges");
  
//   IntegerVector start = ranges.slot("start");
//   IntegerVector width = ranges.slot("width");
//   vector<string> exname = as< vector<string>  >(metadata["exon_name"]);
//   vector<string> txname = as< vector<string>  >(metadata["tx_name"]);

//   for(size_t i; i != start.size();++i)
//     {
//       int stop = start(i) + width(i) - 1;
//       int id = strandValues[i] - 1;
//       vector<string> vec{txname[i], exname[i], strandLevs[id]};
//       intervals.push_back(interval(start(i),stop,vec));
//     }


//   XPtr< intervalTree > tree_ptr(new intervalTree,true);
//   *tree_ptr =  intervalTree(intervals);

//   return(tree_ptr);
//   // END_RCPP
// }



// Templated function to map match results from right & left mates.
// Key = tx_name
// Values = ex_names

template<typename KeyType, typename LeftValue, typename RightValue>
map<KeyType, pair<LeftValue, RightValue> > IntersectMaps(const map<KeyType, LeftValue> & left, 
							 const map<KeyType, RightValue> & right)
{
  map<KeyType, pair<LeftValue, RightValue> > result;
  typename map<KeyType, LeftValue>::const_iterator il = left.begin();
  typename map<KeyType, RightValue>::const_iterator ir = right.begin();
  while (il != left.end() && ir != right.end())
    {
      if (il->first < ir->first)
	++il;
      else if (ir->first < il->first)
	++ir;
      else
        {
	  result.insert(make_pair(il->first, make_pair(il->second, ir->second)));
	  ++il;
	  ++ir;
        }
    }
  return result;
}


// find overlaps for gapped reads
MapVector matchGapped(IRanges_holder* irange, string strand, intervalTree& tree, 
			       int lenElt)
{

  int id = 0;
  MapVector results;
  int i_start, i_width, i_end;
  rangeVector queries;
  
  // create a vector of intervals. They all have to map to the same 'junction' exon
  // junction exons are coded as 2 or more intervals, same ex_name
  for(int id=0; id != lenElt; ++id)
    {
      i_start = get_start_elt_from_IRanges_holder(irange,id);
      i_width = get_width_elt_from_IRanges_holder(irange,id);
      i_end = i_start + i_width - 1;
      
      
      queries.push_back(iRange(i_start,i_end,strand));
    }

  // We expect ALL element of a gapped range to overlap EXACTLY to the same tx-exon pairs
  // We check this iteratively setting the first (the left-most element on + strand) overlap as pivot:
  // all others must be equal to that

  for (rangeVector::iterator q = queries.begin(); q != queries.end(); ++q) {
    MapVector cache_search;
    if(q == queries.begin())
      {
	tree.findEnd(q->start, q->stop, cache_search);

	if(cache_search.size() == 0) // no match...stop!
	  return(results);
	else
	  results = cache_search;
      }
    else if(q < (queries.end()-1)) // all the chunks but first and last
      {
	tree.findWithin(q->start, q->stop, cache_search);
	
	if(cache_search.size() == 0 || cache_search != results)
	  return(MapVector()); // no match here
      }
    else {

      tree.findStart(q->start, q->stop, cache_search);

      if(cache_search.size() == 0 || cache_search != results)
	return(MapVector()); // no match here
    }
  }

  // all 'chunks' must map to the same exon (in Sequgio's txdb junction regions are coded as a single region)
  return(results);
  
}

// find overlaps for ungapped reads

MapVector matchSimple(IRanges_holder* irange, string strand, intervalTree& tree)
{
  int id = 0;
  int i_start = get_start_elt_from_IRanges_holder(irange,id);
  int i_width = get_width_elt_from_IRanges_holder(irange,id);
  int i_end = i_start + i_width - 1;
  
  rangeVector queries;
  queries.push_back(iRange(i_start,i_end,strand));

  MapVector results;

  for (rangeVector::iterator q = queries.begin(); q != queries.end(); ++q) {
    tree.findWithin(q->start, q->stop, results);
  }

  return(results);

}


// [[Rcpp::export]]
void getOverlaps(SEXP r_forest, SEXP r_reads, SEXP r_exnames, SEXP r_colID, 
		 XPtr<BigMatrix> pBigMat)
{

  // Create the matrix accessor so we can get at the elements of the matrix.
  MatrixAccessor<int> ma(*pBigMat);

  const unsigned int colID = as<unsigned int>(r_colID);


  IntegerVector exrow(r_exnames); // NAMED integer vector with count matrix row # and ex names as rownames attribute
  vector<string> exname = as< vector<string> >(exrow.attr("names")); 

  // map count matrix rows to exon names: will provide an index to identify the row to add counts
  tr1::unordered_map<string,int> exmap;

  for(size_t i=0;  i < exrow.size(); i++)
    {
      exmap[exname[i]] = exrow[i]; 
    }



  S4 forest_obj(r_forest);
  SEXP forest_ptr;
  forest_ptr = forest_obj.slot("ptr");

  XPtr< IntervalForest > forest_p(forest_ptr);
  IntervalForest&  forest = *forest_p;

  List reads(r_reads);






  IntegerVector groupid = reads["groupid"];
  IntegerVector strand_v = reads["strand"]; //1=+; 2=-; *=3
  CharacterVector cigar = reads["cigar"];
  int nreads = groupid.size();



  // create a string vector of strand indicator
  const vector<string> strandLev{"+","-","*"}; 
  vector < string > strand;
  for(intIter itr = strand_v.begin(); itr != strand_v.end(); ++itr)
    {
      string tmpVal = strandLev[(*itr)-1];
      strand.push_back(tmpVal);
    }


  // create a string vector of seqnames (chromosomes)
  IntegerVector seqname_v = reads["rname"];
  vector<string> seq_lev = as< vector<string> >(seqname_v.attr("levels"));
  vector<string> seqname;
  for(intIter itr = seqname_v.begin(); itr != seqname_v.end(); ++itr)
    {
      string tmpVal = seq_lev[(*itr)-1];
      seqname.push_back(tmpVal);
    }
  

  
  Environment GA("package:GenomicAlignments");
  CharacterVector CIGAR_OPS = GA["CIGAR_OPS"];

  SEXP flag = R_NilValue;
  SEXP f = R_NilValue;
  LogicalVector with_ops(1), drop_empty_ranges(1), reduce_ranges(1);
  with_ops[0] = 0;
  drop_empty_ranges[0] = 1;
  reduce_ranges[0] = 0;

  IntegerVector pos, space(1);
  pos = reads["pos"];

  space[0] = 3;


  // "CompressedIRangesList": every element is a IRanges object
  SEXP cigar_rng = cigar_ranges(cigar, flag, space, pos, f, CIGAR_OPS, drop_empty_ranges, 
  				reduce_ranges,with_ops);

  CompressedIRangesList_holder cigar_hold = hold_CompressedIRangesList(cigar_rng);


  vector<int> lenOut;

  vector< string > collect_results;

  string i_seqname = "";
  intervalTree tree;
  
  for(size_t i=0; i != nreads; i+=2) // increment by 2 as reads are paired
    {

      string i_strand = strand[i];
      
      if(i_seqname != seqname[i])
	{
	  i_seqname = seqname[i];
	  IntervalForest::const_iterator map_itr = forest.find(i_seqname);
	  if(map_itr == forest.end())
	    {
	      throw std::range_error("Seqname not found in annotation");
	    }
	  tree = map_itr->second;
	}

      MapVector result_plus, result_neg;
      int len_elt;
      IRanges_holder elt;
      
      // exon names pairs in counting are listed in direction + -> -, so we always start from "+" mate
      int j = (i_strand == "+") ? i : (i+1); // recall: condition ? value_if_true : value_if_false

      // This should go in a function that process the i-th read

      // place 1st mate ("+")
      elt = get_elt_from_CompressedIRangesList_holder(&cigar_hold,j);
      len_elt = get_length_from_IRanges_holder(&elt);
      
      if(len_elt == 1)
	{
	  result_plus = matchSimple(&elt, i_strand, tree);
	  // some reads might fall in exons that are shared by different tx, so appear more than once!
	  // if(result_plus.size() > 1)
	  //   {
	  //     copy(result_plus.begin(), result_plus.end(), ostream_iterator<string>(cout, "-"));
	  //     cout << j << endl;
	  //   }
	}
      else
	{
	  result_plus = matchGapped(&elt, i_strand,tree,len_elt);
	}
 
      
      // if no overlap for + read
      if(result_plus.size() == 0)
      	continue;

      // get next read, - strand. Is it listed before or after + strand? -> k
      int k = (j > i) ? i : i+1;

      elt = get_elt_from_CompressedIRangesList_holder(&cigar_hold,k);
      len_elt = get_length_from_IRanges_holder(&elt);
      
      if(len_elt == 1)
	{
	  result_neg = matchSimple(&elt, strand[k],tree);
	}
      else
	{
	  result_neg = matchGapped(&elt, strand[k],tree,len_elt);
	}

      
      // no match from - strand read
      if(result_neg.size() == 0)
      	continue;

      // cout << i << " | " << result_plus.size()<< "/" << result_neg.size();

      MapPairVector merge_results = IntersectMaps(result_plus,result_neg);
      // we might get no intersection!!!
      // cout << "=> " << merge_results.size() << endl;
      
      if(merge_results.size() > 0)
	{
	  for(MapPairVector::iterator iter=merge_results.begin();
	      iter != merge_results.end(); ++iter) {

	    string iKey = iter->second.first + "." + iter->second.second;
	    // collect_results.push_back(iKey);
	    
	    if(exmap.count(iKey) > 0)
	      {
		int rowID = exmap[iKey];
		ma[colID][rowID] =  1 + ma[colID][rowID];
	      }
	    
	  }
	  
	}

    }

}

