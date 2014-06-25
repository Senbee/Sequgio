#ifndef __INTERVALFOREST_H
#define __INTERVALFOREST_H

// [[Rcpp::depends(IRanges,S4Vectors,BH,bigmemory)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <Rdefines.h>

#include <IRanges_interface.h>
#include <_IRanges_stubs.c>

#include <vector>
#include <algorithm>

// bigmatrix related
#include <tr1/unordered_map>
#include <bigmemory/MatrixAccessor.hpp>

using namespace Rcpp;
using namespace std;

typedef IntegerVector::iterator intIter;
typedef CharacterVector::iterator charIter;
typedef vector<int>::iterator vIntIter;
typedef vector<string>::iterator vStrIter;
typedef map< string, string> MapVector;
typedef map< string, pair<string,string> > MapPairVector;



template <class T, typename K = int>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};


class IntervalRange {
public:
  int start;
  int stop;
    IntervalRange(int s, int e)
      : start(s)
      , stop(e)
  { }
};



template <class T, typename K>
int intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
int intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
ostream& operator<<(ostream& out, Interval<T,K>& i) {
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = int>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};


// IntervalTree class

template <class T, typename K = int>
class IntervalTree {
  
public:
  typedef Interval<T,K> interval;
  typedef vector<interval> intervalVector;
  typedef IntervalTree<T,K> intervalTree;


  intervalVector intervals;
  intervalTree* left;
  intervalTree* right;
  int center;
  
  IntervalTree<T,K>(void)
  : left(NULL)
  , right(NULL)
  , center(0)
  { }
  
  IntervalTree<T,K>(const intervalTree& other) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = (intervalTree*) malloc(sizeof(intervalTree));
      *left = *other.left;
    } else {
      left = NULL;
    }
    if (other.right) {
      right = new intervalTree();
      *right = *other.right;
    } else {
      right = NULL;
    }
  }
  
  IntervalTree<T,K>& operator=(const intervalTree& other) {
    center = other.center;
    intervals = other.intervals;
    if (other.left) {
      left = new intervalTree();
      *left = *other.left;
    } else {
      left = NULL;
    }
    if (other.right) {
      right = new intervalTree();
      *right = *other.right;
    } else {
      right = NULL;
    }
    return *this;
  }
  
  IntervalTree<T,K>(
		    intervalVector& ivals,
		    unsigned int depth = 16,
		    unsigned int minbucket = 64,
		    int leftextent = 0,
		    int rightextent = 0,
		    unsigned int maxbucket = 512
		    )
  : left(NULL)
    , right(NULL)
  {
    
    --depth;
    IntervalStartSorter<T,K> intervalStartSorter;
    if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
      sort(ivals.begin(), ivals.end(), intervalStartSorter);
      intervals = ivals;
    } else {
      if (leftextent == 0 && rightextent == 0) {
	// sort intervals by start
	sort(ivals.begin(), ivals.end(), intervalStartSorter);
      }
      
      int leftp = 0;
      int rightp = 0;
      int centerp = 0;
      
      if (leftextent || rightextent) {
	leftp = leftextent;
	rightp = rightextent;
      } else {
	leftp = ivals.front().start;
	vector<K> stops;
	stops.resize(ivals.size());
	transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
	rightp = *max_element(stops.begin(), stops.end());
      }
      
      //centerp = ( leftp + rightp ) / 2;
      centerp = ivals.at(ivals.size() / 2).start;
      center = centerp;
      
      intervalVector lefts;
      intervalVector rights;
      
      for (typename intervalVector::iterator i = ivals.begin(); i != ivals.end(); ++i) {
	interval& interval = *i;
	if (interval.stop < center) {
	  lefts.push_back(interval);
	} else if (interval.start > center) {
	  rights.push_back(interval);
	} else {
	  intervals.push_back(interval);
	}
      }
      
      if (!lefts.empty()) {
	left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
      }
      if (!rights.empty()) {
	right = new intervalTree(rights, depth, minbucket, centerp, rightp);
      }
    }
  }
  
  void findWithin(K start, K stop, MapVector& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start >= interval.start && stop <= interval.stop) {
	  vector<string> insVal = interval.value;
	  // Insert tx_name & ex_name. tx_name will be used to match pair ends (must fall within same tx!)
	  contained.insert(MapVector::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findWithin(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findWithin(start, stop, contained);
    }
    
  }
  
  // find query to have same start as tree interval
  void findStart(K start, K stop, MapVector& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start == interval.start && stop <= interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(MapVector::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findStart(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findStart(start, stop, contained);
    }
  }
  
  
  
  // find query to have same end as tree interval
  void findEnd(K start, K stop, MapVector& contained) {
    if (!intervals.empty() && ! (stop < intervals.front().start)) {
      for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
	interval& interval = *i;
	if (start >= interval.start && stop == interval.stop) {
	  vector<string> insVal = interval.value;
	  contained.insert(MapVector::value_type(insVal[0],insVal[1]));
	}
      }
    }
    
    if (left && start <= center) {
      left->findEnd(start, stop, contained);
    }
    
    if (right && stop >= center) {
      right->findEnd(start, stop, contained);
    }
  }
  
  // find query to have exact overlap as tree interval
  void findExact(K start, K stop, MapVector& contained) {
  if (!intervals.empty() && ! (stop < intervals.front().start)) {
    for (typename intervalVector::iterator i = intervals.begin(); i != intervals.end(); ++i) {
      interval& interval = *i;
      if (start == interval.start && stop == interval.stop) {
	vector<string> insVal = interval.value;
	contained.insert(MapVector::value_type(insVal[0],insVal[1]));
      }
    }
  }
  
  if (left && start <= center) {
    left->findExact(start, stop, contained);
  }
  
  if (right && stop >= center) {
    right->findExact(start, stop, contained);
  }
  }
  
  ~IntervalTree(void) {
    // traverse the left and right
    // delete them all the way down
    if (left) {
      delete left;
    }
    if (right) {
      delete right;
    }
  }
  
};



// IntervalForest Class: should hold as many trees as chromosomes. Use a map to retrive them with 'string'

// template <class T, typename K = int>
// class IntervalForest {
  
// public:

//   typedef Interval< vector<string> > interval;
//   typedef vector<interval> intervalVector;
//   typedef IntervalTree< vector<string> > intervalTree;
//   typedef IntegerVector::iterator intIter;
//   typedef CharacterVector::iterator charIter;
//   typedef vector<int>::iterator vIntIter;
//   typedef vector<string>::iterator vStrIter;

//   map<string,intervalTree> forest;
  
//   // Constructors

//   IntervalForest(void) 
//   : forest()
//   {}

//   IntervalForest(SEXP& r_listData) {

//     S4 listData(r_listData);

//     S4 partitioning = listData.slot("partitioning");
//     S4 unlistData = listData.slot("unlistData");

//     vector<string> seqNames = as< vector<string> >(partitioning.slot("NAMES"));

//     // partitioning@end = stores the end of the 'partition'
//     vector<int> partEnd = as< vector<int> >( partitioning.slot("end") );
//     vIntIter partEnd_itr = partEnd.begin();
//     partEnd.insert(partEnd_itr,0); // partEnd will store start at position i and end (not included) in position i+1 

//     S4 strandRle = unlistData.slot("strand");
//     IntegerVector strandVal = strandRle.slot("values");
//     IntegerVector strandLen = strandRle.slot("lengths");
    
//     vector< int > strandValues;
//     int counter = 0;
//     for(pair<intIter, intIter> itr(strandLen.begin(),strandVal.begin()); itr.first != strandLen.end(); 
// 	++itr.first, ++itr.second)
//       {
// 	unsigned int uInt = (unsigned int)*itr.first;
// 	strandValues.insert(strandValues.end(),uInt,*itr.second);
//       }
    
//     vector< string > strandLevs = as< vector< string > >(strandVal.attr("levels"));
    
//     DataFrame metadata = unlistData.slot("elementMetadata");

//     S4 ranges = unlistData.slot("ranges");
    
//     IntegerVector start = ranges.slot("start");

//     IntegerVector width = ranges.slot("width");

//     vector<string> exname = as< vector<string>  >(metadata["exon_name"]);
//     vector<string> txname = as< vector<string>  >(metadata["tx_name"]);

//     // loop along chromosomes
//     cout << start.size()  << endl;
//     for(pair<vStrIter, vIntIter> itr(seqNames.begin(),partEnd.begin()); itr.first != seqNames.end(); 
// 	++itr.first, ++itr.second)
//       {
// 	string chrName = *itr.first;
// 	int iStart = *itr.second;
// 	int iEnd = *(itr.second+1);

// 	// create an intervalVector
// 	intervalVector intervals;

// 	for(int j=iStart; j != iEnd; ++j)
// 	  {
// 	    int stop = start(j) + width(j) - 1;
// 	    int id = strandValues[j] - 1;
// 	    vector<string> vec{txname[j], exname[j], strandLevs[id]};
// 	    intervals.push_back(interval(start(j),stop,vec));
// 	  }

// 	intervalTree tree = intervalTree(intervals);
// 	forest[chrName]  = tree;
//       }
//   }


//   // Assignment/copy operator

//   IntervalForest& operator=(const IntervalForest& other) {
//     forest = other.forest;
//     return *this;
//   }


//   // Destructor: clear the map
//   ~IntervalForest(void) {
//     forest.clear();
//   }
  
// };





///////////////////////////////////////////////////////////////

typedef Interval< vector<string> > interval;
typedef Interval<string> iRange;
typedef vector<interval> intervalVector;
typedef vector<iRange> rangeVector;
typedef IntervalTree< vector<string> > intervalTree;
typedef map<string, intervalTree> IntervalForest;


SEXP cigar_ranges(SEXP cigar, SEXP flag, SEXP space, SEXP pos, SEXP f,
		  SEXP ops, SEXP drop_empty_ranges, SEXP reduce_ranges,
		  SEXP with_ops);


#endif
