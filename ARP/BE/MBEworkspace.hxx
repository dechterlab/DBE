#ifndef MBEworkspace_HXX_INCLUDED
#define MBEworkspace_HXX_INCLUDED

#include <inttypes.h>

#include "AVLtreeInt64.hxx"
#include "Function.hxx"
#include "Problem.hxx"
#include "Bucket.hxx"
#include "MiniBucket.hxx"
#include "Workspace.hxx"

namespace ARE { class Workspace ; }

namespace BucketElimination
{

class MBEworkspace : public ARE::Workspace
{
public:
#if defined (LINUX)
  static pthread_mutex_t stopSignalMutex;
#endif

protected :

	int32_t _nVars ; // number of variables
	int32_t *_VarList ; // list of variables in order, as an array of variable indeces; in reverse elimination order : [0] is to be eliminated last, [N-1] is to be eliminated first.
	int32_t *_VarPos ; // positions of variables in the order.
	Bucket **_Var2BucketMapping ; // bucket for each variable
//	// when keeping track of the children of each buckets, we need space; all the space for keeping track of all children of all buckets, 
//	// takes no more than _n space. allocate it at once.
//	int32_t *_BTchildlistStorage = NULL ;

	bool _DeleteUsedTables ; // delete tables used (child tables when parent bucket is computed)
	FILE *_fpLOG ;

	int32_t _VarOrdering_InducedWidth ; // induced width of the given ordering
	int32_t _PseudoWidth ; // induced width of the Mini-Bucket Elimination given current i-Bound

	int32_t _iBound ; // iBound is the max num of variables in a mini-bucket, including ones being eliminated and ones remaining; default is 1000000 = meaning infinite
	int32_t _nBucketsWithPartitioning ; // number of buckets with >1 minibuckets
	int32_t _maxDepthPartitionedBucket ; // maximum depth of a bucket with partitioning
	int32_t _MaxNumMiniBucketsPerBucket ; // maximum number of minibuckets in a bucket

	int32_t _MaxTreeHeight_BranchingVars_Limit ; // limit on the MaxTreeHeight_BranchingVars; no path from root to leaf should have more than this many branching variables. default is INT_MAX.

	int32_t _nEvidenceVars ; // number of variables with fixed value

	// memory bounds to run MBE
	double _MaxSpaceAllowed_Log10 ;

	// if the treetype is OR, we need to serialize the pseudotree; many different ways of doing that.
	bool _ForORtreetype_use_DFS_order ; // default is false

public :

	inline int32_t N(void) const { return _nVars ; }
	inline const int32_t *VarList(void) const { return _VarList ; }
	inline const int32_t *VarPos(void) const { return _VarPos ; }
	inline Bucket *MapVar2Bucket(int32_t v) const { return v >= 0 && v < _nVars ? _Var2BucketMapping[v] : NULL ; }
	inline double & MaxSpaceAllowed_Log10(void) { return _MaxSpaceAllowed_Log10 ; }
	inline bool ProblemIsLogScale(void) const { return _Problem->FunctionsAreConvertedToLogScale() ; }
	inline bool & ForORtreetype_use_DFS_order(void) { return _ForORtreetype_use_DFS_order ; }

	inline int32_t VarOrdering_InducedWidth(void) const { return _VarOrdering_InducedWidth ; }
	inline int32_t PseudoWidth(void) const { return _PseudoWidth ; }
	inline int32_t & iBound(void) { return _iBound ; }
	inline int32_t & MaxTreeHeight_BranchingVars_Limit(void) { return _MaxTreeHeight_BranchingVars_Limit ; }
	inline int32_t nEvidenceVars(void) const { return _nEvidenceVars ; }

	inline int32_t nBucketsWithPartitioning(void) const { return _nBucketsWithPartitioning ; }
	inline int32_t maxDepthPartitionedBucket(void) const { return _maxDepthPartitionedBucket ; }
	inline int32_t MaxNumMiniBucketsPerBucket(void) const { return _MaxNumMiniBucketsPerBucket ; }

	inline void SetVarOrdering_InducedWidth(int32_t w) {  _VarOrdering_InducedWidth = w ; }
	inline void SetPseudoWidth(int32_t w) {  _PseudoWidth = w ; }
	inline bool DeleteUsedTables(void) const { return _DeleteUsedTables ; }
	inline void SetDeleteUsedTables(bool v) { _DeleteUsedTables = v ; }

	inline FILE * & logFile(void) { return _fpLOG ; }

	int32_t GetFirstQueryRootVariable(void)
	{
		std::vector<int32_t> & qv = _Problem->QueryVariables() ;
		for (int32_t i = 0 ; i < qv.size() ; ++i) {
			int32_t v = qv[i] ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			if (NULL == b->ParentBucket()) 
				return v ;
			}
		return -1 ;
	}

public :

	// stop/exit computation thread
	long volatile _StopAndExit ;

#if defined WINDOWS || _WINDOWS
	uintptr_t _ThreadHandle ;
#elif defined (LINUX)
	pthread_t _ThreadHandle ;
#endif 

	int32_t CreateThread(void) ;
	int32_t StopThread(void) ;

	// ****************************************************************************************
	// Query.
	// ****************************************************************************************

protected :

	ARE_Function_TableType _AnswerFactor ; // it is in the same normal/log scale as the problem.
	int32_t _FnCombinationType ; // 0 = undef, 1 = product, 2 = sum
	int32_t _VarEliminationType ; // 0 = undef, 1 = sum, 2 = max, 3 = min

	// this is the answer when all variables are eliminated; it is in the same normal/log scale as the problem.
	ARE_Function_TableType _CompleteEliminationResult ;
	ARE_Function_TableType _CompleteEliminationResult_ub ;
	ARE_Function_TableType _CompleteEliminationResult_lb ;

	signed char _CurrentComputationBound ; // -1=min, 0=unknown, +1=max; set when ComputeOutputFunctions() is executed.

	// this is the answer when distribution on 1 variable is computed; it is in the same normal/log scale as the problem.
	int32_t _MarginalSingleVariableDistributionVar ;
	std::vector<ARE_Function_TableType> _MarginalSingleVariableDistribution ;

	int64_t _tStart, _tEnd, _tToStop, _RunTimeInMilliseconds ;

public :

	inline ARE_Function_TableType AnswerFactorEx(void) const { return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _AnswerFactor) : _AnswerFactor ; }
	inline ARE_Function_TableType AnswerFactor(void) const { return _AnswerFactor ; }
	inline ARE_Function_TableType CompleteEliminationResultEx(signed char bound) const
	{
		if (bound < 0) 
			return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _CompleteEliminationResult_lb) : _CompleteEliminationResult_lb ;
		if (bound > 0) 
			return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _CompleteEliminationResult_ub) : _CompleteEliminationResult_ub ;
		return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _CompleteEliminationResult) : _CompleteEliminationResult ;
	}
	inline ARE_Function_TableType CompleteEliminationResult(signed char bound) const 
	{ 
		if (bound < 0) 
			return _CompleteEliminationResult_lb ;
		if (bound > 0) 
			return _CompleteEliminationResult_ub ;
		return _CompleteEliminationResult ; 
	}
	inline ARE_Function_TableType MarginalSingleVariableDistributionEx(int32_t idx) const
	{
		if (idx < 0 || idx >= _MarginalSingleVariableDistribution.size()) 
			return DBL_MAX ;
		return _Problem->FunctionsAreConvertedToLogScale() ? pow(10.0, _MarginalSingleVariableDistribution[idx]) : _MarginalSingleVariableDistribution[idx] ;
	}
	inline ARE_Function_TableType MarginalSingleVariableDistribution(int32_t idx) const
	{
		if (idx < 0 || idx >= _MarginalSingleVariableDistribution.size()) 
			return DBL_MAX ;
		return _MarginalSingleVariableDistribution[idx] ;
	}
	inline std::vector<ARE_Function_TableType> & MarginalSingleVariableDistribution(void) { return _MarginalSingleVariableDistribution ; }
	inline int32_t & MarginalSingleVariableDistributionVar(void) { return _MarginalSingleVariableDistributionVar ; }
	inline int32_t FnCombinationType(void) const { return _FnCombinationType ; }
	inline int32_t VarEliminationType(void) const { return _VarEliminationType ; }
	inline signed char CurrentComputationBound(void) const { return _CurrentComputationBound ; }

	inline void AddAnswerFactor(ARE_Function_TableType & F)
	{
		ApplyFnCombinationOperator(_AnswerFactor, F) ;
	}
	inline void ApplyFnCombinationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				V += v ;
			else 
				V *= v ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V += v ;
			}
	}
	inline void ApplyFnDivisionOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				V -= v ;
			else 
				V /= v ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUB_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V -= v ;
			}
	}
	inline void ApplyVarEliminationOperator(ARE_Function_TableType & V, ARE_Function_TableType v)
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v)
			else 
				V += v ;
			}
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) 
			{ if (v > V) V = v ; }
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) 
			{ if (v < V) V = v ; }
	}
	inline ARE_Function_TableType FnCombinationNeutralValue(void) const
	{
		if (FN_COBINATION_TYPE_PROD == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return 0.0 ;
			return 1.0 ;
			}
		else if (FN_COBINATION_TYPE_SUM == _FnCombinationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return ARP_nInfinity ;
			return 0.0 ;
			}
		return DBL_MAX ;
	}
	inline ARE_Function_TableType VarEliminationDefaultValue(void) const
	{
		if (VAR_ELIMINATION_TYPE_SUM == _VarEliminationType) {
			if (_Problem->FunctionsAreConvertedToLogScale()) 
				return ARP_nInfinity ;
			return 0.0 ;
			}
		else if (VAR_ELIMINATION_TYPE_MAX == _VarEliminationType) {
			return ARP_nInfinity ;
			}
		else if (VAR_ELIMINATION_TYPE_MIN == _VarEliminationType) {
			return ARP_pInfinity ;
			}
		return DBL_MAX ;
	}
	inline void SetFnCombinationType(int32_t v) { _FnCombinationType = v ; }
	inline void SetVarEliminationType(int32_t v) { _VarEliminationType = v ; }
	inline void SetCompleteEliminationResult(ARE_Function_TableType v, signed char bound)
	{
		if (bound < 0) 
			_CompleteEliminationResult_lb = v ;
		else if (bound > 0) 
			_CompleteEliminationResult_ub = v ;
		else
			_CompleteEliminationResult = v ;
	}

	inline int64_t & tStart(void) { return _tStart ; }
	inline int64_t & tEnd(void) { return _tEnd ; }
	inline int64_t & tToStop(void) { return _tToStop ; }
	inline int64_t & RunTimeInMilliseconds(void) { return _RunTimeInMilliseconds ; }

	// ****************************************************************************************
	// Bucket-tree structure.
	// ****************************************************************************************

protected :

	int32_t _nBuckets ;
	Bucket **_Buckets ; // note that _Buckets[] array size may be smaller than _Var2BucketMapping[], but the order must be the same.

	// an ordered list of buckets to compute; typically this list is ordered in the order of decreasing height (of the bucket), 
	// and the computation is carried out from last-to-first.
	// the allocated size of this array is _nVars.
	int32_t _lenComputationOrder ;
	int32_t *_BucketOrderToCompute ;

	int32_t _nBucketsWithSingleChild_initial ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int32_t _nBucketsWithSingleChild_final ; // number of buckets with a single child; these buckets can be eliminated by using superbuckets.
	int32_t _nBuckets_initial ;
	int32_t _nBuckets_final ;
	int32_t _nVarsWithoutBucket ;
	int32_t _nConstValueFunctions ;

	int32_t _MaxNumChildren ; // maximum number of children of any bucket
	int32_t _MaxNumVarsInBucket ; // maximum number of variables in any bucket
	int32_t _MaxTreeHeight ; // maximum height of the bucket tree; as number of edges to travel.
	int32_t _MaxTreeHeight_BranchingVars ; // maximum height of the bucket tree when counting branching variables; as number of branching variables on the path.
	int32_t _MaxBucketFunctionWidth ; // maximum size of any bucket output function
	int32_t _nBucketsWithSingleChild ; // number of buckets with single child
	int32_t _nBucketsWithNoChildren ; // number of leaf buckets
	int32_t _nRoots ; // number of roots in the bucket tree
	int64_t _TotalOriginalFunctionSize ; // size of original functions in existance during any point during BE execution
	int64_t _TotalOriginalFunctionSpace ; // space of original functions in existance during any point during BE execution
	double _TotalNewFunctionSize_Log10 ; // return number of entries over all tables generated during BE execution
	double _TotalNewFunctionSpace_Log10 ; // return space (in bytes) over all tables generated during BE execution
	double _TotalNewFunctionComputationComplexity_Log10 ; // return bucket_width*bucket_nFNs over all buckets; this is the number of operations it takes to compute BE execution.
	double _MaxSimultaneousNewFunctionSize_Log10 ; // max size of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousNewFunctionSpace_Log10 ; // max space of new functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousTotalFunctionSize_Log10 ; // max size of all (old/new) functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.
	double _MaxSimultaneousTotalFunctionSpace_Log10 ; // max space of all (old/new) functions in existance during any point during BE execution; this is determined by DeleteUsedTables flag.

	double _TotalNewFunctionSizeComputed_Log10 ; // statistics computed during the execution

	// for all pairs of variables that are on some path from root to a leaf, store a key (ancestor<<32) + descendant.
	// we can use this to check if some linearization of the BT is correct.
	// this is computed from original bucket tree corresponding to the given variable order.
	CMauiAVL64Tree _BTpartialorder ;

public :

	inline int32_t nBuckets(void) const { return _nBuckets ; }
	inline BucketElimination::Bucket *getBucket(int32_t IDX) const { return NULL != _Buckets ? _Buckets[IDX] : NULL ; }
	inline int32_t lenComputationOrder(void) const { return _lenComputationOrder ; }
	inline int32_t BucketOrderToCompute(int32_t IDX) const { return _BucketOrderToCompute[IDX] ; }
	inline int32_t Height(void)
	{
		int32_t height = -1 ;
		if (_nBuckets <= 0 || NULL == _Buckets) 
			return -1 ;
		for (int32_t i = 0 ; i < _nBuckets ; ++i) {
			if (_Buckets[i]->Height() > height) 
				height = _Buckets[i]->Height() ;
			}
		return height ;
		// return _Buckets[0]->Height() ;
	}

	inline int32_t MaxNumChildren(void) const { return _MaxNumChildren ; }
	inline int32_t MaxNumVarsInBucket(void) const { return _MaxNumVarsInBucket ; }
	inline int32_t MaxTreeHeight(void) const { return _MaxTreeHeight ; }
	inline int32_t nBucketsWithSingleChild_initial(void) const { return _nBucketsWithSingleChild_initial ; }
	inline int32_t nBucketsWithSingleChild_final(void) const { return _nBucketsWithSingleChild_final ; }
	inline int32_t nBuckets_initial(void) const { return _nBuckets_initial ; }
	inline int32_t nBuckets_final(void) const { return _nBuckets_final ; }
	inline int32_t nVarsWithoutBucket(void) const { return _nVarsWithoutBucket ; }
	inline int32_t nConstValueFunctions(void) const { return _nConstValueFunctions ; }
	inline int32_t MaxBucketFunctionWidth(void) const { return _MaxBucketFunctionWidth ; }
	inline int32_t nBucketsWithSingleChild(void) const { return _nBucketsWithSingleChild ; }
	inline int32_t nBucketsWithNoChildren(void) const { return _nBucketsWithNoChildren ; }
	inline int32_t nRoots(void) const { return _nRoots ; }
	inline int64_t TotalOriginalFunctionSize(void) const { return _TotalOriginalFunctionSize ; }
	inline int64_t TotalOriginalFunctionSpace(void) const { return _TotalOriginalFunctionSpace ; }
	inline double TotalNewFunctionSize_Log10(void) const { return _TotalNewFunctionSize_Log10 ; }
	inline double TotalNewFunctionSpace_Log10(void) const { return _TotalNewFunctionSpace_Log10 ; }
	inline double TotalNewFunctionComputationComplexity_Log10(void) const { return _TotalNewFunctionComputationComplexity_Log10 ; }
	inline double MaxSimultaneousNewFunctionSize_Log10(void) const { return _MaxSimultaneousNewFunctionSize_Log10 ; }
	inline double MaxSimultaneousNewFunctionSpace_Log10(void) const { return _MaxSimultaneousNewFunctionSpace_Log10 ; }
	inline double MaxSimultaneousTotalFunctionSize_Log10(void) const { return _MaxSimultaneousTotalFunctionSize_Log10 ; }
	inline double MaxSimultaneousTotalFunctionSpace_Log10(void) const { return _MaxSimultaneousTotalFunctionSpace_Log10 ; }

	inline double & TotalNewFunctionSizeComputed_Log10(void) { return _TotalNewFunctionSizeComputed_Log10 ; }
	inline double GetSolutionCompletionPercentage(void) { return _TotalNewFunctionSize_Log10 >= 0 ? 100.0*pow(10.0, _TotalNewFunctionSizeComputed_Log10 - _TotalNewFunctionSize_Log10) : -1.0 ; }

	// these functions can be executed immediately, as soon as ws is initialled.
	int64_t ComputeTotalOriginalFunctionSizeAndSpace(void) ;

	// this fn can be run any time, whether MB partitioning is done or not
	int32_t ComputeMaxNumVarsInBucket(bool ForceComputeSignature = false) ;

	// this fn can be run any time, whether MB partitioning is done or not
	BucketElimination::Bucket *GetBucketWithMostVariables(void) ;

	// Compute num of roots
	int32_t ComputeNumRoots(void)
	{
		_nRoots = _nEvidenceVars = 0 ;
		for (int32_t i = 0 ; i < _Problem->N() ; ++i) 
			{ if (_Problem->IsEvidenceVar(i)) _nEvidenceVars++ ; }
		for (int32_t i = 0 ; i < _nBuckets ; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (NULL == b) continue ;
//			if (_Problem->IsEvidenceVar(b->V())) 
//				continue ; // don't count evidence variables
			if (NULL == b->ParentBucket()) 
				_nRoots++ ;
			}
		return 0 ;
	}
	int32_t GetRoots(std::vector<BucketElimination::Bucket*> & Roots)
	{
		for (int32_t i = 0 ; i < _nBuckets ; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (NULL == b) continue ;
			if (NULL == b->ParentBucket()) 
				Roots.push_back(b) ;
			}
		return 0 ;
	}

	int32_t GetPseudotreeChildrenOfVariable(int32_t V, std::vector<BucketElimination::Bucket*> & Children)
	{
		Children.clear() ;
		if (V >= 0) {
			BucketElimination::Bucket *b = MapVar2Bucket(V) ;
			if (NULL != b) 
				b->GetChildren(Children) ;
			}
		else {
			GetRoots(Children) ;
			}
		return 0 ;
	}

	// usually below the root there is a chain, length of which is roughly w*
	int32_t GetRootChainLength(void)
	{
		int32_t H = 0 ;
		std::vector<BucketElimination::Bucket*> roots ;
		GetRoots(roots) ;
		for (int32_t i = 0 ; i < roots.size() ; ++i) {
			int32_t h = 0 ;
			for (BucketElimination::Bucket *b = roots[i] ; NULL != b ? 1 == b->nChildren() : false ; b = MapVar2Bucket(b->ChildVar(0)), ++h) ;
			if (h > H) H = h ;
			}
		return H ;
	}

	int32_t ComputeDistanceToNextDescendantBranchingVariable(int32_t Var) 
	// number of edges we need to travel to reach closest branching variable descendant; if return value is <0, then there is no descendant branching variable (we have a chain that ends in a leaf node).
	{
		BucketElimination::Bucket *b = MapVar2Bucket(Var) ;
		if (NULL == b) 
			return INT_MAX ;
		int32_t n = 0 ;
		for (; NULL != b ? 1 == b->nChildren() : false ; ++n) {
			int32_t v = b->ChildVar(0) ;
			if (v < 0) 
				return -n ;
			b = MapVar2Bucket(v) ;
			}
		if (NULL == b) return -n ; // this should not happen
		return b->nChildren() <= 0 ? -n : n ;
	}

	// these functions can be executed once bucket tree structure is set up; no MB partitioning is required.
	// these values are wrt original (full) bucket tree.
	int32_t ComputeMaxNumChildren(void) ;
	int32_t ComputeNBucketsWithSingleChild(void) ;

	// these functions can be executed once MB partitioning is done.
	int32_t ComputeMaxBucketFunctionWidth(void) ;
	void ComputeTotalNewFunctionSizeAndSpace(void) ;
	void ComputeTotalNewFunctionComputationComplexity(void) ;

	// check that current BT satisfies stored BT partial order constraints.
	int32_t CheckBTpartialorder(void) 
	{
		for (int32_t i = 0 ; i < _nBuckets ; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (NULL == b) continue ;
			int32_t vC = b->V() ;
			for (BucketElimination::Bucket *b_ = b->ParentBucket() ; NULL != b_ ; b_ = b_->ParentBucket()) {
				int32_t vA = b_->V() ;
				int64_t key = vA ; key <<= 32 ; key += vC ;
				int res1 = _BTpartialorder.Find(key, NULL ) ;
				int64_t key2 = vC ; key2 <<= 32 ; key2 += vA ;
				int res2 = _BTpartialorder.Find(key2, NULL ) ;
				if (0 != res2) 
					return 1 ;
				}
			}
		return 0 ;
	}

	// **************************************************************************************************
	// Setting up workspace/bucket-tree
	// **************************************************************************************************

public :

	virtual int32_t Initialize(ARE::ARP & Problem, bool UseLogScale, const int32_t *VarOrderingAsVarList, int32_t DeleteUsedTables) ;

	// allocate buckets; assign given/input/original functions to their buckets.
	// this function does no MB partitioning!
	virtual int32_t CreateBuckets(bool ANDORtree, bool KeepBTsignature, bool SimplifyBTstructure, bool CreateSuperBuckets) ;

	// compute various quantities, e.g. height of each bucket, distance to root, number of subtree buckets, etc.
	virtual int32_t ComputeBucketTreeStatistics(void)
	{
		int32_t i ;

		_MaxTreeHeight = 0 ;
		_MaxTreeHeight_BranchingVars = 0 ;

		INT64 nAdded = 0 ;
		BucketElimination::Bucket *BLhead = NULL, *BLtail = NULL ;

		// collect all roots; reset stuff
		for (i = 0 ; i < _nBuckets; i++) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			BucketElimination::Bucket *B = b->ParentBucket() ;
			b->SetMaxDescendantNumVars(-1) ;
			b->SetMaxTreeHeight_BranchingVars(0) ;
			b->SetHeight(-1) ;
			b->SetnSubtreeBuckets(1) ;
			if (NULL == B) {
				if (NULL == BLhead) 
					{ BLhead = BLtail = b ; }
				else 
					{ BLtail->NextInOrder() = b ; BLtail = b ; }
				b->NextInOrder() = NULL ;
				++nAdded ;
				}
			}
		int32_t nRoots = nAdded ;

		// add all other buckets to the list, so that parent is before children
		BucketElimination::Bucket *B = BLhead ;
		while (NULL != B) {
			for (int32_t j = 0 ; j < B->nChildren() ; j++) {
				int32_t childvar = B->ChildVar(j) ;
				BucketElimination::Bucket *b = MapVar2Bucket(childvar) ;
				if (NULL == b) continue ;
				BLtail->NextInOrder() = b ; BLtail = b ; b->NextInOrder() = NULL ;
				++nAdded ;
				}
			B = B->NextInOrder() ;
			}
		if (nAdded != _nBuckets) 
			return 1 ;

		// set : root ptr; distance2root
		for (BucketElimination::Bucket *b = BLhead ; NULL != b ; b = b->NextInOrder()) {
			B = b->ParentBucket() ;
			if (NULL == B) {
				b->SetDistanceToRoot(0) ;
				b->SetRootBucket(b) ;
				}
			else {
				b->SetDistanceToRoot(B->DistanceToRoot() + 1) ;
				b->SetRootBucket(B->RootBucket()) ;
				}
			}

		// reverse the list so that we can go from leaves to the root
		BLtail = BLhead ;
		BucketElimination::Bucket *bnext = NULL ;
		for (BucketElimination::Bucket *b = BLhead->NextInOrder() ; NULL != b ; b = bnext) {
			bnext = b->NextInOrder() ;
			b->NextInOrder() = BLhead ;
			BLhead = b ;
			}
		BLtail->NextInOrder() = NULL ;

		// set : height; nSubTreeBuckets; MaxTreeHeight_BranchingVars; MaxDescendantNumVars
		for (BucketElimination::Bucket *b = BLhead ; NULL != b ; b = b->NextInOrder()) {
			B = b->ParentBucket() ;
			if (b->Height() < 0) 
				b->SetHeight(0) ;
			if (b->nChildren() > 1) 
				b->SetMaxTreeHeight_BranchingVars(1 + b->MaxTreeHeight_BranchingVars()) ;
			if (_MaxTreeHeight_BranchingVars < b->MaxTreeHeight_BranchingVars()) 
				_MaxTreeHeight_BranchingVars = b->MaxTreeHeight_BranchingVars() ;
			if (NULL != B) {
				int32_t h = b->Height() + 1 ;
				if (B->Height() < h) 
					B->SetHeight(h) ;
				if (_MaxTreeHeight < h) 
					_MaxTreeHeight = h ;
				if (B->MaxTreeHeight_BranchingVars() < b->MaxTreeHeight_BranchingVars()) 
					B->SetMaxTreeHeight_BranchingVars(b->MaxTreeHeight_BranchingVars()) ;
				B->SetnSubtreeBuckets(B->nSubtreeBuckets() + b->nSubtreeBuckets()) ;
				int32_t mdnv = b->MaxDescendantNumVarsEx() ;
				if (B->MaxDescendantNumVars() < mdnv) 
					B->SetMaxDescendantNumVars(mdnv) ;
				}
			}

		_nBucketsWithSingleChild_final = ComputeNBucketsWithSingleChild() ;

		ComputeMaxNumChildren() ;

		return 0 ;
	}

	int32_t Enfoce_MaxTreeHeight_BranchingVars_Limit(void)
	{
		if (_MaxTreeHeight_BranchingVars_Limit < 0) 
			return 0 ; // limit not specified; do nothing.
		if (_MaxTreeHeight_BranchingVars < 0) {
			if (0 != ComputeBucketTreeStatistics()) 
				return 1 ;
			}

		// get roots
		std::vector<BucketElimination::Bucket*> roots ;
		GetRoots(roots) ;
		if (0 == roots.size()) return 0 ;
		// filter out those roots that are already fine
		for (int32_t i = roots.size()-1 ; i >= 0 ; --i) {
			if (roots[i]->MaxTreeHeight_BranchingVars() <= _MaxTreeHeight_BranchingVars_Limit) 
				roots.erase(roots.begin() + i) ;
			}

		// temp array to compute DFS ordering of subtrees
		std::vector<BucketElimination::Bucket*> dfs_reverse_chain ; dfs_reverse_chain.reserve(_nBuckets) ;
		std::vector<BucketElimination::Bucket*> dfs_stack_b ; dfs_stack_b.reserve(_nBuckets) ;
		std::vector<int32_t> dfs_stack_i ; dfs_stack_i.reserve(_nBuckets) ;

		// need to bring each disconnected component (tree) into compliance
		std::vector<BucketElimination::Bucket*> chain ; chain.reserve(_Problem->N()) ;
		for (int32_t i = 0 ; i < roots.size() ; ++i) {
			BucketElimination::Bucket *root = roots[i] ;
			while (root->MaxTreeHeight_BranchingVars() > _MaxTreeHeight_BranchingVars_Limit) {
				BucketElimination::Bucket *b = root ;
				chain.clear() ; chain.push_back(b) ;
				int nBranchings = 0 ;
				// within this tree, build a chain from root to a leaf that contains maxNumBranchingVars
				while (NULL != b) {
					// assume b is already in the chain; pick which of its children goes in the chain
					if (b->nChildren() <= 0) 
						break ;
					if (b->nChildren() <= 1) {
						BucketElimination::Bucket *b_ = MapVar2Bucket(b->ChildVar(0)) ;
						if (NULL == b_) 
							break ; // error
						chain.push_back(b_) ; b = b_ ; }
					else { // b is a branching variable
						++nBranchings ;
						// if we decide to eliminate this branching, decide which one to keep
						BucketElimination::Bucket *b_c = NULL ;
						for (int32_t j = 0 ; j < b->nChildren() ; ++j) {
							BucketElimination::Bucket *b_ = MapVar2Bucket(b->ChildVar(j)) ;
							if (NULL == b_) continue ;
							if (NULL == b_c) { b_c = b_ ; continue ; }
							if (b_->MaxTreeHeight_BranchingVars() < b_c->MaxTreeHeight_BranchingVars()) 
								{ continue ; }
							if (b_->MaxTreeHeight_BranchingVars() > b_c->MaxTreeHeight_BranchingVars()) 
								{ b_c = b_ ; continue ; }
							// keep the one with more total nodes
							if (b_c->nSubtreeBuckets() < b_->nSubtreeBuckets()) 
								{ b_c = b_ ; continue ; }
							}
						if (NULL == b_c) 
							break ; // error
						chain.push_back(b_c) ;
						b = b_c ;
						}
					}
				// 'chain' is a chain of nodes, from root to a leaf. there must be exactly root->MaxTreeHeight_BranchingVars() branching variables in this chain.
				if (nBranchings != root->MaxTreeHeight_BranchingVars()) {
					int error = 1 ;
					}
				// pick the branching var with fewest nodes to be eliminated.
				int32_t idxChainToEliminate = -1, sum = -1 ;
				for (int32_t j = 0 ; j < chain.size() ; ++j) {
					BucketElimination::Bucket *b = chain[j] ;
					if (b->nChildren() <= 1) continue ;
					int32_t sumBuckets = 0 ;
					for (int32_t k = 0 ; k < b->nChildren() ; ++k) {
						BucketElimination::Bucket *b_ = MapVar2Bucket(b->ChildVar(k)) ;
						if (NULL == b_) 
							continue ; // error
						if (b_ == chain[j+1]) continue ;
						sumBuckets += b_->nSubtreeBuckets() ;
						}
					if (idxChainToEliminate < 0 ? true : sum > sumBuckets) 
						{ idxChainToEliminate = j ; sum = sumBuckets ; }
					}
				// we have picked the branching variable to eliminate; merge subtrees outside the chain into the chain
				int chain_len = chain.size() ;
				BucketElimination::Bucket *b2Process = chain[idxChainToEliminate] ;
				BucketElimination::Bucket *b2Process_child = chain[idxChainToEliminate+1] ;
				dfs_reverse_chain.clear() ; // collect all other descendants into this array
				for (int32_t k = 0 ; k < b2Process->nChildren() ; ++k) {
					BucketElimination::Bucket *b_ = MapVar2Bucket(b2Process->ChildVar(k)) ;
					if (NULL == b_) continue ;
					if (b_ == b2Process_child) continue ;
					int32_t res = b_->GetDFSorderofsubtree(dfs_reverse_chain, dfs_stack_b, dfs_stack_i) ;
					}
				// insert dfs_reverse_chain nodes in same order between b2Process and b2Process_child
				int32_t children[1] ;
				if (dfs_reverse_chain.size() != sum) {
					return 1 ;
					}
				for (int32_t k = 0 ; k <= dfs_reverse_chain.size() ; ++k) {
					BucketElimination::Bucket *bC = k > 0 ? dfs_reverse_chain[k-1] : b2Process_child ;
					BucketElimination::Bucket *bP = k < dfs_reverse_chain.size() ? dfs_reverse_chain[k] : b2Process ;
					children[0] = bC->V() ;
					bP->SetChildren(1, children) ;
					bC->SetParentBucket(bP) ;
					}
				// TODO : there may be more efficient way of recomputing MaxTreeHeight_BranchingVars
				if (0 != ComputeBucketTreeStatistics()) 
					return 1 ;
				}
			}
		// need to update bucket tree stats
		if (0 != ComputeBucketTreeStatistics()) 
			return 1 ;

		return 0 ;
	}

	// if ReducedProblemType==1 :
	//		apply BE for all buckets, bottom up, recursively, for which the bucket size is 3. i.e. resulting output function is over 2 variables.
	//		at the end, all leaves of the resulting BT have width > 3.
	//		save this reduced (but equivalent to the original problem) in a uai format file; also new_var->old_var mapping.
	//		note : this fn will destroy current MB partitioning.
	// if ReducedProblemType==2 : 
	//		keep query variables
	// if ReducedProblemType<=0 :
	//		copy input to output
	int32_t SaveReducedProblem(
		// IN
		signed char ReducedProblemType, signed char ApproximationBound, std::string & fn, 
		// OUT
		int32_t & nKeepVariables, std::vector<int32_t> & Old2NewVarMap, std::vector<ARE::Function*> & ReducedProblemFunctions) ;

	// find largest i-bound so that new function space is within given space limit.
	// if the function succeeds, upon return, there is a mb partitioning with the i-bound found.
	int32_t FindIBoundForSpaceAllowed(int32_t MinIBound, int32_t & bestIboundFound, double & NewFnSpaceUsed_Log10, int32_t & nBucketsPartitioned, int32_t & maxDepthPartitionedBucket) ;

	// compute output functions of all buckets
	int32_t ComputeOutputFunctions(bool DoMomentMatching, signed char Bound, int64_t & TotalSumOutputFunctionsNumEntries) ;

	// delete all MB generated tables
	int32_t DeleteMBEgeneratedTables(void) ;

	// create an order in which buckets should be processed
	// algorithm=0 means from leaves to root in terms of uniform height, i.e. one level must be finished before starting next level. this is default.
	// algorithm=1 means in the order that minimizes space, e.g. finish computing one bucket before starting next sibling.
	// MB partitioning should be done when this fn is called.
	virtual int32_t CreateComputationOrder(int32_t algorithm) ;

	// compute/set the ArgumentsPermutationList of all functions, such that any arguments index is its distance from root.
	// this is useful when processing a search space, where each node carries along its path assignment (from root to itself);
	// then to evaluate the function, we can use the path-assignment and ArgumentsPermutationList[] of each function in the bucket.
	int32_t ComputeFunctionArgumentsPermutationList(void) ;

	// create/destroy MB partitionin.
	int32_t CreateMBPartitioning(bool CreateTables, bool doMomentMatching, signed char Bound, int32_t ComputeComputationOrder) ;
	int32_t DestroyMBPartitioning(void) ;

	int32_t CheckBucketTreeIntegrity(void) ;
	
	// assuming entire MBE elimination is completed, compute some data
	int32_t PostComputationProcessing(void) ;

	// check that each augmented functions contains the buckets variable
	int32_t VerifyOriginalFunctions(void)
	{
		int32_t i ;
		// create DFS order of the whole problem
		std::vector<BucketElimination::Bucket*> dfs_reverse_chain ; dfs_reverse_chain.reserve(_nBuckets) ;
		std::vector<BucketElimination::Bucket*> dfs_stack_b ; dfs_stack_b.reserve(_nBuckets) ;
		std::vector<int32_t> dfs_stack_i ; dfs_stack_i.reserve(_nBuckets) ;
		std::vector<BucketElimination::Bucket*> roots ;
		for (i = 0 ; i < _nBuckets ; i++) 
			{ BucketElimination::Bucket *b = _Buckets[i] ; if (NULL == b) continue ; if (NULL == b->ParentBucket()) roots.push_back(b) ; }
		for (i = 0 ; i < roots.size() ; ++i) 
			roots[i]->GetDFSorderofsubtree(dfs_reverse_chain, dfs_stack_b, dfs_stack_i) ;
		// check all functions against the DFS order
		for (i = 0 ; i < _Problem->nFunctions() ; ++i) {
			ARE::Function *f = _Problem->getFunction(i) ;
			if (NULL == f) continue ;
			int32_t seen_this_fn = 0 ;
			for (int32_t j = 0 ; j < dfs_reverse_chain.size() ; ++j) {
				BucketElimination::Bucket *b = dfs_reverse_chain[j] ;
				int32_t v = b->V() ;
				if (f->ContainsVariable(v)) {
					int32_t idx = b->ContainsOriginalFunction(*f) ;
					if (0 == seen_this_fn) {
						if (idx < 0) 
							return 1 ; // this is the earliest bucket that contains a var in the fn; this bucket should have this fn in it.
						++seen_this_fn ;
						}
					else {
						if (idx >= 0) 
							return 1 ;
						}
					}
				else {
					int32_t idx = b->ContainsOriginalFunction(*f) ;
					if (idx >= 0) 
						return 1 ; // this bucket's var is not in the fn; it should not contain the fn.
					}
				}
			}
		return 0 ;
	}

	// check that each augmented functions contains the buckets variable
	int32_t VerifyAugmentedFunctions(void)
	{
		int32_t i ;
		// check each bucket is ok
		for (i = 0 ; i < _nBuckets ; ++i) {
			BucketElimination::Bucket *b = _Buckets[i] ;
			if (NULL == b) continue ;
			if (0 != b->VerifyAugmentedFunctions()) 
				return 1 ;
			}
		return 0 ;
	}

	// **************************************************************************************************
	// Executing (M)BE.
	// **************************************************************************************************

public :

	// simulate execution and compute the minimum amount of size/space required (for new functions).
	// this depends on the DeleteUsedTables flag.
	void SimulateComputationAndComputeMinSpace(bool IgnoreInMemoryTables) ;

	// run regular (M)BE on the bucket-tree
	virtual int32_t RunSimple(void) ;

	// extract solution assignment from current MBE execution; operator (min/max) will be obtained from the problem.
	// solution will be stored in the problem.
	int32_t BuildSolution(void) ;

	// **************************************************************************************************
	// Problem generation
	// **************************************************************************************************

public :

	// this function generates a random uniform Bayesian network for the given parameters.
	// it will not fill in any function tables.
	int32_t GenerateRandomBayesianNetworkStructure(
		int32_t N, // # of variables
		int32_t K, // same domain size for all variables
		int32_t P, // # of parents per CPT
		int32_t C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
		int32_t ProblemCharacteristic // 0 = totally random, 1 = 1 leaf node
		) ;

public :

	virtual int32_t Destroy(void) ;
	virtual int32_t DestroyBucketPartitioning(void) ;

	MBEworkspace(const char *BEEMDiskSpaceDirectory = NULL) ;
	virtual ~MBEworkspace(void)
	{
		Destroy() ;
	}
} ;

// this function generates given number of random uniform Bayesian networks for the given parameters, 
// within specified complexity bounds.
// problems will be saved in uai format.
// return value is number of problem generated/saved.
int32_t GenerateRandomBayesianNetworksWithGivenComplexity(
	int32_t nProblems, 
	int32_t N, // # of variables
	int32_t K, // same domain size for all variables
	int32_t P, // # of parents per CPT
	int32_t C,  // # of CPTs; variables that are not a child in a CPT will get a prior.
	int32_t ProblemCharacteristic, // 0 = totally random, 1 = 1 leaf node
	int64_t MinSpace, 
	int64_t MaxSpace
	) ;

} // namespace BucketElimination

#endif // MBEworkspace_HXX_INCLUDED
