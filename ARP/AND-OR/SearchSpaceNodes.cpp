#include "MBEworkspace.hxx"

#include <SearchSpaceNodes.hxx>
#include <AbsSamplingWorkspace.hxx>

//int64_t AndOrSearchSpace::nNODES__ = 0 ;

AndOrSearchSpace::SearchAndNode::SearchAndNode(void) : 
	SearchNode(NULL, -1), _Value(-1), _VarPosInOrder(-1), _ClosestBranchingAncestor(NULL), _ClosestBranchingAncestorVarChildIdx(-1), _ExpansionQueueKey(0), 
	_h(ARP_DBL_MAX), _Cost(ARP_DBL_MAX), _AccumulatedCost(ARP_DBL_MAX), _AccumulatedCostSinceLastBranchingVariable(ARP_DBL_MAX), _ISw(ARP_DBL_MAX), _AccumulatedISw(ARP_DBL_MAX), _ISq(ARP_DBL_MAX), _AccumulatedISwSinceLastBranchingVariable(ARP_DBL_MAX), 
	_nISmerges(1)
{
//	if (NULL != _WS) _SubTreeValue = _WS->FnCombinationNeutralValue() ;
}

AndOrSearchSpace::SearchAndNode::SearchAndNode(AbsSamplingWorkspace *WS, int Var, int Value) : 
	SearchNode(WS, Var), _Value(Value), _ClosestBranchingAncestor(NULL), _ClosestBranchingAncestorVarChildIdx(-1), _ExpansionQueueKey(0), 
	_h(WS->FnCombinationNeutralValue()), _Cost(WS->FnCombinationNeutralValue()), _AccumulatedCost(WS->FnCombinationNeutralValue()), _AccumulatedCostSinceLastBranchingVariable(_WS->FnCombinationNeutralValue()), _ISw(WS->FnCombinationNeutralValue()), _AccumulatedISw(WS->FnCombinationNeutralValue()), _AccumulatedISwSinceLastBranchingVariable(WS->FnCombinationNeutralValue()), _ISq(WS->FnCombinationNeutralValue()), 
	_nISmerges(1)
{
	if (NULL != WS && Var >= 0) 
		_VarPosInOrder = WS->Problem()->VarOrdering_VarPos()[Var] ;
	else 
		_VarPosInOrder = -1 ;
//	if (NULL != _WS) _SubTreeValue = _WS->FnCombinationNeutralValue() ;
}

AndOrSearchSpace::SearchAndNode::SearchAndNode(SearchAndNode & N) : 
	SearchNode(N), _Value(N._Value), 
	_ClosestBranchingAncestor(N._ClosestBranchingAncestor), _ClosestBranchingAncestorVarChildIdx(N._ClosestBranchingAncestorVarChildIdx), _ExpansionQueueKey(N._ExpansionQueueKey), 
	_h(N._h), _Cost(N._Cost), _AccumulatedCost(N._AccumulatedCost), _AccumulatedCostSinceLastBranchingVariable(N._AccumulatedCostSinceLastBranchingVariable), _ISw(N._ISw), _AccumulatedISw(N._AccumulatedISw), _AccumulatedISwSinceLastBranchingVariable(N._AccumulatedISwSinceLastBranchingVariable), _ISq(N._ISq), 
	_nISmerges(0)
{
	if (NULL != _WS && _V >= 0) 
		_VarPosInOrder = _WS->Problem()->VarOrdering_VarPos()[_V] ;
	else 
		_VarPosInOrder = -1 ;
//	if (NULL != _WS) _SubTreeValue = _WS->FnCombinationNeutralValue() ;
}

AndOrSearchSpace::SearchAndNode_WithPath::SearchAndNode_WithPath(AndOrSearchSpace::AbsSamplingWorkspace *WS, int32_t Var, int32_t Value) 
	: SearchAndNode(WS, Var, Value), _Qpath(_WS->FnCombinationNeutralValue())  { /* ++nNODES__ ; */ }
AndOrSearchSpace::SearchAndNode_WithPath::SearchAndNode_WithPath(AndOrSearchSpace::SearchAndNode_WithPath & N) 
	: SearchAndNode(N), _Qpath(_WS->FnCombinationNeutralValue())  { /* ++nNODES__ ; */ }

BucketElimination::Bucket *AndOrSearchSpace::SearchNode::Bucket(void) const
	{ return _V >= 0 ? _WS->MapVar2Bucket(_V) : NULL ; }

void AndOrSearchSpace::SearchAndNode_WithPath::ApplyPathAssignment(bool SetFullPathAssignment)
{
	SearchAndNode_WithPath *node = this ;
again :
	std::vector<int32_t> & pa = node->_PathAssignment ;
	int32_t n = pa.size() ;
	BucketElimination::Bucket *b = ((BucketElimination::MBEworkspace *) _WS)->MapVar2Bucket(node->_V) ;
	for (; n > 0 && NULL != b ;) {
		int32_t value = pa[--n] ;
		_WS->SetCurrentContextValue(b->V(), value) ;
		b = b->ParentBucket() ;
		}
	if (SetFullPathAssignment) {
		if (NULL != node->_ClosestBranchingAncestor) { 
			node = (SearchAndNode_WithPath *) node->_ClosestBranchingAncestor ;
			goto again ; }
		}
}


int32_t AndOrSearchSpace::SearchAndNode_WithPath::InitializeBranchingVarChildrenSubTreeValues(void)
{
	_BranchingVarChildrenSubTreeValues.clear() ;
	BucketElimination::Bucket *B = Bucket() ;
	if (NULL == B) 
		return 0 ;
	int32_t n = B->nChildren() ;
	if (n <= 1) // this can only be computed if nChildren>1
		return 0 ;
	// allocate space/initialize
	if (_BranchingVarChildrenSubTreeValues.capacity() < n) {
		_BranchingVarChildrenSubTreeValues.reserve(n) ;
		if (_BranchingVarChildrenSubTreeValues.capacity() < n) 
			return 1 ;
		}
	_BranchingVarChildrenSubTreeValues.resize(B->nChildren()) ;
	for (int32_t i = 0 ; i < B->nChildren() ; ++i) 
		_BranchingVarChildrenSubTreeValues[i] = _WS->FnCombinationNeutralValue() ;
//	// set path assignment so that function computations have correct context
//	ApplyPathAssignment(true) ;
	// process each function that goes into h; they are Int functions of each child, as well as output functions of all minibuckets in each child
	for (int32_t i = 0 ; i < B->nChildren() ; ++i) {
		BucketElimination::Bucket *b = B->ChildBucket(i) ;
		if (NULL == b) continue ;
		for (int32_t j = 0 ; j < b->nIntermediateFunctions() ; ++j) {
			ARE::Function *f = b->IntermediateFunction(j) ;
			if (NULL == f) continue ;
			ARE_Function_TableType v = f->TableEntry(_PathAssignment.data(), _WS->Problem()->K()) ;
			_WS->ApplyFnCombinationOperator(_BranchingVarChildrenSubTreeValues[i], v) ;
			}
		std::vector<BucketElimination::MiniBucket*> & MBs = b->MiniBuckets() ;
		for (BucketElimination::MiniBucket *mb : MBs) {
			if (NULL == mb) continue ;
			ARE::Function & f = mb->OutputFunction() ;
			ARE_Function_TableType v = f.TableEntry(_PathAssignment.data(), _WS->Problem()->K()) ;
			_WS->ApplyFnCombinationOperator(_BranchingVarChildrenSubTreeValues[i], v) ;
			}
		}
	return 0 ;
}


double AndOrSearchSpace::SearchAndNode_WithPath::ComputeQPathLocal(int idxToExclude) 
{
	double v = _WS->FnCombinationNeutralValue() ;
	for (int32_t i = 0 ; i < _BranchingVarChildrenSubTreeValues.size() ; ++i) {
		if (i == idxToExclude) 
			continue ;
		_WS->ApplyFnCombinationOperator(v, _BranchingVarChildrenSubTreeValues[i]) ;
		}
	return v ;
}


int32_t AndOrSearchSpace::SearchAndNode_WithPath::PrepForDFSExpansion(AndOrSearchSpace::SearchAndNode_WithPath & N, int32_t ChildIdx)
{
	_SubTreeValue = ARP_DBL_MAX ;
	_ClosestBranchingAncestor = &N ;
	_ClosestBranchingAncestorVarChildIdx = ChildIdx ;

	_AccumulatedCostSinceLastBranchingVariable = _WS->FnCombinationNeutralValue() ;
	_AccumulatedISwSinceLastBranchingVariable = _WS->FnCombinationNeutralValue() ;
	_ISq = ARP_nInfinity ;

	_PathAssignment = N._PathAssignment ;
//	_PathAssignment.clear() ;

/*		BucketElimination::Bucket *b = N.Bucket()->ChildBucket(ChildIdx) ;
	int32_t path_len_ = _WS->ComputeDistanceToNextDescendantBranchingVariable(b->V()) ;
	int32_t path_len = 1 + (path_len_ < 0 ? -path_len_ : path_len_) ;
	if (_PathAssignment.capacity() < path_len) 
		_PathAssignment.reserve(path_len) ;*/

	_BranchingVarChildrenSubTreeValues.clear() ;
	_Qpath = N.QPath() ;
	double qp = N.ComputeQPathLocal(ChildIdx) ;
	_WS->ApplyFnCombinationOperator(_Qpath, qp) ;
// DEBUGGGGG
//_Qpath = 0.0 ;
	return 0 ;
}
