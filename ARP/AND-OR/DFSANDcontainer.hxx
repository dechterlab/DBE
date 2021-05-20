#ifndef DFSANDcontainer_HXX_INCLUDED
#define DFSANDcontainer_HXX_INCLUDED

#include <float.h>
#include <stdint.h>
#include <inttypes.h>

#include <MBEworkspace.hxx>
#include <SearchSpaceNodes.hxx>

namespace AndOrSearchSpace
{

class AbsSamplingWorkspace ;
class SearchAndNode ;

inline int32_t SortAndNodesWrtClosestBrVar(std::vector<AndOrSearchSpace::SearchAndNode_WithPath*> & ANDnodes, std::vector<int64_t> & sorthelper)
{
	// sorting helper
	int32_t left[32], right[32] ;

	sorthelper.clear() ;
	if (sorthelper.capacity() < ANDnodes.size()) {
		sorthelper.reserve(ANDnodes.size()) ;
		if (sorthelper.capacity() < ANDnodes.size()) 
			return 1 ;
		}
	sorthelper.clear() ;
	for (AndOrSearchSpace::SearchAndNode_WithPath *N : ANDnodes) {
		AndOrSearchSpace::SearchAndNode *n = N->ClosestBranchingAncestor() ;
		sorthelper.push_back((int64_t) n) ;
		}
	if (sorthelper.size() > 1) 
		QuickSorti64_i64(sorthelper.data(), sorthelper.size(), (int64_t*) ANDnodes.data(), left, right) ;

	return 0 ;
}


class DFSANDcontainer
{
protected :
	AbsSamplingWorkspace *_WS ; // BucketElimination::MBEworkspace *_WS ;
	int _Depth ; // depth of DFS
	DFSANDcontainer *_Parent ; // parent in the DFS tree
	DFSANDcontainer *_Child ; // child node; here so that we can track root->descendats
	int32_t _V ; // variable associated with the node
	BucketElimination::Bucket *_B ; // bucket of the variable
#ifdef _DEBUG
public :
	int32_t _nExpansions ;
#endif
public :
	inline AbsSamplingWorkspace *WS(void) const { return _WS ; }
	inline int32_t V(void) const { return _V ; }
	inline BucketElimination::Bucket *B(void) const { return _B ; }
	inline DFSANDcontainer *Parent(void) const { return _Parent ; }
	inline DFSANDcontainer *Child(void) const { return _Child ; }
	inline void AttachChild(DFSANDcontainer *C) { _Child = C ; }
protected :
	ARE_Function_TableType _SubTreeValue ; // value of this DFS node; as sum of values of all AND nodes in this container
public :
	inline ARE_Function_TableType & SubTreeValue(void) { return _SubTreeValue ; }
protected :
	std::vector<BucketElimination::Bucket*> _Children ; // children of '_V' in the bucket tree
protected :
	std::vector<SearchAndNode_WithPath*> _OpenList ; // a list of AND nodes; their values should be OR'ed together
	int32_t _idxCurrentANDnode ; // idx into '_OpenList' array; when doing non-proper DFS, this is not used, since we expand all OpenList[] at once.
	int32_t _idxCurrentChildnode ; // idx into '_Children' array.
public :
	inline int32_t nChildren(void) const { return _Children.size() ; }
	inline int32_t OpenListSize(void) const { return _OpenList.size() ; }
	inline std::vector<SearchAndNode_WithPath*> & OpenList(void) { return _OpenList ; }
	inline SearchAndNode_WithPath *OpenListNode(int32_t idx) const { return idx < 0 || idx >= _OpenList.size() ? NULL : _OpenList[idx] ; }
	inline SearchAndNode_WithPath *CurrentANDnode(void) const { return _idxCurrentANDnode < 0 || _idxCurrentANDnode >= _OpenList.size() ? NULL : _OpenList[_idxCurrentANDnode] ; }
	inline int32_t idxCurrentChildnode(void) const { return _idxCurrentChildnode ; }
	inline BucketElimination::Bucket *CurrentChildnode(void) const { return _idxCurrentChildnode < 0 || _idxCurrentChildnode >= _Children.size() ? NULL : _Children[_idxCurrentChildnode] ; }
	void SerializeCurrentState(std::string & s) 
	{
		char S[64] ;
		sprintf(S, "{V=%d;nC=%d;nOL=%d;idxC=%d;idxOL=%d}",(int)_V,(int)_Children.size(),(int)_OpenList.size(),_idxCurrentChildnode,_idxCurrentANDnode) ;
		s += S ;
	}
	inline int32_t ImportOpenList(std::vector<SearchAndNode_WithPath*> & OL)
	{
		_OpenList.clear() ;
		if (_OpenList.capacity() < OL.size()) {
			_OpenList.reserve(OL.size()) ;
			if (_OpenList.capacity() < OL.size()) 
				return 1 ;
			}
		for (SearchAndNode_WithPath *n : OL) {
			_OpenList.push_back(n) ;
			if (n->SubTreeValue() >= ARP_DBL_MAX) 
				n->SubTreeValue() = _WS->FnCombinationNeutralValue() ;
			}
		return 0 ;
	}
	inline void MarkCurrentANDNodeProcessed(void)
	{
		_idxCurrentChildnode = _Children.size() ;
	}
	inline BucketElimination::Bucket *FetchNextChildVar2Process_Proper(bool & FreshANDNode)
	{
		FreshANDNode = false ;
		SearchAndNode_WithPath *current_an = _idxCurrentChildnode < 0 ? NULL : CurrentANDnode() ;
		if (NULL != current_an) {
			if (ARP_nInfinity == current_an->SubTreeValue()) {
				FreshANDNode = true ;
				// value of current AND node is -infinity (or 0 on normal scale); there is no point to compute other children of the current AND node
				if (++_idxCurrentANDnode >= _OpenList.size()) 
					return NULL ;
				_idxCurrentChildnode = 0 ;
				SearchAndNode_WithPath *new_an = CurrentANDnode() ;
				new_an->ApplyPathAssignment() ;
				return _Children[_idxCurrentChildnode] ;
				}
			}
		if (++_idxCurrentChildnode >= _Children.size()) {
			FreshANDNode = true ;
			if (++_idxCurrentANDnode >= _OpenList.size()) 
				return NULL ;
			_idxCurrentChildnode = 0 ;
			}
		SearchAndNode_WithPath *new_an = CurrentANDnode() ;
		if (current_an != new_an) 
			new_an->ApplyPathAssignment() ;
		return _Children[_idxCurrentChildnode] ;
	}
/*	inline int32_t FetchNextChildVar2Process(void) 
	{
		if (_nChildrenProcessed >= _Children.size()) 
			return -1 ;
		return _Children[_nChildrenProcessed++]->V() ;
	}
	inline BucketElimination::Bucket *FetchNextChildBucket2Process(void) 
	{
		if (_nChildrenProcessed >= _Children.size()) 
			return NULL ;
		return _Children[_nChildrenProcessed++] ;
	}
	inline SearchOrNode *FetchOrNode2Process(void) 
	{
		if (_nChildrenProcessed >= _Children.size()) 
			return NULL ;
		SearchOrNode *on = new SearchOrNode(_WS, _Children[_nChildrenProcessed++]->V()) ;
		return on ;
	}*/

	/*
	Non-proper DFS functions.
	*/
public :
	// given a DFS branching var node, some of _Children[] have been processed; get next.
	inline BucketElimination::Bucket *FetchNextChildVar2Process_nonProper(void)
	{
		if (++_idxCurrentChildnode >= _Children.size()) 
			return NULL ;
		return _Children[_idxCurrentChildnode] ;
	}
	// import a set of AND nodes as the OpenList[] of this DFS node; upon exit, 'Frontier' will be empty - its contents will now belong to this DFSnode.
	// we assume the Qpath of each Frontier[] node is up-to-date - i.e. 
	inline int32_t ImportOpenListEx(std::vector<SearchAndNode_WithPath*> & Frontier, std::vector<int64_t> *sorthelper)
	{
		_OpenList.clear() ;
		if (_OpenList.capacity() < Frontier.size()) {
			_OpenList.reserve(Frontier.size()) ;
			if (_OpenList.capacity() < Frontier.size()) 
				return 1 ;
			}
		// copy over all nodes; clear Frontier.
		_OpenList.resize(Frontier.size(), NULL) ;
		for (int32_t i = 0 ; i < Frontier.size() ; ++i) 
			_OpenList[i] = Frontier[i] ;
		Frontier.clear() ;
		// sort openlist wrt ClosestBranchingVars of each AND node in openlist; when this DFSnode is processed, this will help pass openlist values up the DFS tree.
		// note that we check if the _OpenList[0] has ClosestBranchingAncestor; if it does not, we assume none of the nodes in _OpenList[] do.
		if (_OpenList.size() > 1 ? NULL != _OpenList[0]->ClosestBranchingAncestor() : false) {
			if (NULL != sorthelper) {
				if (0 != SortAndNodesWrtClosestBrVar(_OpenList, *sorthelper)) 
					return 1 ;
				}
			else {
				std::vector<int64_t> sorthelper_ ;
				if (0 != SortAndNodesWrtClosestBrVar(_OpenList, sorthelper_)) 
					return 1 ;
				}
			}
		// compute 'BranchingVarChildrenSubTreeValues' for each AND node; i.e. the h it receives from each child.
		// Note we don't have to compute Qpath since Qpath of each AND node in 'Frontier' should have been up-to-date.
		if (NULL != _B ? _B->nChildren() > 1 : false) {
			for (SearchAndNode_WithPath *n : _OpenList) {
				if (0 != n->InitializeBranchingVarChildrenSubTreeValues()) 
					return 1 ;
				}
			}
		return 0 ;
	}
	// export current OpenList[] as a copy; set Qpath of each node; set ClosestBranchingAncestorVarChild and its Idx.
	// note OpenList[] will stay as-is. Qpath of the new nodes is OpenList::Qpath * QpathLocal - this is ok, since 
	// we will do DFS down the pseudo(bucket)-tree path of 'ChildIdx'.
	// NOTE : OpenList[] each AND node does not include its _BranchingVarChildrenSubTreeValues.
	//		but exported nodes do; this is ok because when DFS is below the AND node, _BranchingVarChildrenSubTreeValues of the AND variables do not change.
	int32_t ExportOpenListEx(int ChildIdx, std::vector<SearchAndNode_WithPath*> & Frontier)
	{
		BucketElimination::Bucket *b = _Children[ChildIdx] ;
		Frontier.clear() ;
		if (Frontier.capacity() < _OpenList.size()) {
			Frontier.reserve(_OpenList.size()) ;
			if (Frontier.capacity() < _OpenList.size()) 
				return 1 ;
			}
		for (SearchAndNode_WithPath *N : _OpenList) {
			SearchAndNode_WithPath *n = new SearchAndNode_WithPath(*N) ;
			if (NULL == n) return 1 ;
			n->PrepForDFSExpansion(*N, ChildIdx) ;
			Frontier.push_back(n) ;
			}
		return 0 ;
	}
	// assuming a list of AND nodes have been processed - i.e. their estimates (Z) have been computed, 
	// propagate values up the DFS tree.
	// arguments :
	//		1) Frontier - either _Openlist[] of this DFS node, or a list of leaf nodes of the AND/OR search tree that are direct descendants of this DFS node.
	//		2) NeedSort - if AndList==_Openlist[] don't need sort, since it is already sorted, otherwise need sort.
	// return :
	//	   -1 = Frontier is empty; entire Frontier can be eliminated if all are 0 (normal-space); this means root node value is also 0.
	//      0 = ok
	//     +1 = failed
	inline int NoteAndListCompletion(std::vector<SearchAndNode_WithPath*> & Frontier, bool NeedSort, std::vector<int64_t> *sorthelper, CMauiAVL64Tree & node_done_helper, bool IsLeafNode)
	{
		if (0 == Frontier.size()) 
			return -1 ;
		SearchAndNode_WithPath *ccba = (SearchAndNode_WithPath *) Frontier[0]->ClosestBranchingAncestor() ;
		if (NULL == ccba) 
			return 0 ; // Frontier should be just 1 node, which is the root AND node of the AND/OR search tree

		// if needed/requested, sort Frontier wrt ClosestBranchingVars of each AND node in Frontier.
		if (NeedSort && Frontier.size() > 1) {
			if (NULL != sorthelper) 
				SortAndNodesWrtClosestBrVar(Frontier, *sorthelper) ;
			else 
				{ std::vector<int64_t> sorthelper_ ; SortAndNodesWrtClosestBrVar(Frontier, sorthelper_) ; }
			}

		// By construction, nodes of _OpenList are the ClosestBranchingAncestors of Frontier.
		// We need to keep track which nodes in _OpenList receive a msg from Frontier; those that don't need to be removed (discarded).
		// collect all current_dfs nodes in 'node_done_helper'; access is lon(N) time; for very small sets this is overkill, but it should be more scalable.
		node_done_helper.Empty() ;
		for (int32_t i = 0 ; i < _OpenList.size() ; ++i) 
			node_done_helper.Insert((int64_t) _OpenList[i], i) ;

		ccba = (SearchAndNode_WithPath *) Frontier[0]->ClosestBranchingAncestor() ; // current closest branching ancestor
		int32_t j = 0 ;
		for (int32_t i = 1 ; i <= Frontier.size() ; ++i) {
			if (i < Frontier.size() ? ccba != Frontier[i]->ClosestBranchingAncestor() : true) {
				// nodes [j,i) is a list of AND nodes belonging to ccba. 

				// compute sum (normal-space) of nodes [j,i); later we will combine the sum with ccba.
				double V = ARP_nInfinity ;
				for (int32_t k = j ; k < i ; ++k) {
					SearchAndNode_WithPath *n = Frontier[k] ;
					ARE_Function_TableType v = _WS->FnCombinationNeutralValue() ;
					if (IsLeafNode) {
						_WS->ApplyFnCombinationOperator(v, n->AccumulatedISwSinceLastBranchingVariable()) ;
						_WS->ApplyFnCombinationOperator(v, n->AccumulatedCostSinceLastBranchingVariable()) ;
						_WS->ApplyFnCombinationOperator(v, n->h()) ;
						}
					else {
						v = n->SubTreeValue() ;
						if (ARP_DBL_MAX == v) 
							continue ; // this (internal) AND node did not get a value, probably because it the sampled AND/OR tree it has no children; skip it.
						_WS->ApplyFnCombinationOperator(v, n->AccumulatedISwSinceLastBranchingVariable()) ;
						_WS->ApplyFnCombinationOperator(v, n->AccumulatedCostSinceLastBranchingVariable()) ;
						}
					LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, v) ;
					}
				// send V to closest ancestor; if V=-infinity process it still - ancestor value will become -infinity, but that is fine.
				if (ARP_DBL_MAX == ccba->SubTreeValue()) // first time value is set
					ccba->SubTreeValue() = V ;
				else // combine with previous value
					_WS->ApplyFnCombinationOperator(ccba->SubTreeValue(), V) ;

				// if ccba became 0, then it should be removed; ignore it for now
				if (ccba->SubTreeValue() == ARP_nInfinity) {
					}
				else {
					// in ccba, replace Z^ of this child (which originally is computed from h of the subtree) with this newly computed value.
					ccba->UpdateBranchingVarChildrenSubTreeValue(Frontier[j]->ClosestBranchingAncestorVarChildIdx(), V) ;

					// remove ccba from 'node_done_helper'
					node_done_helper.Remove((int64_t) ccba, NULL) ;
					}

				if (i >= Frontier.size()) 
					break ;
				j = i ;
				ccba = (SearchAndNode_WithPath *) Frontier[j]->ClosestBranchingAncestor() ;
				}
			}

		// dump all _OpenList nodes that did not get updated
		if (node_done_helper.GetSize() > 0) {
			if (node_done_helper.GetSize() == _OpenList.size()) {
				return -1 ; // entire OpenList is all 0's
				}
			// first just delete _OpenList[] nodes that did not get updated; later compact the _OpenList[] array so that order stays the same
			long pos = -1 ; int64_t ptr, idx ;
			while (node_done_helper.GetNext(&ptr, &idx, pos) > 0) {
				SearchAndNode_WithPath *n = (SearchAndNode_WithPath *) ptr ;
				_OpenList[idx] = NULL ;
				delete n ;
				}
			int32_t j = 0 ;
			for (int32_t i = 0 ; i < _OpenList.size() ; ++i) {
				if (NULL != _OpenList[i]) 
					_OpenList[j++] = _OpenList[i] ;
				}
			_OpenList.resize(j) ;
			}
		if (0 == _OpenList.size()) 
			return -1 ; // this should not happen

		return 0 ;
	}

public :
	inline int32_t InitializeOpenListSubTreeValue(void)
	{
		for (SearchAndNode_WithPath *n : _OpenList) {
			n->SubTreeValue() = n->AccumulatedISwSinceLastBranchingVariable() ;
			_WS->ApplyFnCombinationOperator(n->SubTreeValue(), n->AccumulatedCostSinceLastBranchingVariable()) ;
			}
		return 0 ;
	}
	int32_t Initialize(AbsSamplingWorkspace*WS, DFSANDcontainer *Parent, int32_t Var)
	{
		_WS = WS ;
		_Depth = NULL != Parent ? 1+Parent->_Depth : 0 ;
		if (_Parent != Parent) {
			if (NULL != _Parent) {
				_Parent->AttachChild(NULL) ;
				}
			_Parent = Parent ;
			if (NULL != _Parent) {
				if (NULL != _Parent->Child()) 
					delete _Parent->Child() ;
				_Parent->AttachChild(this) ;
				}
			}
		_Child = NULL ;
		_V = Var ;
		_B = _V >= 0 ? WS->MapVar2Bucket(_V) : NULL ;
		_SubTreeValue = _WS->FnCombinationNeutralValue() ;
		_Children.clear() ;
		if (NULL != _WS) {
			_WS->GetPseudotreeChildrenOfVariable(_V, _Children) ;
			}
#ifdef _DEBUG
		_nExpansions = 0 ;
#endif
		return 0 ;
	}
	DFSANDcontainer(AbsSamplingWorkspace *WS = NULL, DFSANDcontainer *Parent = NULL, int32_t Var = -1) :
		_WS(WS), 
		_Depth(-1), 
		_Parent(NULL), 
		_Child(NULL), 
		_V(Var), 
		_B(NULL), 
		_idxCurrentANDnode(0), 
		_idxCurrentChildnode(-1)
	{
		Initialize(WS, Parent, Var) ;
	}
	virtual ~DFSANDcontainer(void)
	{
		for (SearchAndNode_WithPath *n : _OpenList) 
			delete n ;
		_OpenList.clear() ;
	}
} ;

} // AndOrSearchSpace

#endif // DFSANDcontainer_HXX_INCLUDED
