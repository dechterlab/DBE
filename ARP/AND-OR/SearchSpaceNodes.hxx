#ifndef SearchSpaceNodes_HXX_INCLUDED
#define SearchSpaceNodes_HXX_INCLUDED

#include <inttypes.h>

#include "MBEworkspace.hxx"

namespace AndOrSearchSpace
{

class AbsSamplingWorkspace ;

class SearchNode
{
protected :
	AbsSamplingWorkspace *_WS ; // BucketElimination::MBEworkspace *_WS ;
	int32_t _V ; // variable associated with the node
public :
	inline AbsSamplingWorkspace *WS(void) const { return _WS ; }
	inline int32_t V(void) const { return _V ; }
	BucketElimination::Bucket *Bucket(void) const ;

protected :
	ARE_Function_TableType _SubTreeValue ;
public :
	inline ARE_Function_TableType & SubTreeValue(void) { return _SubTreeValue ; }

protected :
	SearchNode *_Parent ;
	int32_t _nChildren ;
	SearchNode *_Children ; // head of a list; each child points to its sibling, which point to next sibling, etc.
	SearchNode *_ChildrenEndOfList ; // tail of a list.
public :
	inline SearchNode *Parent(void) const { return _Parent ; }
	inline int32_t nChildren(void) const { return _nChildren ; }
	inline SearchNode *Children(void) const { return _Children ; }
	inline void AttachChildren(SearchNode *N, int32_t n) { _Children = N ; _nChildren = n ; if (NULL == N) _ChildrenEndOfList = NULL ; }
	inline void AddChild(SearchNode & N)
	{
		if (NULL == _Children) 
			{ _Children = _ChildrenEndOfList = &N ; _nChildren = 1 ; }
		else 
			{ _ChildrenEndOfList->AttachSibling(&N) ; _ChildrenEndOfList = &N ; ++_nChildren ; }
		N.AttachSibling(NULL) ;
	}
	inline void RemoveChild(SearchNode & N)
	{
		SearchNode *nPrevious = NULL ;
		for (SearchNode *n = _Children ; NULL != n ; n = n->Sibling()) {
			if (&N == n) {
				if (NULL != nPrevious){ 
					nPrevious->AttachSibling(n->Sibling()) ;
					if (NULL == n->Sibling()) _ChildrenEndOfList = nPrevious; // Changed Filjor - 170629
				}
				else { 
					_Children = n->Sibling() ;
				}
				n->AttachSibling(NULL) ;
				--_nChildren ;
				if (NULL == _Children) 
					_ChildrenEndOfList = NULL ;
				return ;
				}
			nPrevious = n ;
			}
	}
	int32_t SetParent(SearchNode *P, int32_t listsize, bool UpdateParent)
	{
		SearchNode *old_parent = _Parent ;
		_Parent = P ;
		if (UpdateParent) {
			// remove this from old parent
			if (NULL != old_parent) {
				SearchNode *nPrevious = NULL ;
				for (SearchNode *n = old_parent->Children() ; NULL != n ; n = n->Sibling()) {
					if (this == n) {
						if (NULL == nPrevious) {
							old_parent->AttachChildren(Sibling(), old_parent->nChildren()-1) ;
							}
						else {
							nPrevious->AttachSibling(Sibling()) ;
							old_parent->AttachChildren(Children(), old_parent->nChildren()-1) ;
							}
						AttachSibling(NULL) ;
						break ;
						}
					nPrevious = n ;
					}
				}
			// add to new parent
			if (NULL != _Parent) {
				if (listsize < 0) 
					{ listsize = 0 ; for (SearchNode *n = this ; NULL != n ; n = n->Sibling(), ++listsize) ; }
				if (NULL != _Parent->Children()) {
					SearchNode *last = NULL ;
					for (SearchNode *n = _Parent->Children() ; NULL != n ; n = n->Sibling()) last = n ;
					last->AttachSibling(this) ;
					_Parent->AttachChildren(Children(), _Parent->nChildren() + listsize) ;
					}
				else 
					_Parent->AttachChildren(this, listsize) ;
				}
			}
		return 0 ;
	}

protected :
	SearchNode *_Sibling ; // used to create a list of children of a parent
public :
	inline SearchNode *Sibling(void) const { return _Sibling ; }
	inline void AttachSibling(SearchNode *S) { _Sibling = S ; }

protected :
	SearchNode *_NextInExpansionQueue ; // used to create an expansion queue of the AND/OR workspace as a list.
public :
	inline SearchNode * & NextInExpansionQueue(void) { return _NextInExpansionQueue ; }

public :
	SearchNode(AbsSamplingWorkspace *WS, int32_t V) :
		_WS(WS), 
		_V(V), 
		_SubTreeValue(ARP_DBL_MAX), 
		_Parent(NULL), 
		_nChildren(0), _Children(NULL), _ChildrenEndOfList(NULL), _Sibling(NULL), _NextInExpansionQueue(NULL)
	{
	}
	SearchNode(SearchNode & N) :
		_WS(N._WS), 
		_V(N._V), 
		_SubTreeValue(N._SubTreeValue), 
		_Parent(N._Parent), 
		_nChildren(0), _Children(NULL), _ChildrenEndOfList(NULL), _Sibling(NULL), _NextInExpansionQueue(NULL)
	{
	}
	virtual ~SearchNode(void)
	{
	}
} ;

class SearchOrNode : public SearchNode
{
public :
	SearchOrNode(AbsSamplingWorkspace *WS, int32_t Var) : SearchNode(WS, Var)
	{
	}
	virtual ~SearchOrNode(void)
	{
	}
} ;

class SearchAndNode : public SearchNode
{
protected :
	int32_t _Value ; // value of the variable associated with the node
public :
	inline int32_t Value(void) const { return _Value ; }

protected :
	int32_t _VarPosInOrder ; // idx of the variable in the bucket tree order; 0 means root of the tree.
public :
	inline int32_t VarPosInOrder(void) const { return _VarPosInOrder ; }
	inline bool IsBranchingVariable(void) const
	{
		BucketElimination::Bucket *b = ((BucketElimination::MBEworkspace*)_WS)->MapVar2Bucket(_V) ;
		if (NULL == b) 
			// this means that this node is the (AND) root of the AND/OR search tree; these are by definition branching variables.
			return true ;
		return b->nChildren() > 1 ;
	}

protected :
	// closest branching ancestor is a closest ancestor AND node in AND/OR search tree of a variable that has >1 children in the pseudo(bucket) tree.
	SearchAndNode *_ClosestBranchingAncestor ;
	// idx of the ClosestBranchingAncestor variable wrt its parent var; this allows us to quickly find this var (or a path to this var) as a child of the ClosestBranchingAncestor.
	int32_t _ClosestBranchingAncestorVarChildIdx ;
public :
	inline SearchAndNode * & ClosestBranchingAncestor(void) { return _ClosestBranchingAncestor ; }
	inline int32_t & ClosestBranchingAncestorVarChildIdx(void) { return _ClosestBranchingAncestorVarChildIdx ; }
	inline void SetClosestBranchingAncestor(SearchAndNode *ClosestAncestorAndNode, int32_t idx)
	{
		if (NULL == ClosestAncestorAndNode) {
			_ClosestBranchingAncestor = NULL ;
			_ClosestBranchingAncestorVarChildIdx = -1 ;
			}
		else if (ClosestAncestorAndNode->IsBranchingVariable()) {
			_ClosestBranchingAncestor = ClosestAncestorAndNode ;
			_ClosestBranchingAncestorVarChildIdx = idx ;
			}
		else {
			_ClosestBranchingAncestor = ClosestAncestorAndNode->ClosestBranchingAncestor() ;
			_ClosestBranchingAncestorVarChildIdx = ClosestAncestorAndNode->ClosestBranchingAncestorVarChildIdx() ;
			}
	}

protected :
	int64_t _ExpansionQueueKey ; // idx (key) of this node in the queue.
public :
	inline int64_t & ExpansionQueueKey(void) { return _ExpansionQueueKey ; }

protected :
	// h value of this node
	ARE_Function_TableType _h ;
public :
	inline ARE_Function_TableType & h(void) { return _h ; }

protected :
	// combined value of all Original functions in the bucket
	ARE_Function_TableType _Cost ;
	// combined total of all AND node labels from root to here
	ARE_Function_TableType _AccumulatedCost ;
	// combined total of all AND node labels from last branching variable to here
	ARE_Function_TableType _AccumulatedCostSinceLastBranchingVariable ;
public :
	inline ARE_Function_TableType & Cost(void) { return _Cost ; }
	inline ARE_Function_TableType & AccumulatedCost(void) { return _AccumulatedCost ; }
	inline ARE_Function_TableType & AccumulatedCostSinceLastBranchingVariable(void) { return _AccumulatedCostSinceLastBranchingVariable ; }

protected :
	// local IS weight
	double _ISw ;
	// combined total of all AND node labels from root to here
	ARE_Function_TableType _AccumulatedISw ;
	// combined total of all AND node labels from last branching variable to here
	ARE_Function_TableType _AccumulatedISwSinceLastBranchingVariable ;
	// IS proposal q
	ARE_Function_TableType _ISq ;
	// number of times this node has been merged (plus 1).
	int32_t _nISmerges ;
public :
	inline double & ISw(void) { return _ISw ; }
	inline double & AccumulatedISw(void) { return _AccumulatedISw ; }
	inline double & AccumulatedISwSinceLastBranchingVariable(void) { return _AccumulatedISwSinceLastBranchingVariable ; }
	inline double & ISq(void) { return _ISq ; }
	inline int32_t & nISmerges(void) { return _nISmerges ; }

public :
	SearchAndNode(void) ;
	SearchAndNode(AbsSamplingWorkspace *WS, int32_t Var, int32_t Value) ;
	SearchAndNode(SearchAndNode & N) ;
	virtual ~SearchAndNode(void)
	{
	}
} ;

//extern int64_t nNODES__ ;

class SearchAndNode_WithPath : public SearchAndNode
{
public :
	// list of var=value assignments, from earliest[0] to latest[size-1]; last value is the assignment to the variable of this AND node. 
	// note that here is the full path assignment, not just context assignment.
	// note that the variables included are (closest branching variable of this AND node, this AND node].
	std::vector<int32_t> _PathAssignment ;
public :
	inline std::vector<int32_t> & PathAssignment(void) { return _PathAssignment ; }
	inline int32_t PathAssignmentValue(int32_t idx) const { return _PathAssignment[idx] ; }
	void ApplyPathAssignment(bool SetFullPathAssignment = false) ;
public :
	// for branching variables, we store here SubTree (i.e. local) values for each child.
	// initially they are filled in with (product of) h-values that this node receives from each child.
	// this var is used when doing non-proper abs sampling (DFS); not used when doing proper abs sampling.
	std::vector<double> _BranchingVarChildrenSubTreeValues ;
public :
	// initialize _BranchingVarChildrenSubTreeValues[] with h values of the chilren. size of _BranchingVarChildrenSubTreeValues[] is Bucket()->nChildren().
	int32_t InitializeBranchingVarChildrenSubTreeValues(void) ;
	// update BranchingVarChildrenSubTreeValue
	int32_t UpdateBranchingVarChildrenSubTreeValue(int32_t idx, double v)
	{
		if (idx < 0 || idx >= _BranchingVarChildrenSubTreeValues.size()) 
			return 1 ;
		_BranchingVarChildrenSubTreeValues[idx] = v ; // 
		return 0 ;
	}
public :
	// let : V(v,j) is the product of subtree values (estimates) of children of variable 'v', except for child 'j'.
	// let : P(this) is the path (this, root] in the AND/OR search tree.
	// let : BV(this) is the set of branching variables on P(this).
	// let : Qpath is the product of V(w,j) where w \in BV(this) and j is a child of w that in in P(this)
	double _Qpath ;
public :
	inline double & QPath(void) { return _Qpath ; }
	// compute the product of all BranchingVarChildrenSubTreeValues, except for [idxToExclude]
	double ComputeQPathLocal(int idxToExclude) ;
public :
	// prepare this node for DFS expansion from N
	int32_t PrepForDFSExpansion(SearchAndNode_WithPath & N, int32_t ChildIdx) ;
	SearchAndNode_WithPath(void) : SearchAndNode(), _Qpath(ARP_DBL_MAX) { /* ++nNODES__ ; */ }
	SearchAndNode_WithPath(AbsSamplingWorkspace *WS, int32_t Var, int32_t Value) ;
	SearchAndNode_WithPath(SearchAndNode_WithPath & N) ;
//	~SearchAndNode_WithPath(void) { --nNODES__ ; }
} ;

} // AndOrSearchSpace

#endif // SearchSpaceNodes_HXX_INCLUDED
