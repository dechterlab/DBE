#include <vector>
#include <cmath>

#include "AVLtreeObjects.hxx"
#include "Globals.hxx"
#include "MBEworkspace.hxx"
#include "SearchSpaceNodes.hxx"
#include <AbsSamplingWorkspace.hxx>
#include <DFSANDcontainer.hxx>

//#define DEBUG_DFS_TRAVERSAL
//#define IGNORE_PROPER_ABS_NOSORTING
//#define IGNORE_CONTEXT_ABS_SPECIALCODE
//#define USE_WMBq

static void DeleteDFSstack(AndOrSearchSpace::DFSANDcontainer *dfs_current)
{
	// leave root node intact; delete everything else.
	AndOrSearchSpace::DFSANDcontainer *p ;
	for (; NULL != dfs_current->Parent() ; dfs_current = p) {
		p = dfs_current->Parent() ;
		delete dfs_current ;
		}
}

int32_t AndOrSearchSpace::AbsSamplingWorkspace::RunDFS_nonProper(int32_t & dDeepestBucketExpanded, int32_t & nDFSBranchingPointsProcessed)
{
	dDeepestBucketExpanded = nDFSBranchingPointsProcessed = 0 ;

#ifdef DEBUG_DFS_TRAVERSAL
FILE *fpLOG = fopen("C:\\UCI\\ARP\\AND-OR\\log.txt", "w") ;
#endif // DEBUG_DFS_TRAVERSAL
	_Root.SubTreeValue() = _AnswerFactor ;
	_nNodesInTree = _nNodesCreated = _nNodesISmerged = _nExpansionCalls = _nTotalSumFrontierSize = _nNumVarExpansions = 0 ;

	for (int32_t i = 0 ; i < _Problem->N() ; ++i) 
		_CurrentContextValues[i] = -1 ;

	bool no_AND_sorting = false, merge_all = false ;
	if (AbsSamplingTwoAndNodeCompare_Unique == _CompFn) 
		no_AND_sorting = true ;
	else if (AbsSamplingTwoAndNodeCompare_Knuth == _CompFn) 
		no_AND_sorting = merge_all = true ;
	else if (AbsSamplingTwoAndNodeCompare_CustomProper_DFS == _CompFn ? MaxContextSize() <= 1 : false) 
		no_AND_sorting = merge_all = true ;
	// for context-based abstraction, the CompFn should be AbsSamplingTwoAndNodeCompare_CustomProper_DFS()

	// temp array for sort keys
	std::vector<int32_t> sort_keys ;
	if (AbsSamplingTwoAndNodeCompare_RandCntxt == _CompFn) 
		sort_keys.reserve(100) ;

	// temporary arrays of AND openlist nodes
	std::vector<SearchAndNode_WithPath*> OL_from, OL_to ;
	OL_from.reserve(100) ; OL_to.reserve(100) ; // allocate some space

	CMauiAVL64Tree node_done_helper(0, NULL, NULL, N(), 1024) ;

	// sorting helper
	int32_t left[32], right[32] ;
	std::vector<int64_t> sorthelper ;

/*	std::vector<BucketElimination::Bucket*> roots ;
	int32_t resRoots = GetRoots(roots) ;

	for (BucketElimination::Bucket *bRoot : roots) {
		int32_t vRoot = bRoot->V() ;
		// check if h is exact; process completely here if yes; this should not happen, except at the root
		if (EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(vRoot) : false) {
			_nNodesExactHComputed += 1 ;
			// get the output value of the bucket
			double v = FnCombinationNeutralValue() ;
			std::vector<BucketElimination::MiniBucket*> & MBs = bRoot->MiniBuckets() ;
			for (BucketElimination::MiniBucket *mb : MBs) {
				ARE::Function & f = mb->OutputFunction() ;
				ApplyFnCombinationOperator(v, f.ConstValue()) ;
				}
			// add v to _Root value
			ApplyFnCombinationOperator(_Root.SubTreeValue(), v) ;
			continue ;
			}

		// push this root node
		DFSANDcontainer *dfs_root = new DFSANDcontainer(this, NULL, vRoot) ;
		if (NULL == dfs_root) return 1 ;
		int32_t idxS = 0, idxE = _Problem->K(vRoot) ;
		if (_Problem->Value(vRoot) >= 0) 
			{ idxS = _Problem->Value(vRoot) ; idxE = idxS+1 ; }
		for (int32_t iV = idxS ; iV < idxE ; ++iV) {
			_CurrentContextValues[vRoot] = iV ;
			SearchAndNode_WithPath *anNew = new SearchAndNode_WithPath(this, vRoot, iV) ;
			if (NULL == anNew) 
				return 1 ;
			++_nNodesCreated ;
			bRoot->ComputeCost(_CurrentContextValues, anNew->Cost()) ;
			bRoot->ComputeHeuristic(_CurrentContextValues, anNew->h()) ;
			anNew->AccumulatedCost() = anNew->Cost() ;
			anNew->AccumulatedCostSinceLastBranchingVariable() = anNew->Cost() ;
			anNew->ISq() = anNew->AccumulatedCost() ;
			ApplyFnCombinationOperator(anNew->ISq(), anNew->h()) ;
			if (ARP_nInfinity == anNew->ISq()) 
				{ delete anNew ; continue ; } // ignore this node; ISq=-inf means c or h is -inf; h is upper boound.
			OL_from.push_back(anNew) ;
			}
		dfs_root->ImportOpenListEx(OL_from, &sorthelper) ; OL_from.clear() ;
// bug	dfs_root->InitializeOpenListSubTreeValue() ;
*/

	DFSANDcontainer *dfs_root = new DFSANDcontainer(this, NULL, -1) ;
	if (NULL == dfs_root) return 1 ;
	SearchAndNode_WithPath *and_root = new SearchAndNode_WithPath(this, -1, -1) ; // dummy root AND node; product of all root buckets.
	if (NULL == and_root) return 1 ;
	++_nNodesCreated ;
	OL_from.push_back(and_root) ;
	dfs_root->ImportOpenListEx(OL_from, &sorthelper) ; OL_from.clear() ;

	// at top level, we expand a frontier of AND-nodes from a branching var to a next descendant branching var.
	DFSANDcontainer *dfs_current = dfs_root ; // current DFS node
	while (NULL != dfs_current) {
		DFSANDcontainer *p = dfs_current->Parent() ;

		// at this level, we pick pick AND-node frontier from the btm node, and expand it in the direction of next branching variable child.
		BucketElimination::Bucket *BC = dfs_current->FetchNextChildVar2Process_nonProper() ;
		// if all children at this level are exhausted, backtrack
		if (NULL == BC) {
			// in the current DFS node, parents of the AND nodes may be different; we need to OR those AND nodes that belong to the same parent.
			// assume that dfs_current->_OpenList is sorted wrt ClosestBranchingVars (this should have been done when the dfs_current was created).
			if (NULL != p) 
				p->NoteAndListCompletion(dfs_current->OpenList(), false, NULL, node_done_helper, false) ;
			// done; backtrack.
			if (dfs_current != dfs_root) delete dfs_current ;
			dfs_current = p ;
			if (NULL != dfs_current) dfs_current->AttachChild(NULL) ;
			continue ;
			}

		int32_t VC = BC->V() ; // VC/BC is a child of dfs_current variable.

		// need to make a local copy of dfs_current OpenList
		if (0 != dfs_current->ExportOpenListEx(dfs_current->idxCurrentChildnode(), OL_from)) 
			return 1 ;

		++nDFSBranchingPointsProcessed ;

		int32_t nLevels = 0 ;
		BucketElimination::Bucket *bC = BC ;

expand_more :
		++nLevels ;
		++_nNumVarExpansions ;
		_nTotalSumFrontierSize += OL_from.size() ;
		int32_t vC = bC->V() ;
		int32_t nC = bC->nChildren() ;
		BucketElimination::Bucket *bParent = bC->ParentBucket() ;
		// expand the frontier down one level from [OL_from, parent_of(vC/bC)]->[OL_to,vC/bC]
		// i.e. we are computing a new frontier corresponding to vC/bC, which is one level down.

		if (bC->DistanceToRoot() > dDeepestBucketExpanded) dDeepestBucketExpanded = bC->DistanceToRoot() ;

		if (0 == OL_from.size()) {
			// Z estimate of the problem is 0 (normal space). done here.
			and_root->SubTreeValue() = ARP_nInfinity ; // dfs_root->SubTreeValue() = ARP_nInfinity ;
			DeleteDFSstack(dfs_current) ;
			break ;
			}

		// if bC has no children, generate new frontier, compute its value and rewind (don't need to apply abstraction)...
		if (0 == nC) {
			// compute AND children of the frontier
			for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
				SearchAndNode_WithPath *an = OL_from[j] ;
				int32_t res_expand = ComputeAndChildrenEx((SearchAndNode_WithPath*)an->ClosestBranchingAncestor(), an->ClosestBranchingAncestorVarChildIdx(), an, bC, 1 == nLevels ? true : false, OL_to) ;
				}
			// delete OL_from nodes
			for (SearchAndNode_WithPath *an : OL_from) delete an ;
			OL_from.clear() ;
			// compute/propagate values
			dfs_current->NoteAndListCompletion(OL_to, true, &sorthelper, node_done_helper, true) ;
			// some cleanup
			for (SearchAndNode_WithPath *an : OL_to) delete an ;
			OL_to.clear() ;
			continue ;
			}

		// BASIC IDEA :
		//		1) expand one level down by generating AND children of current frontier (OL_from)
		//		2) check if h is exact at this level; if yes, compute exact and be done with this DFS branch
		//		3) sort the level (except when Knuth's/unique abstraction, in which case we don't need to sort since everything/nothing is merged).
		//		4) do the merging
		//		5) repeat the path down (vC/bC is non-branching variable) or create dfs node (vC/bC is a branching variable).

		// compute AND children of the frontier
		for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
			SearchAndNode_WithPath *an = OL_from[j] ;
			double newAndNodeISw = 1 == nLevels ? FnCombinationNeutralValue() : an->ISw() ;
			int32_t res_expand = ComputeAndChildrenEx((SearchAndNode_WithPath*)an->ClosestBranchingAncestor(), an->ClosestBranchingAncestorVarChildIdx(), an, bC, 1 == nLevels ? true : false, OL_to) ;
			}
		// delete OL_from nodes they are no longer needed
		for (SearchAndNode_WithPath *an : OL_from) delete an ;
		OL_from.clear() ;
		// if OL_to is empty, then that means all AND children had h/isq equal to -inf (0).
		// since entire frontier is wiped out, we are done.
		// Z estimate of the problem is 0 (normal space). done here.
		if (0 == OL_to.size()) {
			and_root->SubTreeValue() = ARP_nInfinity ; // dfs_root->SubTreeValue() = ARP_nInfinity ;
			DeleteDFSstack(dfs_current) ;
			break ;
			}

		// if bC(vC) has exact h, compute values here
		if (EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(bC->V()) : false) {
			_nNodesExactHComputed += OL_to.size() ;
			dfs_current->NoteAndListCompletion(OL_to, true, &sorthelper, node_done_helper, true) ;
			// delete OL_from nodes; at this point, OL_from contains newly generated nodes
			for (SearchAndNode_WithPath *an : OL_to) delete an ;
			OL_to.clear() ;
			continue ;
			}

		if (no_AND_sorting) {
			if (merge_all) {
				SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, 0, OL_to.size()) ;
				OL_from.push_back(an_merged) ;
				}
			// else leave all as is
			else {
				OL_from = OL_to ; OL_to.clear() ;
				}
			}
		// randomized abstractions; compute abs index for each node, sort and merge
		else if (AbsSamplingTwoAndNodeCompare_RandCntxt == _CompFn) {
			sort_keys.clear() ;
			if (sort_keys.capacity() < OL_to.size()) sort_keys.reserve(OL_to.size()) ;
			for (int32_t idx_ = 0 ; idx_< OL_to.size() ; ++idx_) {
				AndOrSearchSpace::SearchAndNode_WithPath *n = OL_to[idx_] ;
				sort_keys.push_back(ComputeNodeAbstrationID(*n)) ;
				}
			QuickSortLong_i64(sort_keys.data(), OL_to.size(), (int64_t*) OL_to.data(), left, right) ;
			// go through the output openlist and do IS merges
			int32_t idxS = 0, idxE = -1 ;
			for (idxE = idxS+1 ; idxE < OL_to.size() ; ++idxE) {
				if (sort_keys[idxS] != sort_keys[idxE]) {
					SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, idxS, idxE) ;
					idxS = idxE ;
					OL_from.push_back(an_merged) ;
					}
				}
			SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, idxS, idxE) ;
			OL_from.push_back(an_merged) ;
			}
		// context-based abstraction
		else {
			// sort output openlist
			if (OL_to.size() > 1) 
				QuickSort((void**) OL_to.data(), OL_to.size(), left, right, _CompFn) ;
			// go through the output openlist and do IS merges
			int32_t idxS = 0, idxE = -1 ;
			for (idxE = idxS+1 ; idxE < OL_to.size() ; ++idxE) {
				int32_t comp_res = (*_CompFn)(OL_to[idxS], OL_to[idxE]) ;
				if (0 != comp_res) {
					SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, idxS, idxE) ;
					idxS = idxE ;
					OL_from.push_back(an_merged) ;
					}
				}
			SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, idxS, idxE) ;
			OL_from.push_back(an_merged) ;
			}

// done with generating next-level frontier; note OL_from contains newly generated nodes
frontier_expansion_done :

		OL_to.clear() ;
		if (0 == OL_from.size()) { // this means frontier would be all 0's
			// Z estimate of the problem is 0 (normal space). done here.
			and_root->SubTreeValue() = ARP_nInfinity ; // dfs_root->SubTreeValue() = ARP_nInfinity ;
			DeleteDFSstack(dfs_current) ;
			break ;
			}

		// if branching variable, create a new DFS node
		if (nC > 1) {
			DFSANDcontainer *dfs_node = new DFSANDcontainer(this, dfs_current, vC) ;
			if (NULL == dfs_node) return 1 ;
			int32_t res = dfs_node->ImportOpenListEx(OL_from, &sorthelper) ;
			// delete OL_from nodes; at this point, OL_from contains newly generated nodes
			for (SearchAndNode_WithPath *an : OL_from) delete an ;
			OL_from.clear() ;
			if (0 != res) {
				delete dfs_node ;
				return 1 ;
				}
			// note : we will keep the subtree value of all open list node as 1; when backtracking, we will incorporate accumulated (since last branching variable) ISw & Cost.
// bug		dfs_node->InitializeOpenListSubTreeValue() ;
			dfs_current = dfs_node ;
			}
		// else continue to next level down
		else {
			int32_t cv = bC->ChildVar(0) ;
			bC = MapVar2Bucket(cv) ;
			goto expand_more ;
			}
		}

	ApplyFnCombinationOperator(_Root.SubTreeValue(), and_root->SubTreeValue()) ;
	delete dfs_root ; dfs_root = NULL ;

done :

#ifdef DEBUG_DFS_TRAVERSAL
fclose(fpLOG) ; fpLOG = NULL ;
#endif // DEBUG_DFS_TRAVERSAL

	return 0 ;
}

