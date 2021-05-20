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

int32_t AndOrSearchSpace::AbsSamplingWorkspace::RunDFS(void)
{
#ifdef DEBUG_DFS_TRAVERSAL
FILE *fpLOG = fopen("C:\\UCI\\ARP\\AND-OR\\log.txt", "w") ;
#endif // DEBUG_DFS_TRAVERSAL
	_Root.SubTreeValue() = _AnswerFactor ;
	_nNodesInTree = _nNodesCreated = _nNodesISmerged = 0 ;

	for (int32_t i = 0 ; i < _Problem->N() ; ++i) 
		_CurrentContextValues[i] = -1 ;

//	if (AbsSamplingTwoAndNodeCompare_CustomProper_DFS == _CompFn ? 0 == abs_context.size() : false) 
//		_CompFn = AbsSamplingTwoAndNodeCompare_Knuth ;
	bool no_AND_sorting = false ;
#if !defined(IGNORE_PROPER_ABS_NOSORTING)
	if (AbsSamplingTwoAndNodeCompare_Knuth == _CompFn) 
		no_AND_sorting = true ;
#endif 

	// temporary arrays of AND openlist nodes
	std::vector<SearchAndNode_WithPath*> OL_from, OL_to ;
	OL_from.reserve(100) ; OL_to.reserve(100) ; // allocate some space

	// push roots in the first node
	DFSANDcontainer *dfs_root = new DFSANDcontainer(this, NULL, -1) ;
	if (NULL == dfs_root) return 1 ;
	SearchAndNode_WithPath *an = new SearchAndNode_WithPath(this, -1, -1) ; // dummy root AND node; product of all root buckets.
	if (NULL == an) return 1 ;
	++_nNodesCreated ;
	OL_from.push_back(an) ;
	dfs_root->ImportOpenList(OL_from) ; OL_from.clear() ;
// bug	dfs_root->InitializeOpenListSubTreeValue() ;

	// sorting helper
	int32_t left[32], right[32] ;

	// when using context based abstraction, we don't store, instead do adr indexing directly
	if (AbsSamplingTwoAndNodeCompare_CustomProper_DFS == _CompFn && _max_abs_context_num_configs > 0 ? OL_to.capacity() < _max_abs_context_num_configs : false) {
		if (OL_to.capacity() < _max_abs_context_num_configs) {
			OL_to.reserve(_max_abs_context_num_configs) ;
			if (OL_to.capacity() < _max_abs_context_num_configs) 
				return 1 ;
			}
		if (OL_from.capacity() < _max_abs_context_num_configs) {
			OL_from.reserve(_max_abs_context_num_configs) ;
			if (OL_from.capacity() < _max_abs_context_num_configs) 
				return 1 ;
			}
		}

	DFSANDcontainer *dfs_current = dfs_root ; // current DFS node
	while (NULL != dfs_current) {
#ifdef DEBUG_DFS_TRAVERSAL
std::string s, S ;
for (DFSANDcontainer *dfs = dfs_root ; NULL != dfs ; dfs = dfs->Child()) {
	if (S.length() > 0) S += ';' ;
	else  S = '\n' ;
	s.clear() ; dfs->SerializeCurrentState(s) ;
	S += s ;
	}
fwrite(S.c_str(), 1, S.length(), fpLOG) ;
fflush(fpLOG) ;
#endif // DEBUG_DFS_TRAVERSAL
		// advance AND/OR nodes on the openlist
		bool freshANDNode ;
		BucketElimination::Bucket *BC = dfs_current->FetchNextChildVar2Process_Proper(freshANDNode) ;

		// if all children at this level are exhausted, backtrack
		if (NULL == BC) {
			// compute value of current DFS node, as sum over OpenList; for each node in OpenList, incorporate accumulated (since last branching variable) ISw & Cost
			dfs_current->SubTreeValue() = ARP_nInfinity ;
			for (int32_t i = 0 ; i < dfs_current->OpenListSize() ; ++i) {
				SearchAndNode_WithPath *and_node = dfs_current->OpenListNode(i) ;
				ARE_Function_TableType v = and_node->SubTreeValue() ;
				ApplyFnCombinationOperator(v, and_node->AccumulatedISwSinceLastBranchingVariable()) ;
				ApplyFnCombinationOperator(v, and_node->AccumulatedCostSinceLastBranchingVariable()) ;
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(dfs_current->SubTreeValue(), dfs_current->SubTreeValue(), v) ;
				}
			// update parent DFS node with the value of this DFS node
			DFSANDcontainer *p = dfs_current->Parent() ;
			if (NULL != p) {
				SearchAndNode_WithPath *pAND = p->CurrentANDnode() ;
				ApplyFnCombinationOperator(pAND->SubTreeValue(), dfs_current->SubTreeValue()) ;
				}
			// TODO MAYBE : erase current context value from 'cn' to 'parent'
			// done; backtrack.
			if (dfs_current != dfs_root) delete dfs_current ;
			dfs_current = p ;
			if (NULL != dfs_current) dfs_current->AttachChild(NULL) ;
			continue ;
			}

		// if fresh AND node; check if h is exact; process completely here if yes.
		if (freshANDNode && EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(dfs_current->B()->V()) : false) {
			SearchAndNode_WithPath *can = dfs_current->CurrentANDnode() ;
			OL_from.clear() ; OL_from.push_back(can) ;
			_nNodesExactHComputed += OL_from.size() ;
			ARE_Function_TableType V = ARP_nInfinity ;
			ComputeFrontierValue(true, OL_from, V) ;
			// multiply with the current AND DFS node
			ApplyFnCombinationOperator(can->SubTreeValue(), V) ;
			OL_from.clear() ;
			dfs_current->MarkCurrentANDNodeProcessed() ; // mark this AND node is processed
			continue ;
			}

		// this is the DFS parent of the coming chain
		int32_t VC = BC->V() ;
		int32_t path_len_ = ComputeDistanceToNextDescendantBranchingVariable(VC) ;
		int32_t path_len = 1 + (path_len_ < 0 ? -path_len_ : path_len_) ;

		SearchAndNode_WithPath *can = dfs_current->CurrentANDnode() ;
		OL_from.clear() ; OL_from.push_back(can) ;
		// given next child, expand iteratively all children (descendants) until next branching variable.
		int32_t nLevels = 0 ;
		BucketElimination::Bucket *bC = BC ;

#ifdef _DEBUG
		if (++dfs_current->_nExpansions > dfs_current->OpenListSize()*dfs_current->nChildren()) {
			int bug_here = 1 ;
			}
#endif
		// expand the frontier down one level
expand_more :
#ifdef _DEBUG
{		int OL_to_size = OL_to.size() ;
		int OL_from_size = OL_from.size() ;
		int here = 1 ; }
#endif
		++nLevels ;
		int32_t vC = bC->V() ;
		int32_t nC = bC->nChildren() ;
		BucketElimination::Bucket *bParent = bC->ParentBucket() ;

		// 'OL_from' is the current frontier; it corresponds to the parent var/bucket of vC/bC; i.e. we are computing a new frontier corresponding to vC/bC, which is one level down

#ifdef _DEBUG
		int OL_to_size = OL_to.size() ;
		for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
			SearchAndNode_WithPath *an = OL_from[j] ;
			int isnan_res_an = isnan(an->SubTreeValue()) ;
			int isnan_res_g = isnan(an->Cost()) ;
			int isnan_res_h = isnan(an->h()) ;
			int isnan_res_ISw = isnan(an->ISw()) ;
			int isnan_res_aISw = isnan(an->AccumulatedISw()) ;
			int isnan_res_aISwLBV = isnan(an->AccumulatedISwSinceLastBranchingVariable()) ;
			int isnan_res_aC = isnan(an->AccumulatedCost()) ;
			int isnan_res_aClbv = isnan(an->AccumulatedCostSinceLastBranchingVariable()) ;
			if (0 != isnan_res_an || 0 != isnan_res_g || 0 != isnan_res_h || 0 != isnan_res_ISw || 0 != isnan_res_aISw || 0 != isnan_res_aISwLBV || 0 != isnan_res_aC || 0 != isnan_res_aClbv) 
				{ int is_bad = 1 ; }
			}
#endif // DEBUG

/*		// check and remove from the frontier all nodes that are -infinity (or 0 on normal scale)
		for (int32_t j = OL_from.size()-1 ; j >= 0 ; --j) {
			SearchAndNode_WithPath *an = OL_from[j] ;
			if (ARP_nInfinity == an->ISq()) { // IS q can be -inf when cost/h are -inf; since h is upper bound, when h=-inf, then the value of all nodes below is -inf.
				OL_from[j] = OL_from[OL_from.size()-1] ;
				OL_from.resize(OL_from.size()-1) ;
				if (nLevels > 1) 
					delete an ;
				}
			}*/
		if (0 == OL_from.size()) {
			// we should multiply 'can' with 0 (normal scale)
			can->SubTreeValue() = ARP_nInfinity ;
			dfs_current->MarkCurrentANDNodeProcessed() ; // done with this AND node; it is 0 and will not change.
			continue ;
			}

		// if bC has no children, generate new frontier, compute its value and rewind
		if (0 == nC) {
			// compute AND children of the frontier
			for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
				SearchAndNode_WithPath *an = OL_from[j] ;
				double newAndNodeISw = 1 == nLevels ? FnCombinationNeutralValue() : an->ISw() ;
				int32_t res_expand = ComputeAndChildren(can, dfs_current->idxCurrentChildnode(), an, bC, 1 == nLevels ? true : false, false, OL_to) ;
				}
			// delete OL_from nodes
			if (nLevels > 1) 
				for (SearchAndNode_WithPath *an : OL_from) delete an ;
			OL_from.clear() ;
			// compute frontier value
			ARE_Function_TableType V = ARP_nInfinity ;
			ComputeFrontierValue(false, OL_to, V) ;
			// multiply with the current AND DFS node
			ApplyFnCombinationOperator(can->SubTreeValue(), V) ;
			// some cleanup
			for (SearchAndNode_WithPath *an : OL_to) delete an ;
			OL_to.clear() ;
			continue ;
			}

		// if there is a limit on a number of branching variables for refined abstractions, switch to Knuth abstraction (i.e. merge all nodes into 1 node)
		if (nLevelsLimit() >= 0 && NumberOfAncestorBranchingVariables(vC) >= nLevelsLimit()) 
			no_AND_sorting = true ;
		else 
			no_AND_sorting = false ;

		// compute if we need to sort of not
		std::vector<int> & abs_context = MapVar2AbstractionContext(vC) ;
#if !defined(IGNORE_CONTEXT_ABS_SPECIALCODE)
		// if proper context-based, do it here; we can be better without sorting
		if (! no_AND_sorting && AbsSamplingTwoAndNodeCompare_CustomProper_DFS == _CompFn ? abs_context.size() > 0 : false) {
			// do : generate nodes, place at right idx, with merging if same abs node already exists
			// set up space where we store the nodes
			int32_t abs_context_num_configs = MapVar2AbstractionContextNumConfigs(vC) ;
			OL_to.resize(abs_context_num_configs, NULL) ;
			// fetch range of values for the child
			int32_t idxS = 0, idxE = _Problem->K(vC) ;
			if (_Problem->Value(vC) >= 0) 
				{ idxS = _Problem->Value(vC) ; idxE = idxS+1 ; }
			// compute AND children of the frontier
			BucketElimination::Bucket *b_an = bC->ParentBucket() ;
			for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
				SearchAndNode_WithPath *an = OL_from[j] ;
				// apply path assignment
				std::vector<int32_t> & pa = an->PathAssignment() ;
				int32_t n = pa.size() ; BucketElimination::Bucket *b = b_an ;
				for (; n > 0 && NULL != b ; b = b->ParentBucket()) 
					SetCurrentContextValue(b->V(), pa[--n]) ;
				// enumerate child domain
#ifdef USE_WMBq
				// define combination of all bChild intermediate functions here so that they get computed only once; this is ok, since all AND children generated here have the same context.
				double fInt = ARP_DBL_MAX ;
#endif
				for (int32_t iV = idxS ; iV < idxE ; iV++) {
					SearchAndNode_WithPath *anNew = new SearchAndNode_WithPath(this, vC, iV) ;
					if (NULL == anNew) 
						return 1 ;
					++_nNodesCreated ;
					anNew->SetClosestBranchingAncestor(can, dfs_current->idxCurrentChildnode()) ;
					anNew->_PathAssignment.reserve(nLevels) ;
					anNew->_PathAssignment = an->_PathAssignment ;
					anNew->_PathAssignment.push_back(iV) ;
					_CurrentContextValues[vC] = iV ;
					// compute label; accumulated cost(label)
					bC->ComputeCost(_CurrentContextValues, anNew->Cost()) ;
					bC->ComputeHeuristic(_CurrentContextValues, anNew->h()) ;
					if (nLevels > 1) {
						anNew->AccumulatedCostSinceLastBranchingVariable() = an->AccumulatedCostSinceLastBranchingVariable() ;
						anNew->AccumulatedISwSinceLastBranchingVariable() = an->AccumulatedISwSinceLastBranchingVariable() ;
						}
					else {
						anNew->AccumulatedCostSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
						anNew->AccumulatedISwSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
						}
					anNew->AccumulatedCost() = an->AccumulatedCost() ;
					ApplyFnCombinationOperator(anNew->AccumulatedCost(), anNew->Cost()) ;
					ApplyFnCombinationOperator(anNew->AccumulatedCostSinceLastBranchingVariable(), anNew->Cost()) ;
					anNew->ISw() = FnCombinationNeutralValue() ;
					anNew->AccumulatedISw() = an->AccumulatedISw() ;
					// compute IS q value for this node
					anNew->ISq() = anNew->AccumulatedCost() ;
					ApplyFnCombinationOperator(anNew->ISq(), anNew->h()) ;
					ApplyFnCombinationOperator(anNew->ISq(), anNew->AccumulatedISw()) ;
					if (ARP_nInfinity == anNew->ISq()) {
						delete anNew ;
						continue ; // ignore this node; ISq=-inf means c or h is -inf; h is upper boound.
						}
#ifdef USE_WMBq
					double q ;
					int wmbeq_res = bC->ComputeWMBEq(_CurrentContextValues, anNew->AccumulatedISw(), an->AccumulatedCost(), fInt, q) ;
					if (0 == wmbeq_res) 
						anNew->ISq() = q ;
#endif
					// compute output idx
					int32_t current_multiplier = 1 ; int32_t adr = 0 ;
					for (int32_t iCS = 0 ; iCS < abs_context.size() ; ++iCS) {
						int32_t vCS = abs_context[iCS] ;
						adr += current_multiplier*_CurrentContextValues[vCS] ;
						current_multiplier *= _Problem->K(vCS) ;
						}
					if (NULL == OL_to[adr]) 
						OL_to[adr] = anNew ;
					else 
						OL_to[adr] = (SearchAndNode_WithPath *) DoISmerge(anNew, OL_to[adr]) ;
					}
				}
			// delete OL_from nodes; at this point OL_from may contain node(s) from dfs_current; don't delete it (them) since they belong to dfs_current
			if (nLevels > 1) 
				for (SearchAndNode_WithPath *an : OL_from) delete an ;
			OL_from.clear() ;
			// copy OL_to to OL_from
			for (int32_t j = 0 ; j < OL_to.size() ; ++j) {
				SearchAndNode_WithPath *an = OL_to[j] ;
				if (NULL != an) OL_from.push_back(an) ;
				}
			// done; OL_from should be filled in with new frontier
#ifdef _DEBUG
int OL_from_size = OL_from.size() ;
int here = 1 ;
#endif
			goto frontier_expansion_done ;
			}
#endif

		// general case : first generate next-level frontier, then sort, then do IS merges.
#if !defined(IGNORE_PROPER_ABS_NOSORTING)
		if (AbsSamplingTwoAndNodeCompare_Knuth == _CompFn || AbsSamplingTwoAndNodeCompare_CustomProper_DFS == _CompFn) 
			no_AND_sorting = true ;
#endif
		// compute AND children of the frontier
		for (int32_t j = 0 ; j < OL_from.size() ; ++j) {
			SearchAndNode_WithPath *an = OL_from[j] ;
			double newAndNodeISw = 1 == nLevels ? FnCombinationNeutralValue() : an->ISw() ;
			int32_t res_expand = ComputeAndChildren(can, dfs_current->idxCurrentChildnode(), an, bC, 1 == nLevels ? true : false, false, OL_to) ;
			}
		// delete OL_from nodes; at this point OL_from may contain node(s) from dfs_current; don't delete it (them) since they belong to dfs_current
		if (nLevels > 1) 
			for (SearchAndNode_WithPath *an : OL_from) delete an ;
		OL_from.clear() ;
		int32_t idxS, idxE ; idxS = 0 ; idxE = -1 ;
		if (0 == OL_to.size()) 
			goto frontier_expansion_done ;
		else if (no_AND_sorting) {
			idxE = OL_to.size() ;
			}
		else {
			// sort output openlist
			if (OL_to.size() > 1) 
				QuickSort((void**) OL_to.data(), OL_to.size(), left, right, _CompFn) ;
			// go through the output openlist and do IS merges
			for (idxE = idxS+1 ; idxE < OL_to.size() ; ++idxE) {
				int32_t comp_res = (*_CompFn)(OL_to[idxS], OL_to[idxE]) ;
				if (0 != comp_res) {
					SearchAndNode_WithPath *an_merged = DoISmerge(OL_to, idxS, idxE) ;
					idxS = idxE ;
					OL_from.push_back(an_merged) ;
					}
				}
			}
		SearchAndNode_WithPath *an_merged ; an_merged = DoISmerge(OL_to, idxS, idxE) ;
		OL_from.push_back(an_merged) ;

// done with generating next-level frontier; note OL_from contains newly generated nodes
frontier_expansion_done :

		OL_to.clear() ;
		if (0 == OL_from.size()) { // this means frontier would be all 0's
			// we should multiply 'can' with 0 (normal scale)
			can->SubTreeValue() = ARP_nInfinity ;
			dfs_current->MarkCurrentANDNodeProcessed() ; // done with this AND node; its value is 0 (normal scale) and will not change
			continue ;
			}

		// if bC(vC) has exact h, compute values here
		if (EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(bC->V()) : false) {
			_nNodesExactHComputed += OL_from.size() ;
			ARE_Function_TableType V = ARP_nInfinity ;
			ComputeFrontierValue(true, OL_from, V) ;
			// multiply with the current AND DFS node
			ApplyFnCombinationOperator(can->SubTreeValue(), V) ;
			// delete OL_from nodes; at this point, OL_from contains newly generated nodes
			for (SearchAndNode_WithPath *an : OL_from) delete an ;
			OL_from.clear() ;
			continue ;
			}

		// if branching variable, create a new DFS node
		if (nC > 1) {
			DFSANDcontainer *dfs_node = new DFSANDcontainer(this, dfs_current, vC) ;
			if (NULL == dfs_node) return 1 ;
			if (0 != dfs_node->ImportOpenList(OL_from)) {
				// delete OL_from nodes; at this point, OL_from contains newly generated nodes
				for (SearchAndNode_WithPath *an : OL_from) delete an ;
				delete dfs_node ;
				return 1 ;
				}
			// note : we will keep the subtree value of all open list node as 1; when backtracking, we will incorporate accumulated (since last branching variable) ISw & Cost.
// bug			dfs_node->InitializeOpenListSubTreeValue() ;
			OL_from.clear() ;
			dfs_current = dfs_node ;
			}
		// else continue to next level down
		else {
			int32_t cv = bC->ChildVar(0) ;
			bC = MapVar2Bucket(cv) ;
			goto expand_more ;
			}
		}
	ApplyFnCombinationOperator(_Root.SubTreeValue(), dfs_root->SubTreeValue()) ;
	delete dfs_root ; dfs_root = NULL ;

#ifdef DEBUG_DFS_TRAVERSAL
fclose(fpLOG) ; fpLOG = NULL ;
#endif // DEBUG_DFS_TRAVERSAL

	return 0 ;
}

