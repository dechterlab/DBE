#include <vector>
#include <cmath>

#include "AVLtreeObjects.hxx"
#include "MBEworkspace.hxx"
#include "SearchSpaceNodes.hxx"
#include <AbsSamplingWorkspace.hxx>
#include <DFSANDcontainer.hxx>

bool AndOrSearchSpace::EXPLOIT_EXACT_H = true ;

static MTRand RNG_abs_cntxt_gen ;

int32_t AndOrSearchSpace::AbsSamplingWorkspace::ComputeRandAbstractionFactors(int n)
{
	_nRandAbs = n ;

	if (_RandAbsFactors.size() != _Problem->N()) {
		_RandAbsFactors.clear() ;
		if (_RandAbsFactors.capacity() < _Problem->N()) {
			_RandAbsFactors.reserve(_Problem->N()) ;
			if (_RandAbsFactors.capacity() < _Problem->N()) 
				return 1 ;
			}
		_RandAbsFactors.resize(_Problem->N(), 0) ;
		}

	if (_RandAbsFactorPerVar.size() != _Problem->N()) {
		_RandAbsFactorPerVar.clear() ;
		if (_RandAbsFactorPerVar.capacity() < _Problem->N()) {
			_RandAbsFactorPerVar.reserve(_Problem->N()) ;
			if (_RandAbsFactorPerVar.capacity() < _Problem->N()) 
				return 1 ;
			}
		_RandAbsFactorPerVar.resize(_Problem->N()) ;
		}

	for (int32_t i = 0 ; i < _Problem->N() ; ++i) {
		int k = _Problem->K(i) ;
		_RandAbsFactors[i] = 1 + (k <= 1 ? 0 : RNG_abs_cntxt_gen.randInt(k-1)) ;
		}

	for (int32_t v = 0 ; v < _Problem->N() ; ++v) {
		std::vector<int32_t> & abs_f = _RandAbsFactorPerVar[v] ;
		std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
		for (int32_t j = 0 ; j < abs_context.size() ; ++j) {
			int32_t u = abs_context[j] ;
			int k = _Problem->K(u) ;
			abs_f.push_back(1 + (k <= 1 ? 0 : RNG_abs_cntxt_gen.randInt(k-1))) ;
			}
		}

	return 0 ;
}


int32_t AndOrSearchSpace::AbsSamplingWorkspace::ComputeRandAbstractionFactors_Indv(void)
{
	for (int32_t v = 0 ; v < _Problem->N() ; ++v) {
		std::vector<int32_t> & abs_f = _RandAbsFactorPerVar[v] ;
		std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
		if (abs_f.size() != abs_context.size()) 
			continue ;
		for (int32_t j = 0 ; j < abs_context.size() ; ++j) {
			int32_t u = abs_context[j] ;
			int k = _Problem->K(u) ;
			abs_f[j] = 1 + (k <= 1 ? 0 : RNG_abs_cntxt_gen.randInt(k-1)) ;
			}
		}

	return 0 ;
}


AndOrSearchSpace::SearchAndNode_WithPath *AndOrSearchSpace::AbsSamplingWorkspace::DoISmerge(std::vector<SearchAndNode_WithPath*> & openList, int32_t idxS, int32_t idxE) 
// do IS merge of nodes [idxS, idxE); return the result node
{
	SearchAndNode_WithPath *anRes = NULL ;
	int32_t k, n = idxE - idxS, idxPicked = -1 ;
	if (n <= 0) {
		return anRes ;
		}
	double qSUM = ARP_DBL_MAX, qPicked = ARP_DBL_MAX, pPicked = ARP_DBL_MAX ;
	// notation : 
	//		Q_i is the (non-log) q-value of each node
	//		q_i if log(Q_i)
//	std::vector<double> qLIST ; qLIST.reserve(n) ; qLIST.resize(n) ;
	if (n > 1) {
		_nNodesISmerged += n-1 ;
		// compute the sum of q's
		for (k = idxS ; k < idxE ; ++k) {
			SearchAndNode_WithPath *a = openList[k] ;
//			double v = pow(10.0, a->ISq()) ;
			if (idxS == k) qSUM = a->ISq() ;
			else LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(qSUM,qSUM,a->ISq())
//			qLIST[k-idxS] = v ;
//			qSUM += v ;
			}
		// notation : qSUM is the log of sum of q_i's
		if (ARP_nInfinity == qSUM) {
			// all nodes have q=0; pick one uniformly randomly
			double n_ = (double) n ;
			pPicked = 1.0/n_ ;
			double r = n_ * _RNG.randExc() ;
			idxPicked = idxS + ((int)r) ; if (idxPicked >= idxE) idxPicked = idxE-1 ;
			qSUM = 0.0 ; qPicked = -log10(n_) ;
			}
		else {
			double r = _RNG.rand() ; // *qSUM ;
			for (k = idxS ; k < idxE ; ++k) {
				SearchAndNode_WithPath *a = openList[k] ;
				double p = pow(10.0, a->ISq() - qSUM) ;
				r -= p ;
				if (r <= 0.0) 
					{ break ; }
				}
			idxPicked = k ; if (idxPicked >= idxE) idxPicked = idxE-1 ;
			qPicked = openList[idxPicked]->ISq() ;
			pPicked = pow(10.0, qPicked - qSUM) ;
			}
		anRes = openList[idxPicked] ;
		// reweight stuff
		if (_Problem->FunctionsAreConvertedToLogScale()) {
			double log_p = log10(pPicked) ;
			anRes->ISw() = anRes->ISw() - log_p ;
			anRes->ISq() = anRes->ISq() - log_p ;
			anRes->AccumulatedISw() = anRes->AccumulatedISw() - log_p ;
			anRes->AccumulatedISwSinceLastBranchingVariable() = anRes->AccumulatedISwSinceLastBranchingVariable() - log_p ;
			}
		else {
			anRes->ISw() = anRes->ISw()/pPicked ;
			anRes->AccumulatedISw() = anRes->AccumulatedISw()/pPicked ;
			anRes->ISq() = anRes->ISq()/pPicked ;
			anRes->AccumulatedISwSinceLastBranchingVariable() = anRes->AccumulatedISwSinceLastBranchingVariable()/pPicked ;
			}
		}
	else 
		anRes = openList[idxS] ;
	for (k = idxS ; k < idxE ; ++k) {
		SearchAndNode_WithPath *a = openList[k] ;
		if (anRes != a) 
			{ delete a ; openList[k] = NULL ; }
		}
	return anRes ;
}


int32_t AndOrSearchSpace::AbsSamplingWorkspace::ComputeAndChildren(SearchAndNode_WithPath *ClosestAncestorAndNode, int32_t ClosestAncestorAndNodeChildIdx, SearchAndNode_WithPath *nParent, BucketElimination::Bucket *bChild, bool IsNewBeginning, bool SetFullPathAssignment, std::vector<SearchAndNode_WithPath*> & openListChild) 
{
	BucketElimination::Bucket *bParent = nParent->Bucket() ;
	int32_t v = NULL != bParent ? bParent->V() : -1 ;
	int32_t child = bChild->V() ;
	// fetch range of values for the child
	int32_t idxS = 0, idxE = _Problem->K(child) ;
	if (_Problem->Value(child) >= 0) 
		{ idxS = _Problem->Value(child) ; idxE = idxS+1 ; }
	// figure out pathassignment
	int32_t path_len = -1 ;
	if (IsNewBeginning) {
		if (SetFullPathAssignment) 
			nParent->ApplyPathAssignment(SetFullPathAssignment) ;
		path_len = 1 ;
		}
	else {
		nParent->ApplyPathAssignment(SetFullPathAssignment) ;
		path_len = nParent->_PathAssignment.size() + 1 ;
		}
	// enumarate the values of the OR node; generate AND children.
#ifdef USE_WMBq
	// define combination of all bChild intermediate functions here so that they get computed only once; this is ok, since all AND children generated here have the same context.
	double fInt = ARP_DBL_MAX ;
#endif
	for (int32_t iV = idxS ; iV < idxE ; iV++) {
		SearchAndNode_WithPath *anNew = new SearchAndNode_WithPath(this, child, iV) ;
		if (NULL == anNew) 
			return 1 ;
		++_nNodesCreated ;
		anNew->QPath() = nParent->QPath() ;
// DEBUGGGGG
//anNew->QPath() = 0.0 ;
		anNew->SetClosestBranchingAncestor(ClosestAncestorAndNode, ClosestAncestorAndNodeChildIdx) ;
		anNew->_PathAssignment.reserve(path_len) ;
		if (! IsNewBeginning) anNew->_PathAssignment = nParent->_PathAssignment ;
		anNew->_PathAssignment.push_back(iV) ;
		_CurrentContextValues[child] = iV ;
		// compute label; accumulated cost(label)
		bChild->ComputeCost(_CurrentContextValues, anNew->Cost()) ;
		bChild->ComputeHeuristic(_CurrentContextValues, anNew->h()) ;
		if (IsNewBeginning) {
			anNew->AccumulatedCostSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
			anNew->AccumulatedISwSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
			}
		else {
			anNew->AccumulatedCostSinceLastBranchingVariable() = nParent->AccumulatedCostSinceLastBranchingVariable() ;
			anNew->AccumulatedISwSinceLastBranchingVariable() = nParent->AccumulatedISwSinceLastBranchingVariable() ;
			}
		anNew->AccumulatedCost() = nParent->AccumulatedCost() ;
		ApplyFnCombinationOperator(anNew->AccumulatedCost(), anNew->Cost()) ;
		ApplyFnCombinationOperator(anNew->AccumulatedCostSinceLastBranchingVariable(), anNew->Cost()) ;
		anNew->ISw() = FnCombinationNeutralValue() ;
		anNew->AccumulatedISw() = nParent->AccumulatedISw() ;
		// compute IS q value for this node
		anNew->ISq() = anNew->AccumulatedCost() ;
		ApplyFnCombinationOperator(anNew->ISq(), anNew->h()) ;
		ApplyFnCombinationOperator(anNew->ISq(), anNew->AccumulatedISw()) ;
		if (anNew->QPath() < ARP_DBL_MAX) 
			ApplyFnCombinationOperator(anNew->ISq(), anNew->QPath()) ;
		if (ARP_nInfinity == anNew->ISq()) 
			{ delete anNew ; continue ; } // ignore this node; ISq=-inf means c or h is -inf; h is upper boound.
#ifdef USE_WMBq
		double q ;
		int wmbeq_res = bChild->ComputeWMBEq(_CurrentContextValues, anNew->AccumulatedISw(), nParent->AccumulatedCost(), fInt, q) ;
		if (0 == wmbeq_res) 
			anNew->ISq() = q ;
#endif
		// push to output list
		openListChild.push_back(anNew) ;
		}
	return 0 ;
}


int32_t AndOrSearchSpace::AbsSamplingWorkspace::ComputeAndChildrenEx(SearchAndNode_WithPath *ClosestAncestorAndNode, int32_t ClosestAncestorAndNodeChildIdx, SearchAndNode_WithPath *nParent, BucketElimination::Bucket *bChild, bool IsNewBeginning, std::vector<SearchAndNode_WithPath*> & openListChild) 
{
	++_nExpansionCalls ;
	BucketElimination::Bucket *bParent = nParent->Bucket() ;
	int32_t v = NULL != bParent ? bParent->V() : -1 ;
	int32_t child = bChild->V() ;
	// fetch range of values for the child
	int32_t idxS = 0, idxE = _Problem->K(child) ;
	if (_Problem->Value(child) >= 0) 
		{ idxS = _Problem->Value(child) ; idxE = idxS+1 ; }
	int32_t path_len = nParent->_PathAssignment.size() + 1 ;
//	// figure out pathassignment
//	if (IsNewBeginning) {
//		if (SetFullPathAssignment) 
//			nParent->ApplyPathAssignment(SetFullPathAssignment) ;
//		path_len = 1 ;
//		}
//	else {
//		nParent->ApplyPathAssignment(SetFullPathAssignment) ;
//		path_len = nParent->_PathAssignment.size() + 1 ;
//		}
	// enumarate the values of the OR node; generate AND children.
	for (int32_t iV = idxS ; iV < idxE ; iV++) {
		SearchAndNode_WithPath *anNew = new SearchAndNode_WithPath(this, child, iV) ;
		if (NULL == anNew) 
			return 1 ;
		++_nNodesCreated ;
		anNew->QPath() = nParent->QPath() ;
// DEBUGGGGG
//anNew->QPath() = 0.0 ;
		anNew->SetClosestBranchingAncestor(ClosestAncestorAndNode, ClosestAncestorAndNodeChildIdx) ;
		anNew->_PathAssignment.reserve(path_len) ;
//		if (! IsNewBeginning) 
		anNew->_PathAssignment = nParent->_PathAssignment ;
		anNew->_PathAssignment.push_back(iV) ;
		// compute label; accumulated cost(label)
		bChild->ComputeCostEx((anNew->_PathAssignment).data(), anNew->Cost()) ;
		bChild->ComputeHeuristicEx((anNew->_PathAssignment).data(), anNew->h()) ;
		if (IsNewBeginning) {
			anNew->AccumulatedCostSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
			anNew->AccumulatedISwSinceLastBranchingVariable() = FnCombinationNeutralValue() ;
			}
		else {
			anNew->AccumulatedCostSinceLastBranchingVariable() = nParent->AccumulatedCostSinceLastBranchingVariable() ;
			anNew->AccumulatedISwSinceLastBranchingVariable() = nParent->AccumulatedISwSinceLastBranchingVariable() ;
			}
		anNew->AccumulatedCost() = nParent->AccumulatedCost() ;
		ApplyFnCombinationOperator(anNew->AccumulatedCost(), anNew->Cost()) ;
		ApplyFnCombinationOperator(anNew->AccumulatedCostSinceLastBranchingVariable(), anNew->Cost()) ;
		anNew->ISw() = FnCombinationNeutralValue() ;
		anNew->AccumulatedISw() = nParent->AccumulatedISw() ;
		// compute IS q value for this node
		anNew->ISq() = anNew->AccumulatedCost() ;
		ApplyFnCombinationOperator(anNew->ISq(), anNew->h()) ;
		ApplyFnCombinationOperator(anNew->ISq(), anNew->AccumulatedISw()) ;
		if (anNew->QPath() < ARP_DBL_MAX) 
			ApplyFnCombinationOperator(anNew->ISq(), anNew->QPath()) ;
		if (ARP_nInfinity == anNew->ISq()) 
			{ delete anNew ; continue ; } // ignore this node; ISq=-inf means c or h is -inf; h is upper boound.
		// push to output list
		openListChild.push_back(anNew) ;
		}
	return 0 ;
}


int32_t AbsSamplingTwoAndNodeCompare_Knuth(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (A1->ClosestBranchingAncestor() < A2->ClosestBranchingAncestor()) 
		return -1 ;
	else if (A1->ClosestBranchingAncestor() > A2->ClosestBranchingAncestor()) 
		return 1 ;
/*	if (A1->Value() < A2->Value()) 
		return -1 ;
	else if (A1->Value() > A2->Value()) 
		return 1 ;*/

/*int32_t rescmp = AbsSamplineTwoAndNodeCompare(Obj1, Obj2) ;
if (0 != rescmp) {
	int bug = 1 ;
	printf("\nRRRR") ;
	}*/

	return 0 ;
}


int32_t AbsSamplingTwoAndNodeCompare_CustomProper(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;

	// check properness
	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;

	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (A1->ClosestBranchingAncestor() < A2->ClosestBranchingAncestor()) 
		return -1 ;
	else if (A1->ClosestBranchingAncestor() > A2->ClosestBranchingAncestor()) 
		return 1 ;
	if (A1->V() < 0) 
		return 0 ;

	// AND nodes have the same variable now.

//	// check if over the limit of number of nodes; if yes, switch to level-0 (Knuth) abstraction.
//	if (ws->nNodesCreated() > ws->nNodesCreatedLimitBeforeLevellingOff()) 
//		return 0 ;

/*	int32_t j = A1->V();
	int32_t nBVC = 0;
	while((j = ws->ClosestAncestorBranchingVariable(j)) >= 0){
		nBVC++;				
		}
	if (nBVC >= ws->nLevelsLimit()){
		return 0 ;
		}*/
	if (ws->nLevelsLimit() >=0 && ws->NumberOfAncestorBranchingVariables(A1->V()) >= ws->nLevelsLimit()) 
		return 0 ;

	// AND nodes have the same closest branching variable in AND/OR search tree; check the values of abs context.
	std::vector<int> & abs_context = A1->WS()->MapVar2AbstractionContext(A1->V()) ;
	int32_t nContext = abs_context.size() ;

	if (nContext <= 0)
		return 0 ;

	AndOrSearchSpace::SearchAndNode_WithPath *AwP1 = dynamic_cast<AndOrSearchSpace::SearchAndNode_WithPath*>(A1) ;
	AndOrSearchSpace::SearchAndNode_WithPath *AwP2 = dynamic_cast<AndOrSearchSpace::SearchAndNode_WithPath*>(A2) ;
	if (NULL != AwP1 && NULL != AwP2) {
		if (AwP1->_PathAssignment.size() != AwP2->_PathAssignment.size()) 
			// error; this should not happen; return based on whichever has smaller path.
			return AwP1->_PathAssignment.size() < AwP2->_PathAssignment.size() ? -1 : 1 ;
		int min = nContext ;
		if (min > AwP1->_PathAssignment.size()) min = AwP1->_PathAssignment.size() ;
		if (min > AwP2->_PathAssignment.size()) min = AwP2->_PathAssignment.size() ;
		std::vector<int32_t> & path1 = AwP1->_PathAssignment ;
		std::vector<int32_t> & path2 = AwP2->_PathAssignment ;
		int32_t i1 = path1.size()-1 ;
		int32_t i2 = path2.size()-1 ;
		for (int32_t ii = 0 ; ii < min ; ++ii) {
			if (path1[i1] < path2[i2]) return -1 ;
			else if (path1[i1] > path2[i2]) return 1 ;
			--i1 ; --i2 ;
			}
		return 0 ;
		}

	AndOrSearchSpace::SearchAndNode *a1 = A1, *a2 = A2 ;
	int32_t idxContext = 0 ;
	// printf(" START Var: %d %d A1: %d, A2: %d\n", A1->V(), A2->V(), A1->Value(), A2->Value());

	while (NULL != a1 && NULL != a2) {
		// check if a1/a2 value are given/set in abs context; if yes, compare values.
		if (a1->V() == abs_context[idxContext]) {
			// printf(" Var: %d %d a1: %d, a2: %d\n", a1->V(), a2->V(), a1->Value(), a2->Value());

			if (a1->Value() < a2->Value()) 
				return -1 ;
			else if (a1->Value() > a2->Value()) 
				return 1 ;
			// else values are same; do continue.
			if (++idxContext >= nContext) 
				break ; // no more context variables
			}
		// move up the tree
		AndOrSearchSpace::SearchOrNode *o1 = (AndOrSearchSpace::SearchOrNode *) a1->Parent(), *o2 = (AndOrSearchSpace::SearchOrNode *) a2->Parent() ;
		if (NULL == o1 || NULL == o2) break ;
		a1 = (AndOrSearchSpace::SearchAndNode *) o1->Parent() ; a2 = (AndOrSearchSpace::SearchAndNode *) o2->Parent() ;
		}
	// if(A1->h() != A2->h()){
	// 	printf("Not equal heuristics %g %g %d\n", A1->h(), A2->h(), nContext);
	// }

	return 0 ;
}


int32_t AbsSamplingTwoAndNodeCompare_CustomProper_DFS(void *Obj1, void *Obj2)
// this is a version of AbsSamplingTwoAndNodeCompare_CustomProper() specialized to DFS
// it should be faster than AbsSamplingTwoAndNodeCompare_CustomProper(); 
// we assume all AND nodes are instances of SearchAndNode_WithPath class.
// we assume properness by design; ie. will not check for it, since two AND nodes on the frontier are by definition/design descendants of the same closest branching variable.
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode_WithPath *A1 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj1 ;
	AndOrSearchSpace::SearchAndNode_WithPath *A2 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;

	if (A1->V() < 0) 
		return 0 ;

	// AND nodes have the same variable now.

	// check if over the limit of number of nodes; if yes, switch to level-0 (Knuth) abstraction.
	if (ws->nLevelsLimit() >=0 && ws->NumberOfAncestorBranchingVariables(A1->V()) >= ws->nLevelsLimit()) 
		return 0 ;

	// check the values of abs context.
	std::vector<int> & abs_context = A1->WS()->MapVar2AbstractionContext(A1->V()) ;
	int32_t nContext = abs_context.size() ;
	if (nContext <= 0)
		return 0 ;

	std::vector<int32_t> & path1 = A1->_PathAssignment ;
	std::vector<int32_t> & path2 = A2->_PathAssignment ;
	if (path1.size() != path2.size()) 
		// error; this should not happen; return based on whichever has smaller path.
		return path1.size() < path2.size() ? -1 : 1 ;

	// go from this var up the pseudo tree; compare values of variables in the context
	int32_t min = A1->_PathAssignment.size() ;
	BucketElimination::Bucket *b = A1->Bucket() ;
	int32_t iContext = 0 ;
	for (int32_t iPath = min-1 ; iPath >= 0 && iContext < nContext ; --iPath, b = b->ParentBucket()) {
		int32_t var = b->V() ;
		if (var != abs_context[iContext]) 
			// this var not in context; skip it
			continue ;
		// this var is in context; check values
		if (path1[iPath] < path2[iPath]) return -1 ;
		else if (path1[iPath] > path2[iPath]) return 1 ;
		++iContext ;
		}
	// all values the same; equal.
	return 0 ;
}


int32_t AbsSamplingTwoAndNodeCompare_Heuristic1(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;
	double hC = ws->heuristicCoefficient();

	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (A1->ClosestBranchingAncestor() < A2->ClosestBranchingAncestor()) 
		return -1 ;
	else if (A1->ClosestBranchingAncestor() > A2->ClosestBranchingAncestor()) 
		return 1 ;
	// AND nodes have the same variable
	if (A1->V() < 0) 
		return 0 ;
	//printf("A1->h: %g, A2->h: %g \n", A1->h(), A2->h());

	double g1_d = A1->h() + A1->AccumulatedCost();
	double g2_d = A2->h() + A2->AccumulatedCost();

	if(g1_d == -std::numeric_limits<float>::infinity() || g2_d == -std::numeric_limits<float>::infinity() ){
		return 0;
	}

	int32_t g1 = int(floor(g1_d*hC));
	int32_t g2 = int(floor(g2_d*hC));

	// printf("A1->f: %d, A2->f: %d  %g , %g \n", g1, g2, A1->h() + A1->AccumulatedCost(), A2->h() + A2->AccumulatedCost());

	if (g1 < g2) 
		return -1 ;
	else if (g1 > g2)
		return 1 ;

	return 0 ;
}

int32_t AbsSamplingTwoAndNodeCompare_Heuristic2(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;
	double hC = ws->heuristicCoefficient();

	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (A1->ClosestBranchingAncestor() < A2->ClosestBranchingAncestor()) 
		return -1 ;
	else if (A1->ClosestBranchingAncestor() > A2->ClosestBranchingAncestor()) 
		return 1 ;
	// AND nodes have the same variable
	if (A1->V() < 0) 
		return 0 ;
	//printf("A1->h: %g, A2->h: %g \n", A1->h(), A2->h());

	double g1_d = A1->h() + A1->AccumulatedCost();
	double g2_d = A2->h() + A2->AccumulatedCost();

	double g1 = floor(pow(10, g1_d-hC));
	double g2 = floor(pow(10, g2_d-hC));

	// int32_t g1 = int(g1_d*hC);
	// int32_t g2 = int(g2_d*hC);

	// printf("A1->f: %g, A2->f: %g  %g , %g \n", g1, g2, A1->h() + A1->AccumulatedCost(), A2->h() + A2->AccumulatedCost());

	if (g1 < g2) 
		return -1 ;
	else if (g1 > g2)
		return 1 ;

	return 0 ;
}

int32_t AbsSamplingTwoAndNodeCompare_Heuristic(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;
	double hC = ws->heuristicCoefficient();
	double hP = ws->heuristicPower();
	int level = A1->VarPosInOrder();
	// printf("\n%d %d %d\n", level, A1->Value(), A2->Value());

	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (A1->ClosestBranchingAncestor() < A2->ClosestBranchingAncestor()) 
		return -1 ;
	else if (A1->ClosestBranchingAncestor() > A2->ClosestBranchingAncestor()) 
		return 1 ;
	// AND nodes have the same variable
	if (A1->V() < 0) 
		return 0 ;
	//printf("A1->h: %g, A2->h: %g \n", A1->h(), A2->h());

	double g1_d = A1->h() + A1->AccumulatedCost();
	double g2_d = A2->h() + A2->AccumulatedCost();

	if(g1_d == -std::numeric_limits<float>::infinity() || g2_d == -std::numeric_limits<float>::infinity() ){
		// printf("Ndodhi\n");
		return 0;
	}

	// level = 0;
	double g1 = floor(pow(10, hP*(g1_d-(hC-level*log10(2)))));
	double g2 = floor(pow(10, hP*(g2_d-(hC-level*log10(2)))));
	// printf("hP: %lf\n", hP);

	// printf("A1->f: %g, A2->f: %g %g , %g \n", g1, g2, A1->h() + A1->AccumulatedCost(), A2->h() + A2->AccumulatedCost());

	if (g1 < g2)
		return -1 ;
	else if (g1 > g2)
		return 1 ;

	return 0 ;
}

int32_t AbsSamplingTwoAndNodeCompare_Unique(void *Obj1, void *Obj2)
{
	if (Obj1 == Obj2) 
		return 0 ;
	AndOrSearchSpace::SearchAndNode *A1 = (AndOrSearchSpace::SearchAndNode*) Obj1 ;
	AndOrSearchSpace::SearchAndNode *A2 = (AndOrSearchSpace::SearchAndNode*) Obj2 ;
	if (A1->VarPosInOrder() < A2->VarPosInOrder()) 
		return -1 ;
	else if (A1->VarPosInOrder() > A2->VarPosInOrder()) 
		return 1 ;
	if (Obj1 < Obj2) 
		return -1 ;
	else if (Obj1 > Obj2) 
		return 1 ;
	return 0 ;
}


int32_t AbsSamplingTwoAndNodeCompare_ContextNonProper(void *Obj1, void *Obj2)
// we assume all AND nodes are instances of SearchAndNode_WithPath class.
// we will compare the two nodes along the path from Obj1/Obj2 to the root; 
// some nodes on the path may not be in the context.
// we also assume that '_PathAssignment' member variable of a1/a2 is full path assigment, from a1/a2 to the root, not from a1/a2 to the closest branching variable.
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode_WithPath *A1 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj1 ;
	AndOrSearchSpace::SearchAndNode_WithPath *A2 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;

	if (A1->V() < 0) 
		return 0 ;
	if (A1->V() < A2->V()) 
		return -1 ;
	if (A1->V() > A2->V()) 
		return 1 ;

	// AND nodes have the same variable now.

	// check the values of abs context.
	std::vector<int> & abs_context = A1->WS()->MapVar2AbstractionContext(A1->V()) ;
	int32_t nContext = abs_context.size() ;
	if (nContext <= 0)
		return 0 ;
	int32_t iContext = 0 ;

	AndOrSearchSpace::SearchAndNode_WithPath *a1 = A1 ;
	AndOrSearchSpace::SearchAndNode_WithPath *a2 = A2 ;

next_block :

	std::vector<int32_t> & path1 = a1->_PathAssignment ;
	std::vector<int32_t> & path2 = a2->_PathAssignment ;
	if (path1.size() != path2.size()) 
		// error; this should not happen; return based on whichever has smaller path.
		return path1.size() < path2.size() ? -1 : 1 ;

	// go from this var up the pseudo tree; compare values of variables in the context
	int32_t path_len = path1.size() ;
	BucketElimination::Bucket *b = a1->Bucket() ;
	for (int32_t iPath = path_len-1 ; iPath >= 0 && iContext < nContext ; --iPath, b = b->ParentBucket()) {
		int32_t var = b->V() ;
		if (var != abs_context[iContext]) 
			// this var not in context; skip it
			continue ;
		// this var is in context; check values
		if (path1[iPath] < path2[iPath]) return -1 ;
		else if (path1[iPath] > path2[iPath]) return 1 ;
		++iContext ;
		}

	// 2018-09-28 KK : nonProper stores in PathAssignment the full assignment and in context the complete context; don't need to do 'next_block'.
	return 0 ;

	// all values the same; equal. 
	if (iContext >= nContext) 
		return 0 ;

	// continue up from next closest ancestor branching variable
	a1 = (AndOrSearchSpace::SearchAndNode_WithPath *) a1->ClosestBranchingAncestor() ;
	a2 = (AndOrSearchSpace::SearchAndNode_WithPath *) a2->ClosestBranchingAncestor() ;
	if (NULL == a1 || NULL == a2) 
		return 0 ;
	goto next_block ;
}


int32_t AbsSamplingTwoAndNodeCompare_RandCntxt(void *Obj1, void *Obj2)
// we assume all AND nodes are instances of SearchAndNode_WithPath class.
// we will compare the two nodes along the path from Obj1/Obj2 to the root; 
// some nodes on the path may not be in the context.
// we also assume that '_PathAssignment' member variable of a1/a2 is full path assigment, from a1/a2 to the root, not from a1/a2 to the closest branching variable.
{
	if (Obj1 == Obj2) 
		return 0 ;

	AndOrSearchSpace::SearchAndNode_WithPath *A1 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj1 ;
	AndOrSearchSpace::SearchAndNode_WithPath *A2 = (AndOrSearchSpace::SearchAndNode_WithPath*) Obj2 ;
	AndOrSearchSpace::AbsSamplingWorkspace *ws = A1->WS() ;

	if (A1->V() < 0) 
		return 0 ;
	if (A1->V() < A2->V()) 
		return -1 ;
	if (A1->V() > A2->V()) 
		return 1 ;

	// AND nodes have the same variable now.

	int a1 = A1->WS()->ComputeNodeAbstrationID(*A1) ;
	int a2 = A2->WS()->ComputeNodeAbstrationID(*A2) ;
	if (a1 == a2) 
		return 0 ;
	return a1 < a2 ? -1 : 1 ;
}

