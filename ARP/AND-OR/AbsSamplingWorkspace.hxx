#ifndef AbsSamplingWorkspace_HXX_INCLUDED
#define AbsSamplingWorkspace_HXX_INCLUDED

#include <vector>
#include <utility>

#include "AVLtreeObjects.hxx"
#include "AVLtreeInt64.hxx"
#include "MersenneTwister.h"
#include "MBEworkspace.hxx"
#include "SearchSpaceNodes.hxx"

//#define DEBUG_AND_OR_CONSTRUCTION
//#define DEBUG_AND_OR_COMPUTATION
//#define ABS_CONTEXT_IS_EVERYTHING

typedef int32_t (AbsSamplingCompFn)(void *Obj1, void *Obj2) ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_Unique ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_Knuth ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_CustomProper ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_CustomProper_DFS ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_Heuristic ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_ContextNonProper ;
AbsSamplingCompFn AbsSamplingTwoAndNodeCompare_RandCntxt ;

namespace AndOrSearchSpace
{

class SearchAndNode_WithPath ;

extern bool EXPLOIT_EXACT_H ;

class AbsSamplingWorkspace : public BucketElimination::MBEworkspace
{
protected :

	// root of the AND/OR search space; var is -1; children are roots of the pseudo-tree as OR nodes.
	SearchAndNode _Root ;

	// comparison fn for comparing two AND nodes; this defines the abstraction.
	AbsSamplingCompFn *_CompFn ;

public :

	inline SearchAndNode & Root(void) { return _Root ; }

protected :

	// limit on the number of nodes created before we switch to level-0 (Knuth) abstraction (this is the abstraction with smallest number of nodes).
	// before this point, the number of nodes created is (worst case) exponential (in some parameters); once this limit is reached, it is linear.
	int64_t _nNodesCreatedLimitBeforeLevellingOff ;
	int32_t _nLevelsLimit ;
	double _heuristicCoefficient;
	double _heuristicPower;

public :

	inline int64_t & nNodesCreatedLimitBeforeLevellingOff(void) { return _nNodesCreatedLimitBeforeLevellingOff ; }	
	inline int32_t & nLevelsLimit(void) { return _nLevelsLimit; }
	inline double & heuristicCoefficient(void) {return _heuristicCoefficient; }
	inline double & heuristicPower(void) {return _heuristicPower; }
protected :

	// number of nodes in the AND/OR graph; including root node.
	int32_t _nNodesInTree ;

	// number of nodes generated during tree construction.
	int32_t _nNodesCreated ;

	// number of nodes Importance Sampling merged.
	int32_t _nNodesISmerged ;

	// number of nodes with exact H computed, without expansion.
	int32_t _nNodesExactHComputed ;

	// number of Expansion fn calls
	int32_t _nExpansionCalls ;

	// number of Expansion fn calls
	int32_t _nTotalSumFrontierSize ;

	// number of variable expansions
	int32_t _nNumVarExpansions ;

	// number of variables in the bucket tree that are partitioned.
	int32_t _nVarsPartitioned ;

public :

	inline int32_t nNodesInTree(void) { return _nNodesInTree ; }
	inline int32_t nNodesCreated(void) { return _nNodesCreated ; }
	inline int32_t nNodesISmerged(void) { return _nNodesISmerged ; }
	inline int32_t nNodesExactHComputed(void) { return _nNodesExactHComputed ; }
	inline int32_t nVarsPartitioned(void) { return _nVarsPartitioned ; }
	inline int32_t nExpansionCalls(void) { return _nExpansionCalls ; }
	inline int32_t nTotalSumFrontierSize(void) { return _nTotalSumFrontierSize ; }
	inline int32_t nNumVarExpansions(void) { return _nNumVarExpansions ; }

protected :

	// a closest branching variable for each variable.
	std::vector<int32_t> _ClosestAncestorBranchingVariable ;

	// number of ancestor branching variables. for each variables.
	std::vector<int32_t> _NumberOfAncestorBranchingVariables ;

public :

	inline int32_t ClosestAncestorBranchingVariable(int32_t var) const { return var >= 0 ? _ClosestAncestorBranchingVariable[var] : -1 ; }
	inline int32_t NumberOfAncestorBranchingVariables(int32_t var) const { return var >= 0 ? _NumberOfAncestorBranchingVariables[var] : -1 ; }

	int32_t ComputeClosestBranchingVariables(void)
	{
		_ClosestAncestorBranchingVariable.clear() ;
		_NumberOfAncestorBranchingVariables.clear() ;
		if (_ClosestAncestorBranchingVariable.capacity() < _Problem->N()) {
			_ClosestAncestorBranchingVariable.reserve(_Problem->N()) ;
			if (_ClosestAncestorBranchingVariable.capacity() < _Problem->N()) 
				return 1 ;
			}
		if (_NumberOfAncestorBranchingVariables.capacity() < _Problem->N()) {
			_NumberOfAncestorBranchingVariables.reserve(_Problem->N()) ;
			if (_NumberOfAncestorBranchingVariables.capacity() < _Problem->N()) 
				return 2 ;
			}
		_ClosestAncestorBranchingVariable.resize(_Problem->N(), -1) ;
		_NumberOfAncestorBranchingVariables.resize(_Problem->N(), -1) ;
		if (_Problem->N() != _ClosestAncestorBranchingVariable.size()) 
			return 3 ;
		if (_Problem->N() != _NumberOfAncestorBranchingVariables.size()) 
			return 4 ;

		for (int32_t i = 0 ; i < _Problem->N() ; ++i) {
			int32_t v = _Problem->VarOrdering_VarList()[i] ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			BucketElimination::Bucket *bParent = b->ParentBucket() ;
			if (NULL == bParent) 
				{ _ClosestAncestorBranchingVariable[v] = -1 ; continue ; }
			int32_t vParent = bParent->Var(0) ; ;
			if (bParent->nChildren() > 1) 
				_ClosestAncestorBranchingVariable[v] = vParent ;
			else 
				_ClosestAncestorBranchingVariable[v] = _ClosestAncestorBranchingVariable[vParent] ;
			}

		for (int32_t v = 0 ; v < _Problem->N() ; ++v) {
			int32_t nBVC = 0 ;
			int32_t j = v ;
			while ((j = ClosestAncestorBranchingVariable(j)) >= 0) {
				nBVC++ ;				
				}
			_NumberOfAncestorBranchingVariables[v] = nBVC ;
			}

		return 0 ;
	}

	int32_t MaxNumBranchingVariablesInChain() const {
		int32_t maxNBVC = -1;
		for(int32_t i = 0; i < _Problem->N(); i++){
			int32_t nBVC = 0;
			int32_t j = i;
			while((j = ClosestAncestorBranchingVariable(j)) >= 0){
				nBVC++;				
			}
			if(nBVC > maxNBVC){
				maxNBVC = nBVC;
			}
		}

		return maxNBVC;
	}

	int32_t NumBranchingVariables() const {
		int32_t NBV = 0;
		for(int32_t i = 0; i < _Problem->N(); i++){
			int32_t v = _Problem->VarOrdering_VarList()[i] ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			if(b->nChildren() > 1) NBV++;
		}

		return NBV;
	}

	void ExportPseudotree(int32_t * & ptArray) const {
		ptArray = new int32_t[_Problem->N()];
		for (int32_t i = 0 ; i < _Problem->N() ; ++i) {
			int32_t v = _Problem->VarOrdering_VarList()[i] ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			BucketElimination::Bucket *bParent = b->ParentBucket() ;
			if (NULL == bParent){
				ptArray[v] = -1 ; 
				continue ; 
			}
			int32_t vParent = bParent->Var(0) ;
			ptArray[v] = vParent;
		}
	}

	int64_t EstimateSampleTreeSize() const {
        int64_t totalSTS = 0;
        
		for(int32_t i = 0; i < _Problem->N(); i++){
			int32_t j = i;
            int64_t layerSize = 1;
            
            int32_t v1 = _Problem->VarOrdering_VarList()[j] ;
            BucketElimination::Bucket *b = MapVar2Bucket(v1) ;

            if(b == NULL) continue;

            if(b->nChildren() > 1){
            	layerSize *= _Problem->K(v1);
        	}            
			while((j = ClosestAncestorBranchingVariable(j)) >= 0){
                int32_t v = _Problem->VarOrdering_VarList()[j] ;
                layerSize *= _Problem->K(v);
                //printf("Layer %lld\n", layerSize);
			}

			totalSTS += layerSize;
			//printf("Total: %lld\n", totalSTS );
		}
        return 2*totalSTS;        
    }


	std::pair <int64_t, int64_t> EstimateSampleTreeSizeFullContext(void) {
        int64_t AND_totalSTS = 0;
        int64_t OR_totalSTS = 0;
        int64_t AND_totalExpanded = 0;
        int64_t OR_totalExpanded = 0;
/* this does not compile with VS2015; array declaration below not ok.
        int64_t AND_layerSize[_Problem->N()];
        int64_t OR_layerSize[_Problem->N()];
        int64_t AND_total_layerSize[_Problem->N()];
        int64_t OR_total_layerSize[_Problem->N()];
        
        int32_t k[_Problem->N()];

        for(int32_t v  = 0; v < _Problem->N(); v++){
        	AND_layerSize[v] = 0;
        	OR_layerSize[v] = 0;
        	AND_total_layerSize[v] = 0;
        	OR_total_layerSize[v] = 0;
        	if(_Problem->Value(v) != -1)
        	{ 
        		k[v] = 1;
        	}
        	else
        	{
        		k[v] = _Problem->K(v);
        	}	            
        }     			

		for(int32_t i = 0; i < _Problem->N(); i++){
            int32_t v = _Problem->VarOrdering_VarList()[i];
	        int32_t cABV = ClosestAncestorBranchingVariable(v);
	        int32_t vParent;
            std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
            BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) 
			{
				continue ;
			}
			BucketElimination::Bucket *bParent = b->ParentBucket() ;
			if(bParent == NULL){
				vParent = -1;
			}
			else
			{
				vParent = bParent->Var(0);
			}

			if(vParent == -1)
			{
					AND_layerSize[v] = k[v];
            		OR_layerSize[v] = 1;

            		AND_total_layerSize[v] = k[v];
            		OR_total_layerSize[v] = 1;
			}
            // else if(abs_context.size() <= 1)
            // {
            // 	if(cABV == -1){
            // 		AND_layerSize[v] = k[v];
            // 		OR_layerSize[v] = 1;
            // 		AND_total_layerSize[v] = AND_layerSize[vParent]*k[v];
            // 		OR_total_layerSize[v] = AND_layerSize[vParent];
            // 	}
            // 	else
            // 	{
            // 		AND_layerSize[v] = AND_layerSize[cABV]*k[v];
            // 		OR_layerSize[v] = AND_layerSize[cABV];
            // 		AND_total_layerSize[v] = AND_layerSize[vParent]*k[v];
            // 		OR_total_layerSize[v] = AND_layerSize[vParent];
            // 	}
            // }
            else
            {	
            	int64_t num_contexts = 1;
            	for(int32_t j = 1; j < abs_context.size(); j++)
            	{
            		num_contexts*= k[abs_context[j]];
            	}

            	if(cABV == -1){
            		AND_layerSize[v] = num_contexts*k[v];
            		OR_layerSize[v] = num_contexts;
            		AND_total_layerSize[v] = AND_layerSize[vParent]*k[v];
            		OR_total_layerSize[v] = AND_layerSize[vParent];
            	}
            	else
            	{
            		AND_layerSize[v] = AND_layerSize[cABV]*num_contexts*k[v];
            		OR_layerSize[v] = AND_layerSize[cABV]*num_contexts;
            		AND_total_layerSize[v] = AND_layerSize[vParent]*k[v];
            		OR_total_layerSize[v] = AND_layerSize[vParent];
            	}
            }

			AND_totalSTS += AND_layerSize[v];
			OR_totalSTS += OR_layerSize[v];
			AND_totalExpanded += AND_total_layerSize[v];
			OR_totalExpanded += OR_total_layerSize[v];
			//printf("\nlayerSize: %lld  %d %d", layerSize[v], v, cABV);
		}

		// printf("\nTotal: %lld  %lld %lld\n", AND_totalSTS, OR_totalSTS, AND_totalSTS + OR_totalSTS + 1);
		// printf("Total Expanded: %lld  %lld %lld\n", AND_totalExpanded, OR_totalExpanded, AND_totalExpanded + OR_totalExpanded + 1);
*/
		std::pair<int64_t, int64_t> return_pair;
		return_pair = std::make_pair(AND_totalExpanded + OR_totalExpanded + 1, AND_totalSTS + OR_totalExpanded + 1);
		return return_pair;
        // return AND_totalSTS + OR_totalExpanded + 1;   
    }

	std::pair <int64_t, int64_t> EstimateTreeSizeContextMinimal() {
        int64_t AND_totalSTS = 0;
        int64_t OR_totalSTS = 0;
        int64_t AND_totalExpanded = 0;
        int64_t OR_totalExpanded = 0;
/* this does not compile with VS2015; array declaration below not ok.
        int64_t AND_layerSize[_Problem->N()];
        int64_t OR_layerSize[_Problem->N()];
        int64_t AND_total_layerSize[_Problem->N()];
        int64_t OR_total_layerSize[_Problem->N()];
        
        int32_t k[_Problem->N()];

        for(int32_t v  = 0; v < _Problem->N(); v++){
        	AND_layerSize[v] = 0;
        	OR_layerSize[v] = 0;
        	AND_total_layerSize[v] = 0;
        	OR_total_layerSize[v] = 0;
        	if(_Problem->Value(v) != -1)
        	{ 
        		k[v] = 1;
        	}
        	else
        	{
        		k[v] = _Problem->K(v);
        	}	            

			}     			

		for(int32_t i = 0; i < _Problem->N(); i++){
            int32_t v = _Problem->VarOrdering_VarList()[i];
	        int32_t vParent;
            std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
            BucketElimination::Bucket *b = MapVar2Bucket(v) ;
            if (NULL == b) continue ;
			BucketElimination::Bucket *bParent = b->ParentBucket() ;
			if(bParent == NULL){
				vParent = -1;
			}
			else
			{
				vParent = bParent->Var(0);
			}

            //b->ComputeANDContextSignature();
			const int32_t *sig = b->ANDContextSignature(); // true AND context
			int nSig = b->ANDContextWidth();
			if (nSig < 0) continue ;
			std::vector<int32_t> context_vars;
			context_vars.clear();
			for (int32_t j = 0 ; j < nSig ; ++j) {
				int32_t u = sig[j] ;
				if (u == v) continue ;
				context_vars.push_back(u) ;
			}			

			if(vParent == -1)
			{
					AND_layerSize[v] = k[v];
            		OR_layerSize[v] = 1;

            		AND_total_layerSize[v] = k[v];
            		OR_total_layerSize[v] = 1;
			}
            else
            {	
            	int64_t num_contexts = 1;
            	for(int32_t j = 0; j < context_vars.size(); j++)
            	{
            		num_contexts*= k[context_vars[j]];
            	}

           		AND_layerSize[v] = num_contexts*k[v];
           		OR_layerSize[v] = num_contexts;
           		AND_total_layerSize[v] = AND_layerSize[vParent]*k[v];
           		OR_total_layerSize[v] = AND_layerSize[vParent];
            }

			AND_totalSTS += AND_layerSize[v];
			OR_totalSTS += OR_layerSize[v];
			AND_totalExpanded += AND_total_layerSize[v];
			OR_totalExpanded += OR_total_layerSize[v];
			//printf("\nlayerSize: %lld  %d %d", layerSize[v], v, cABV);
		}

		// printf("\nTotal ContextMinimal: %lld  %lld %lld\n", AND_totalSTS, OR_totalSTS, AND_totalSTS + OR_totalSTS + 1);
		// printf("Total ContextMinimal Expanded: %lld  %lld %lld\n", AND_totalExpanded, OR_totalExpanded, AND_totalExpanded + OR_totalExpanded + 1);
		// printf("\nTotal: %lld  %lld\n", totalSTS, OR_totalSTS);
*/
		std::pair<int64_t, int64_t> return_pair;
		return_pair = std::make_pair(AND_totalExpanded + OR_totalExpanded + 1, AND_totalSTS + OR_totalExpanded + 1);
		return return_pair;
        // return AND_totalSTS + OR_totalExpanded + 1;        
    }

protected :

	// a boolean flag indicating whether any of its descendants have MB partitioning.
	// if none of its descendats has partitioning, then the h computed at this node is exact, assuming this node and all ancestors are instantiated.
	// this means we don't have to expand anything below this variable, since h at this variable is exact.
	std::vector<signed char> _SubtreeHasPartitioning ;

public :

	inline bool SubtreeHasPartitioning(int32_t var) const { return var >= 0 ? 0 != _SubtreeHasPartitioning[var] : false ; }

	int32_t ComputeSubtreePartitioning(void)
	{
		_nVarsPartitioned = -1 ;
		_SubtreeHasPartitioning.clear() ;
		_SubtreeHasPartitioning.reserve(_Problem->N()) ;
		if (_SubtreeHasPartitioning.capacity() < _Problem->N()) 
			return 1 ;
		_SubtreeHasPartitioning.resize(_Problem->N(), -1) ;
		if (_Problem->N() != _SubtreeHasPartitioning.size()) 
			return 2 ;
		_nVarsPartitioned = 0 ;
		// traverse the bucket tree in the computation_order; we assume one is computed...
		for (int32_t i = _lenComputationOrder - 1 ; i >= 0 ; --i) {
			int32_t idxToCompute = BucketOrderToCompute(i) ;
			if (idxToCompute < 0) continue ;
			BucketElimination::Bucket *b = getBucket(idxToCompute) ;
			if (NULL == b) continue ;
			int32_t v = b->V() ;
//		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
//			int32_t v = _Problem->VarOrdering_VarList()[i] ;
//			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
//			if (NULL == b) continue ;
			if (b->nMiniBuckets() > 1) 
				++_nVarsPartitioned ;
			BucketElimination::Bucket *bParent = b->ParentBucket() ;
			int32_t vParent = NULL != bParent ? bParent->V() : -1 ;
			if (_SubtreeHasPartitioning[v] > 0) {
				if (vParent >= 0) 
					_SubtreeHasPartitioning[vParent] = _SubtreeHasPartitioning[v] ;
				continue ;
				}
			if (b->nMiniBuckets() > 1) {
//				_SubtreeHasPartitioning[v] = 1 ;
				if (vParent >= 0 ? _SubtreeHasPartitioning[vParent] <= 0 : false) 
					_SubtreeHasPartitioning[vParent] = 1 ;
				}
			else 
				_SubtreeHasPartitioning[v] = 0 ;
			}
#ifdef DEBUG_AND_OR_CONSTRUCTION
		printf("\nnVarsPartitioned=%d", _nVarsPartitioned) ;
#endif
		return 0 ;
	}

protected :

	// for each variable in the bucket tree, its abstraction context.
	// abs context is a vector of variables, along the path from this variable to the root. the order is very important; note that some variables along the path to the root may not be in the context, e.g. if the var is not in the scope of the h function of the node.
	// e.g. this variable, and then all variables from closest ancestor to the root.
	std::vector<std::vector<int32_t>> _MapVar2AbstractionContext ;
	std::vector<int32_t> _MapVar2AbstractionContextNumConfigs ; // number of different configurations of the abs context; equals the product of domain sizes of all the context variables.
	int32_t _maxContextSize, _max_abs_context_num_configs ;
	int32_t _PTInducedWidth ;

	// for randomized abstractions
	std::vector<int32_t> _RandAbsFactors ; // for each var, its weight
	std::vector<std::vector<int32_t>> _RandAbsFactorPerVar ;
	int32_t _nRandAbs ; // number randomized abstractions per variable

public :

	inline std::vector<int32_t> & MapVar2AbstractionContext(int var) { return _MapVar2AbstractionContext[var] ; }
	inline int32_t MapVar2AbstractionContextNumConfigs(int var) { return _MapVar2AbstractionContextNumConfigs[var] ; }
	inline int32_t MaxContextSize(void) { return _maxContextSize; }
	inline int32_t MaxAbsContextNumConfigs(void) { return _max_abs_context_num_configs; }
	inline int32_t PTInducedWidth(void) { return _PTInducedWidth; }

	int32_t ComputePTInducedWidth(){
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			b->ComputeSignature(true) ; // true OR context + current variable; all functions included
			const int32_t *sig = b->Signature() ; 
			int nSig = b->Width() ;
			if (nSig < 0) 
				continue ;

			std::vector<int32_t> context_vars;
			for (int32_t j = 0 ; j < nSig ; ++j) {
				int32_t u = sig[j] ;
				if (u == v) continue ;
				context_vars.push_back(u) ; 
			}

			int32_t nContext = context_vars.size();
		
			if(nContext > _PTInducedWidth){
				_PTInducedWidth = nContext;
			}
		}	

	}

	// n is number of variables from union(given_var,context(given_var)) to include in the abs context
	int32_t ComputeRandAbstractionFactors(int n) ;
	// generate random values for '_RandAbsFactorPerVar'; assume size is correct.
	int32_t ComputeRandAbstractionFactors_Indv(void) ;

	int32_t ComputeNodeAbstrationID(AndOrSearchSpace::SearchAndNode_WithPath & N)
	{
		int32_t v = N.V() ;

		// check if over the limit of number of nodes; if yes, switch to level-0 (Knuth) abstraction.
		if (_nLevelsLimit >= 0 && _NumberOfAncestorBranchingVariables[v] >= _nLevelsLimit) 
			return 0 ;

		std::vector<int32_t> & abs_context = MapVar2AbstractionContext(v) ;
		std::vector<int32_t> & abs_f = _RandAbsFactorPerVar[v] ;
		int32_t nContext = abs_context.size() ;
		if (nContext <= 0)
			return 0 ;
		std::vector<int32_t> & path = N._PathAssignment ;
		int32_t path_len = path.size() ;
		BucketElimination::Bucket *b = N.Bucket() ;
		int32_t iContext = 0 ;
		int64_t x = 0 ;
		for (int32_t iPath = path_len-1 ; iPath >= 0 && iContext < nContext ; --iPath, b = b->ParentBucket()) {
			int32_t var = b->V() ;
			if (var != abs_context[iContext]) 
				// this var not in context; skip it
				continue ;
//			x += _RandAbsFactors[var]*path[iPath] ;
			x += abs_f[iContext]*path[iPath] ;
			++iContext ;
			}
		int32_t abs = x % _nRandAbs ;
		return abs ;
	}

	int32_t ComputeAbstractionContextReal(int n, bool ignoreclosestbranchingvariable)
	// n is number of variables from union(given_var,context(given_var)) to include in the abs context
	{
		_MapVar2AbstractionContext.clear() ;
		_MapVar2AbstractionContext.reserve(_Problem->N()) ;
		if (_MapVar2AbstractionContext.capacity() < _Problem->N()) 
			return 1 ;
		_MapVar2AbstractionContext.resize(_Problem->N()) ;
		if (_Problem->N() != _MapVar2AbstractionContext.size()) 
			return 2 ;

		_MapVar2AbstractionContextNumConfigs.clear() ;
		_MapVar2AbstractionContextNumConfigs.reserve(_Problem->N()) ;
		if (_MapVar2AbstractionContextNumConfigs.capacity() < _Problem->N()) 
			return 1 ;
		_MapVar2AbstractionContextNumConfigs.resize(_Problem->N(), 1) ;
		if (_Problem->N() != _MapVar2AbstractionContextNumConfigs.size()) 
			return 2 ;

		// if n=0 we are done; but create/default_fill _MapVar2AbstractionContext/_MapVar2AbstractionContextNumConfigs
		if (n <= 0) 
			return 0 ;

#ifdef ABS_CONTEXT_IS_EVERYTHING
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
			std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
			abs_context.clear() ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			for (; NULL != b ; b = b->ParentBucket()) {
				for (int32_t j = 0 ; j < b->nVars() ; j++) {
					int32_t u = b->Var(j) ;
					abs_context.push_back(u) ;
					}
				}
			}
#endif 

		// compute max signature over all buckets
		int32_t maxNSig = -1 ;
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
//			int32_t v = (_Problem->VarOrdering_VarList())[i] ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			// get signature
			b->ComputeSignature(0x06) ; // 0x06 = AND context (aug+int functions included); 0x07 = true context (all functions included). AND context is the scope of the h of the node.
			const int32_t *sig = b->Signature() ;
			int nSig = b->Width() ;
			if (nSig < 0) 
				continue ;
			if (nSig > maxNSig) 
				maxNSig = nSig ;
			}
		int32_t w1 = VarOrdering_InducedWidth() ;
		int32_t max_context_size = 1 + (maxNSig > w1 ? maxNSig : w1) ;

		// two helper arrays
		std::vector<int32_t> context_vars, context_var_indeces ;
		context_vars.reserve(max_context_size) ;
		context_var_indeces.reserve(max_context_size) ;
		if (context_vars.capacity() != context_var_indeces.capacity()) 
			return 1 ;
		if (context_vars.capacity() != max_context_size) 
			return 1 ;

		int32_t left[32], right[32] ;
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
//			int32_t v = (_Problem->VarOrdering_VarList())[i] ;
			std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
			abs_context.clear() ;
			int32_t & abs_context_num_configs = _MapVar2AbstractionContextNumConfigs[v] ;
			abs_context_num_configs = 1 ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			// get signature; should be computed in previous paragraph...
			const int32_t *sig = b->Signature() ;
			int nSig = b->Width() ;
			if (nSig < 0) 
				continue ;
			// build context; build sort keys for context vars
			int32_t vClosestAncestor = _ClosestAncestorBranchingVariable[v] ;
			int32_t posClosestAncestor = vClosestAncestor >= 0 ? _Problem->VarOrdering_VarPos()[vClosestAncestor] : -1 ;
			context_vars.clear() ; context_var_indeces.clear() ;
			for (int32_t j = 0 ; j < nSig ; ++j) {
				int32_t u = sig[j] ;
				if (u == v) continue ;
				int32_t pos = _Problem->VarOrdering_VarPos()[u] ;
				if (! ignoreclosestbranchingvariable && pos <= posClosestAncestor) 
					continue ; // skip variables at or before closest branching ancestor
				context_vars.push_back(u) ;
				context_var_indeces.push_back(-pos) ;
				}

			if (context_var_indeces.size() != context_vars.size()) {
				printf("ERROR - Different sizes \n");
				printf("\nCxtSize: %d CtxIdxS: %d\n", (int) context_vars.size(), (int) context_var_indeces.size());
				return 1 ;
				}

			if (context_var_indeces.size() > 1) {
				QuickSortLong(context_var_indeces.data(), context_var_indeces.size(), context_vars.data(), left, right) ;
				}

			// add this var to context
			abs_context.push_back(v) ;
			abs_context_num_configs *= _Problem->K(v) ;
			int n_ = n-1 ; // number of variables remaining
			for (int32_t j = 0 ; j < n_ && j < context_vars.size() ; ++j) {
				int32_t u = context_vars[j] ;
				abs_context_num_configs *= _Problem->K(u) ;
				abs_context.push_back(u) ;
				}

			int32_t nAbsContext = abs_context.size();
			if (nAbsContext > _maxContextSize)
				_maxContextSize = nAbsContext;
			if (_max_abs_context_num_configs < abs_context_num_configs) 
				_max_abs_context_num_configs = abs_context_num_configs ;
			}

		// collect/print stats for debug purposes
		std::vector<int32_t> context_size_list ; context_size_list.reserve(_Problem->N()+1) ; context_size_list.resize(_Problem->N()+1, 0) ;
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
			std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
			int32_t context_size = abs_context.size() ;
			++context_size_list[context_size] ;
			}
/*		printf("\nDEBUG : context_size_list (size, num) :") ;
		for (int32_t i = 0 ; i <= _Problem->N() ; ++i) {
			int32_t num = context_size_list[i] ;
			if (0 == num) continue ;
			printf("\n   %d %d", i, num) ;
			}*/

		return 0 ;
	}

	int32_t ComputeAbstractionContext(int n)
	// n is number of variables from union(given_var,context(given_var)) to include in the abs context
	{
		_MapVar2AbstractionContext.clear() ;
		_MapVar2AbstractionContext.reserve(_Problem->N()) ;
		if (_MapVar2AbstractionContext.capacity() < _Problem->N()) 
			return 1 ;
		_MapVar2AbstractionContext.resize(_Problem->N()) ;
		if (_Problem->N() != _MapVar2AbstractionContext.size()) 
			return 2 ;

		if (n <= 0) 
			return 0 ;

#ifdef ABS_CONTEXT_IS_EVERYTHING
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
			std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
			abs_context.clear() ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			for (; NULL != b ; b = b->ParentBucket()) {
				for (int32_t j = 0 ; j < b->nVars() ; j++) {
					int32_t u = b->Var(j) ;
					abs_context.push_back(u) ;
					}
				}
			}
#endif 

		int32_t left[32], right[32] ;
		for (int32_t i = _Problem->N()-1 ; i >= 0 ; --i) {
			int32_t v = i ;
//			int32_t v = (_Problem->VarOrdering_VarList())[i] ;
			std::vector<int32_t> & abs_context = _MapVar2AbstractionContext[v] ;
			abs_context.clear() ;
			abs_context.reserve(n) ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;
			if (NULL == b) continue ;
			// add this var to context
			abs_context.push_back(v) ;
			// collect another n-1 variables up the pseudotree
			int n_ = n - 1 ;
			BucketElimination::Bucket *b_ = b->ParentBucket() ;
			for (; n_ > 0 && NULL != b_ ; --n_, b_ = b_->ParentBucket()) {
				abs_context.push_back(b_->V()) ;
				}
			}
		return 0 ;
	}

protected :

	// current context, used to compute h values of AND/OR search space nodes.
	int32_t *_CurrentContextValues ;

public :

	inline int32_t *CurrentContextValues(void) { return _CurrentContextValues ; }
	inline void SetCurrentContextValue(int32_t Var, int32_t Value) { _CurrentContextValues[Var] = Value ; }

	int32_t SerializeCurrentContextValues(std::string & S)
	{
		char s[32] ;
		S.clear() ;
		for (int32_t i = 0 ; i < _Problem->N() ; i++) {
			if (_CurrentContextValues[i] >= 0) {
				if (S.length() > 0) S += ';' ;
				sprintf(s, "%d=%d", i, _CurrentContextValues[i]) ;
				S += s ;
				}
			}
		return 0 ;
	}

	inline int32_t GenerateCurrentContext(SearchNode & N, int32_t *context_values)
	{
		SearchNode *n = &N ;
		do {
			SearchAndNode *an = dynamic_cast<SearchAndNode*>(n) ;
			if (NULL != an ? an->V() >= 0 : false) 
				context_values[an->V()] = an->Value() ;
			n = n->Parent() ;
		}
		while (NULL != n) ;
		return 0 ;
	}

protected :

	MTRand _RNG ;

protected :

	// open list has a sorted (total order) nodes in increasing abstraction.
	// order is imposed by user-provided comparison function.
	// for proper abstraction : key is vector <var-pos-in-bt-order,ptr-to-closest-branching-node-AND-ansector>
	// for unique abstraction : key is vector <var-pos-in-bt-order,ptr-to-node>
	CMauiAVLTreeObj _OpenList ;

public :

	int32_t QueueExpansionNode(SearchAndNode & Node)
	{
		int32_t v = Node.V() ;
		if (0 == _OpenList.Insert(&Node)) 
			return 1 ;
		Node.NextInExpansionQueue() = NULL ;
		Node.ExpansionQueueKey() = _I64_MAX ;
		return 0 ;
	}

	inline SearchAndNode *PopNextNodeToExpandFromQueue(void) 
	{
		SearchAndNode *n = NULL ;
		_OpenList.RemoveFirst((void**)&n) ;
		if (NULL == n) 
			return NULL ;
		n->ExpansionQueueKey() = 0 ; // this node no longer in the queue
		return n ;
	}

public :

	int32_t PrepSampledSearchTree(int nAbsContext, bool ignoreclosestbranchingvariable)
	{
		if (0 != ComputeFunctionArgumentsPermutationList()) 
			return 1 ;
		if (0 != ComputeSubtreePartitioning()) 
			return 2 ;
		if (0 != ComputeClosestBranchingVariables()) 
			return 3 ;
		if (0 != ComputeAbstractionContextReal(nAbsContext, ignoreclosestbranchingvariable)) 
			return 4 ;
		return 0 ;	
	}

	// this fn is used by DFS traversal
	int32_t ComputeFrontierValue(/* BucketElimination::Bucket *B, */ bool include_h, std::vector<SearchAndNode_WithPath*> & Frontier, ARE_Function_TableType & V)
	{
		V = ARP_nInfinity ;
		if (0 == Frontier.size()) 
			return 0 ;
		for (int32_t j = 0 ; j < Frontier.size() ; ++j) {
			SearchAndNode_WithPath *a = Frontier[j] ;
			a->SubTreeValue() = a->AccumulatedISwSinceLastBranchingVariable() ;
			ApplyFnCombinationOperator(a->SubTreeValue(), a->AccumulatedCostSinceLastBranchingVariable()) ;
			if (include_h) 
//				{ a->ApplyPathAssignment() ;
//				ARE_Function_TableType v = ARP_DBL_MAX ;
//				B->ComputeHeuristic(_CurrentContextValues, v) ;
				ApplyFnCombinationOperator(a->SubTreeValue(), a->h()) ;
//				}
			if (ARP_DBL_MAX == V) 
				V = a->SubTreeValue() ;
			else 
				LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(V, V, a->SubTreeValue()) ;
			}
		return 0 ;
	}

	// do IS merge of two AND nodes [idxS,idxE)
	SearchAndNode_WithPath *DoISmerge(std::vector<SearchAndNode_WithPath*> & openList, int32_t idxS, int32_t idxE) ;

	// do IS merge of two AND nodes
	SearchAndNode *DoISmerge(SearchAndNode *A1, SearchAndNode *A2)
	{
		++_nNodesISmerged ;
		double qThis = A2->ISq() ;
		double qExisting = A1->ISq() ;
		// decide which one to keep; with probablity qThis/(qThis + qExisting) keep new and discard old.
		double sum = ARP_DBL_MAX, p ;
		if (ARP_nInfinity == qThis && ARP_nInfinity == qExisting) {
			p = 0.5 ;
			sum = ARP_nInfinity ;
			}
		else {
			LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(sum,qThis,qExisting)
			p = pow(10.0, qThis - sum) ;
			}
		double r = _RNG.rand() ;
		bool keep_existing = r >= p ;
		if (keep_existing) {
			// keep existing; update its weights.
			double p_ = 1.0 - p ; double log_p_ = log10(p_) ;
			A1->ISw() = _Problem->FunctionsAreConvertedToLogScale() ? A1->ISw() - log_p_ : A1->ISw()/p_ ;
			A1->AccumulatedISw() = _Problem->FunctionsAreConvertedToLogScale() ? A1->AccumulatedISw() - log_p_ : A1->AccumulatedISw()/p_ ;
			A1->ISq() = _Problem->FunctionsAreConvertedToLogScale() ? A1->ISq() - log_p_ : A1->ISq()/p_ ;
			A1->AccumulatedISwSinceLastBranchingVariable() = _Problem->FunctionsAreConvertedToLogScale() ? A1->AccumulatedISwSinceLastBranchingVariable() - log_p_ : A1->AccumulatedISwSinceLastBranchingVariable()/p_ ;
			delete A2 ;
			return A1 ;
			}
		else {
			// keep new; update its weights.
			double log_p = log10(p) ;
			A2->ISw() = _Problem->FunctionsAreConvertedToLogScale() ? A2->ISw() - log_p : A2->ISw()/p ;
			A2->AccumulatedISw() = _Problem->FunctionsAreConvertedToLogScale() ? A2->AccumulatedISw() - log_p : A2->AccumulatedISw()/p ;
			A2->ISq() = _Problem->FunctionsAreConvertedToLogScale() ? A2->ISq() - log_p : A2->ISq()/p ;
			A2->AccumulatedISwSinceLastBranchingVariable() = _Problem->FunctionsAreConvertedToLogScale() ? A2->AccumulatedISwSinceLastBranchingVariable() - log_p : A2->AccumulatedISwSinceLastBranchingVariable()/p ;
			delete A1 ;
			return A2 ;
			}
		return NULL ;
	}

	// this fn is used by DFS traversal
	int32_t ComputeAndChildren(
		// IN : closest ancestor branching AND node
		SearchAndNode_WithPath *ClosestAncestorAndNode, int32_t ClosestAncestorAndNodeChildIdx, 
		// IN : AND node for whom we want OR(AND) children generated
		SearchAndNode_WithPath *nParent, 
		// IN : bucket (var) that is a child of the given AND var, for which child OR(AND) nodes should be generated
		BucketElimination::Bucket *bChild, 
		// IN : if a beginning of a new DFS branch
		bool IsNewBeginning, 
		// IN : set values (in CurrentContextValues[]) of the entire path from this var to the root
		bool SetFullPathAssignment, 
		// OUT
		std::vector<SearchAndNode_WithPath*> & openListChild) ;
	// this fn is mostly the same as ComputeAndChildren(), except it accumulates full path assignment.
	int32_t ComputeAndChildrenEx(
		// IN : closest ancestor branching AND node
		SearchAndNode_WithPath *ClosestAncestorAndNode, int32_t ClosestAncestorAndNodeChildIdx, 
		// IN : AND node for whom we want OR(AND) children generated
		SearchAndNode_WithPath *nParent, 
		// IN : bucket (var) that is a child of the given AND var, for which child OR(AND) nodes should be generated
		BucketElimination::Bucket *bChild, 
		// IN : if a beginning of a new DFS branch
		bool IsNewBeginning, 
		// OUT
		std::vector<SearchAndNode_WithPath*> & openListChild) ;

	// generate AS tree and computeZ-estimate using DFS.
	// estimate (value) is stored in the _Root node.
	int32_t RunDFS(void) ;
	int32_t RunDFS_nonProper(int32_t & dDeepestBucketExpanded, int32_t & nDFSBranchingPointsProcessed) ;

	int32_t GenerateSampledSearchTree(void)
	{
		_Root.AccumulatedCost() = _Root.Cost() = _Root.h() = _Root.ISw() = _Root.AccumulatedISw() = FnCombinationNeutralValue() ;
		QueueExpansionNode(_Root) ;

		CMauiAVL64Tree orNodesToCheck ;

#ifdef DEBUG_AND_OR_CONSTRUCTION
		CMauiAVL64Tree nodes_expanded ;
		int32_t idxLastVar = -1 ;
#endif
		// expand queue until empty
		std::vector<BucketElimination::Bucket *> children_as_buckets ;
		while (true) {
			SearchAndNode *n = PopNextNodeToExpandFromQueue() ;

			if (NULL == n) 
				break ;
			int32_t v = n->V() ;
			BucketElimination::Bucket *b = MapVar2Bucket(v) ;

#ifdef DEBUG_AND_OR_CONSTRUCTION
			int32_t idxThisVar = v >= 0 ? (_Problem->VarOrdering_VarPos())[v] : -1 ;
			if (idxThisVar < idxLastVar) {
				printf("\n   Processing AND var %d (value=%d) nOpenList=%lld : error 1", v, n->Value(), (int64_t) _OpenList.GetSize()) ;
				exit(111) ;
				}
/*			int32_t r = nodes_expanded.Find((int64_t)n, NULL) ;
			if (0 != r) {
				printf("\n   Processing AND var %d (value=%d) nOpenList=%lld : error : this node already expanded", v, n->Value(), (int64_t) _OpenList.GetSize()) ;
				exit(112) ;
				}*/
//printf("\n   Processing AND var %d (value=%d) nOpenList=%lld", v, n->Value(), (int64_t) _OpenList.GetSize()) ;
			idxLastVar = idxThisVar ;
			nodes_expanded.Insert((int64_t)n, NULL) ;
#endif
			if (NULL != b) {
				if (EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(v) : false) {
					// if subtree of the bucket (pseudo) tree has no partitioning at this var, we don't need to expand it any more here. leave as is.
					continue ;
					}
				}
			// get child buckets in the pseudo tree
			children_as_buckets.clear() ;
			if (v < 0) {
				GetRoots(children_as_buckets) ;
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n nRoots=%d", children_as_buckets.size()) ;
#endif
				}
			else {
				if (NULL != b) b->GetChildren(children_as_buckets) ;
				}

			// generate current context; child AND nodes need it; do it here instead of while @ AND nodes.
			GenerateCurrentContext(*n, _CurrentContextValues) ;

			// generate OR children; for each generate AND children.
			for (int32_t iC = 0 ; iC < children_as_buckets.size() ; iC++) {
				BucketElimination::Bucket *bC = children_as_buckets[iC] ; if (NULL == bC) continue ;
				int32_t child = bC->V() ;
				BucketElimination::Bucket *bChild = MapVar2Bucket(child) ; if (NULL == bChild) continue ;
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n   Adding OR parent var %d (value=%d) this %d", v, n->Value(), child) ;
#endif
				SearchOrNode *onNew = new SearchOrNode(this, child) ;
				if (NULL == onNew) 
					return 1 ;
				 ++_nNodesCreated ;
				onNew->SetParent(n, -1, false) ;
				n->AddChild(*onNew) ;
				++_nNodesInTree ;
				// fetch range of values for the child
				int32_t idxS = 0, idxE = _Problem->K(child) ;
				if (_Problem->Value(child) >= 0) 
					{ idxS = _Problem->Value(v) ; idxE = idxS+1 ; }
				// generate AND children.
				for (int32_t iV = idxS ; iV < idxE ; iV++) {
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n   Adding AND var %d value %d", child, iV) ;
#endif
					SearchAndNode *anNew = new SearchAndNode(this, child, iV) ;
					if (NULL == anNew) 
						return 1 ;
					++_nNodesCreated ;
					anNew->SetParent(onNew, -1, false) ;
					anNew->SetClosestBranchingAncestor(n, -1) ;
					// add this AND nodes value to context
					_CurrentContextValues[child] = iV ;
					// compute label; accumulated cost(label)
					bChild->ComputeCost(_CurrentContextValues, anNew->Cost()) ;
					bChild->ComputeHeuristic(_CurrentContextValues, anNew->h()) ;
/*					if (EXPLOIT_EXACT_H ? ! SubtreeHasPartitioning(child) : false) {
						onNew->AddChild(*anNew) ;
						++_nNodesInTree ;
						continue ;
						}*/
					anNew->AccumulatedCost() = n->AccumulatedCost() ;
					ApplyFnCombinationOperator(anNew->AccumulatedCost(), anNew->Cost()) ;
					// for newly generated AND node, its own IS weight is 1 (normal scale, not log-scale)
					anNew->AccumulatedISw() = n->AccumulatedISw() ;
					// compute IS q value for this node
					anNew->ISq() = anNew->AccumulatedCost() ;
					ApplyFnCombinationOperator(anNew->ISq(), anNew->h()) ;
					ApplyFnCombinationOperator(anNew->ISq(), anNew->AccumulatedISw()) ;
					// check if an AND node with the abstraction already exists; if no, add this node; if yes, decide if need to merge.
					void *anExisting_ = NULL ;
					_OpenList.Find(anNew, &anExisting_) ;
					SearchAndNode *anExisting = (SearchAndNode*) anExisting_ ;
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n    andNew g=%g G=%g h=%g ISW=%g", anNew->Cost(), anNew->AccumulatedCost(), anNew->h(), anNew->AccumulatedISw()) ;
#endif
					if (NULL != anExisting_) {
						++_nNodesISmerged ;
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n    AND var=%d value=%d has existing", (int) child, (int) anNew->Value()) ;
#endif
						double qThis = anNew->ISq() ;
						double qExisting = anExisting->ISq() ;
						// decide which one to keep; with probablity qThis/(qThis + qExisting) keep new and discard old.
						double sum = ARP_DBL_MAX, p ;
double inf__ = -std::numeric_limits<double>::infinity() ;
						if (ARP_nInfinity == qThis && ARP_nInfinity == qExisting) {
							p = 0.5 ;
							sum = ARP_nInfinity ;
							}
						else {
							LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(sum,qThis,qExisting)
							p = pow(10.0, qThis - sum) ;
							}
						double r = _RNG.rand() ;
						bool keep_existing = r >= p ;
/*						if (_Problem->FunctionsAreConvertedToLogScale()) {
							qThis = pow(10.0, qThis) ;
							qExisting = pow(10.0, qExisting) ;
						}
						double sum = qThis + qExisting ;
						if (sum < 1.0e-128) 
							{ qThis = qExisting = 0.5 ; sum = 1.0 ; }
						double p = qThis/sum ;
						double r = _RNG.randExc(sum) ;
						bool keep_existing = r >= qThis ;*/

						if (keep_existing) {
							// keep existing; update its weights.
							double p_ = 1.0 - p ; double log_p_ = log10(p_) ;
double xyz = _Problem->FunctionsAreConvertedToLogScale() ? anExisting->ISw() - log_p_ : anExisting->ISw()/p_ ;
if (std::isnan(xyz)) {
	int bug = 1 ;
	}
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n    q=(%g,%g) p=%g keep existing wnew=%g/%g", qThis, qExisting, p, anExisting->ISw(), log_p_) ;
#endif
							anExisting->ISw() = _Problem->FunctionsAreConvertedToLogScale() ? anExisting->ISw() - log_p_ : anExisting->ISw()/p_ ;
							anExisting->AccumulatedISw() = _Problem->FunctionsAreConvertedToLogScale() ? anExisting->AccumulatedISw() - log_p_ : anExisting->AccumulatedISw()/p_ ;
							anExisting->ISq() = _Problem->FunctionsAreConvertedToLogScale() ? anExisting->ISq() - log_p_ : anExisting->ISq()/p_ ;
							delete anNew ; anNew = NULL ;
							// NOTE : orNew may have no children; but we will check and clean up at the end of the AND for-loop.
							++anExisting->nISmerges() ;
							}
						else {
							// keep new; update its weights.
							double log_p = log10(p) ;
double xyz = _Problem->FunctionsAreConvertedToLogScale() ? anNew->ISw() - log_p : anNew->ISw()/p ;
if (std::isnan(xyz)) {
	int bug = 1 ;
	}
#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n    q=(%g,%g) p=%g keep new      wnew=%g/%g", qThis, qExisting, p, anNew->ISw(), log_p) ;
#endif
							anNew->ISw() = _Problem->FunctionsAreConvertedToLogScale() ? anNew->ISw() - log_p : anNew->ISw()/p ;
							anNew->AccumulatedISw() = _Problem->FunctionsAreConvertedToLogScale() ? anNew->AccumulatedISw() - log_p : anNew->AccumulatedISw()/p ;
							anNew->ISq() = _Problem->FunctionsAreConvertedToLogScale() ? anNew->ISq() - log_p : anNew->ISq()/p ;
							// remove existing from open list
							int32_t resRemove = _OpenList.Remove(anExisting) ;
							// remove existing node from parent(its parent may have no other children; need to clean up);
							SearchOrNode *or_ = (SearchOrNode *) anExisting->Parent() ;
							or_->RemoveChild(*anExisting) ;
							--_nNodesInTree ;
							if (0 == or_->nChildren()) 
								orNodesToCheck.Insert((int64_t) or_, NULL) ;
							delete anExisting ; anExisting_ = NULL ;
							++anNew->nISmerges() ;
							}
						}
					else {
						// new AND node
					}

					if (NULL != anNew) {
						onNew->AddChild(*anNew) ;
						++_nNodesInTree ;
						int32_t resQue = QueueExpansionNode(*anNew) ;
						if (0 != resQue) 
							{ onNew->RemoveChild(*anNew) ; delete anNew ; return 1 ; }
						}
					}
					// orNew may have no children; need to clean up in this case.
					if (0 == onNew->nChildren()) 
						orNodesToCheck.Insert((int64_t) onNew, NULL) ;
					}
				}
			// clean up; check if affected OR nodes have any children; if no, do cleanup.
			SearchOrNode *or__ = NULL ;
			while (orNodesToCheck.RemoveFirst((int64_t*) &or__, NULL) > 0) {
				if (0 == or__->nChildren()) {
					SearchAndNode *an = (SearchAndNode *) or__->Parent() ;
					an->RemoveChild(*or__) ;
					--_nNodesInTree ;
					delete or__ ;
					// NOTE : an should be an expanded AND node, i.e. not on open list; if it has no children, remove it.
					if (0 == an->nChildren() ? &_Root != an : false) {
						SearchOrNode *or_ = (SearchOrNode *) an->Parent() ;
						or_->RemoveChild(*an) ;
						--_nNodesInTree ;
						if (0 == or_->nChildren()) 
							orNodesToCheck.Insert((int64_t) or_, NULL) ;
						// NOTE : an should not be in expansion queue
						delete an ;
					}
				}
			}

//if (0 == _nNodesInTree%1000000) 
#ifdef DEBUG_AND_OR_CONSTRUCTION
//printf("\nnodes current summary %d", _nNodesInTree) ;
#endif

#ifdef DEBUG_AND_OR_CONSTRUCTION
printf("\n   nodes summary nInTree=%d nOpenList=%d nNodesCreated=%d nNodesISmerged=%d", (int) _nNodesInTree, (int) _OpenList.GetSize(), (int) _nNodesCreated, (int) _nNodesISmerged) ;
#endif
		return 0 ;
		}

	int32_t ComputeSearchTreeValue(ARE_Function_TableType & TreeValue)
	{
		int32_t res = 1 ;

		TreeValue = _AnswerFactor ;

		for (int i = 0 ; i < _Problem->N() ; ++i) 
			_CurrentContextValues[i] = -1 ;

		int32_t DFS_height = 1 + 2*(1+MaxTreeHeight()) ;
		SearchNode **stack_node = new SearchNode*[DFS_height] ;
		SearchNode **stack_next_child2process = new SearchNode*[DFS_height] ;
		if (NULL == stack_node || NULL == stack_next_child2process) 
			goto done ;
		int32_t stack_depth ; stack_depth = 0 ;
		stack_node[0] = &_Root ;
		stack_next_child2process[0] = _Root.Children() ;
		
		while (stack_depth >= 0) {
			if (NULL == stack_next_child2process[stack_depth]) {
				SearchOrNode *on = dynamic_cast<SearchOrNode*>(stack_node[stack_depth]) ;
				SearchAndNode *an = dynamic_cast<SearchAndNode*>(stack_node[stack_depth]) ;
				if (NULL != on) {
					SearchNode *nC = on->Children() ;
					on->SubTreeValue() = nC->SubTreeValue() ;
#ifdef DEBUG_AND_OR_COMPUTATION
					printf("\nComputing OR (var=%d) child value=%g", on->V(), nC->SubTreeValue()) ;
#endif
					for (nC = nC->Sibling() ; NULL != nC ; nC = nC->Sibling()) {
						LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(on->SubTreeValue(), on->SubTreeValue(), nC->SubTreeValue()) ;
#ifdef DEBUG_AND_OR_COMPUTATION
						printf("\nComputing OR (var=%d) child value=%g", on->V(), nC->SubTreeValue()) ;
#endif
						}
#ifdef DEBUG_AND_OR_COMPUTATION
					std::string s ;
					SerializeCurrentContextValues(s) ;
					printf("\nComputing OR (var=%d) context={%s} final value=%g", on->V(), s.c_str(), on->SubTreeValue()) ;
#endif
					}
				else if (NULL != an) {
					BucketElimination::Bucket *b = MapVar2Bucket(an->V()) ;
					// start with IS weight of this node; could be 1 (or 0 in log scale).
					an->SubTreeValue() = an->ISw() ;
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing AND (var=%d;value=%d) ISw=%g (v so far %g)", an->V(), an->Value(), an->ISw(), an->SubTreeValue()) ;
#endif
					// fetch cost of 'b'
					if (NULL != b) {
						// 2016-11-10 KK : we assume cost is computed.
						ARE_Function_TableType v = an->Cost() ;
						// b->ComputeCost(_CurrentContextValues, v) ;
						ApplyFnCombinationOperator(an->SubTreeValue(), v) ;
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing AND (var=%d;value=%d) Cost=%g (v so far %g)", an->V(), an->Value(), v, an->SubTreeValue()) ;
#endif
						if (0 == an->nChildren()) {
							// if AND node has no children, it may have h value defined (e.g. when this subtree has no partitioning); include this value.
							b->ComputeHeuristic(_CurrentContextValues, v) ;
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing AND (var=%d;value=%d) h=%g nBfns=%d", an->V(), an->Value(), v, b->nOriginalFunctions()) ;
#endif
							ApplyFnCombinationOperator(an->SubTreeValue(), v) ;
							}
						}
					// add values of children
					for (SearchNode *nC = an->Children() ; NULL != nC ; nC = nC->Sibling()) {
						ApplyFnCombinationOperator(an->SubTreeValue(), nC->SubTreeValue()) ;
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing AND (var=%d;value=%d) child_value=%g (v so far %g)", an->V(), an->Value(), nC->SubTreeValue(), an->SubTreeValue()) ;
#endif
						}
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing AND (var=%d;value=%d) final_value=%g", an->V(), an->Value(), an->SubTreeValue()) ;
					std::string s ;
					SerializeCurrentContextValues(s) ;
					printf("\nComputing AND (var=%d;value=%d) context={%s} nBfns=%d value=%g", an->V(), an->Value(), s.c_str(), (int) (NULL != b) ? b->nOriginalFunctions() : -1, an->SubTreeValue()) ;
#endif
#ifdef DEBUG_AND_OR_COMPUTATION
					// DEBUGGG; if h is exact, compute value
					if (an->V() >= 0 ? ! SubtreeHasPartitioning(an->V()) : false) {
						ARE_Function_TableType local_value = an->Cost() ;
						ApplyFnCombinationOperator(local_value, an->h()) ;
						std::string s ;
						SerializeCurrentContextValues(s) ;
						printf("\n  Computing *** AND (var=%d;value=%d) ISvalue=%g LOCALvalue=%g(%g + %g) nIS=%d context={%s}", an->V(), an->Value(), an->SubTreeValue(), local_value, an->Cost(), an->h(), an->nISmerges(), s.c_str()) ;
						}
#endif
					// erase current context value
					if (an->V() >= 0) 
						_CurrentContextValues[an->V()] = -1 ;
					}
				// backtrack
				--stack_depth ;
				}
			else {
#ifdef DEBUG_AND_OR_COMPUTATION
SearchNode *nP = stack_next_child2process[stack_depth] ;
#endif
				stack_node[stack_depth+1] = stack_next_child2process[stack_depth] ;
				stack_next_child2process[stack_depth] = stack_next_child2process[stack_depth]->Sibling() ;
				stack_next_child2process[++stack_depth] = stack_node[stack_depth]->Children() ;
				// fill in context values array
				SearchAndNode *an = dynamic_cast<SearchAndNode*>(stack_node[stack_depth]) ;
				if (NULL != an) {
					_CurrentContextValues[an->V()] = an->Value() ;
					}
#ifdef DEBUG_AND_OR_COMPUTATION
SearchOrNode *onP = dynamic_cast<SearchOrNode*>(nP) ;
SearchAndNode *anP = dynamic_cast<SearchAndNode*>(nP) ;
if (NULL != onP) {
	printf("\n  DFS expand : parent OR %d; children :",nP->V()) ;
	for (SearchNode *n_ = stack_next_child2process[stack_depth] ; NULL != n_ ; n_ = n_->Sibling()) {
		SearchAndNode *an_ = (SearchAndNode *) n_ ;
		printf(" %d=%d",an_->V(),an_->Value()) ;
		}
	}
else {
	printf("\n  DFS expand : parent AND %d=%d; children :",nP->V(), anP->Value()) ;
	for (SearchNode *n_ = stack_next_child2process[stack_depth] ; NULL != n_ ; n_ = n_->Sibling()) {
		SearchOrNode *on_ = (SearchOrNode *) n_ ;
		printf(" %d",on_->V()) ;
		}
	}
#endif
				}
			}
#ifdef DEBUG_AND_OR_COMPUTATION
printf("\n  Computing root AND AnswerFactor()=%g", AnswerFactor()) ;
#endif
		ApplyFnCombinationOperator(_Root.SubTreeValue(), AnswerFactor()) ;
		res = 0 ;

		// DEBUGGGGGING : log root values from sampling and MBE tree (forest)
		/*
		{
		std::vector<BucketElimination::Bucket *> children_as_buckets ;
		GetRoots(children_as_buckets) ;
		//printf("\n-----ABS-ROOTS-----") ;
		for (SearchNode *n_ = _Root.Children() ; NULL != n_ ; n_ = n_->Sibling()) {
			printf("\nv=%d value=%g", n_->V(), n_->SubTreeValue()) ;
			}
		//printf("\n-----MBE-ROOTS-----") ;
		for (int32_t iC = 0 ; iC < children_as_buckets.size() ; iC++) {
			BucketElimination::Bucket *bC = children_as_buckets[iC] ; if (NULL == bC) continue ;
			if (bC->nMiniBuckets() <= 0) continue ;
			printf("\nv=%d value=%g (nSTbuckets=%d,nChildren=%d)", bC->V(), bC->MiniBuckets()[0]->OutputFunction().ConstValue(), bC->nSubtreeBuckets(), bC->nChildren()) ;
			}
		}
		*/

done :
		if (NULL != stack_node) 
			delete stack_node ;
		if (NULL != stack_next_child2process) 
			delete stack_next_child2process ;
		return res ;
	}

	int32_t ComputeTreeSize(void)
	{
		int32_t n = 1 ;

		int32_t DFS_height = 1 + 2*(1+MaxTreeHeight()) ;
		SearchNode **stack_node = new SearchNode*[DFS_height] ;
		SearchNode **stack_next_child2process = new SearchNode*[DFS_height] ;
		if (NULL == stack_node || NULL == stack_next_child2process) 
			goto done ;
		int32_t stack_depth ; stack_depth = 0 ;
		stack_node[0] = &_Root ;
		stack_next_child2process[0] = _Root.Children() ;
		
		while (stack_depth >= 0) {
			if (NULL == stack_next_child2process[stack_depth]) {
				--stack_depth ;
				}
			else {
				++n ;
				stack_node[stack_depth+1] = stack_next_child2process[stack_depth] ;
				stack_next_child2process[stack_depth] = stack_next_child2process[stack_depth]->Sibling() ;
				stack_next_child2process[++stack_depth] = stack_node[stack_depth]->Children() ;
				}
			}

done :
		if (NULL != stack_node) 
			delete stack_node ;
		if (NULL != stack_next_child2process) 
			delete stack_next_child2process ;
		return n ;
	}

public :

	virtual int32_t Initialize(ARE::ARP & Problem, bool UseLogScale, const int32_t *VarOrderingAsVarList, int32_t DeleteUsedTables)
	{
		int32_t res = BucketElimination::MBEworkspace::Initialize(Problem, UseLogScale, VarOrderingAsVarList, DeleteUsedTables) ;
		if (0 != res) 
			return res ;

		if (Problem.N() < 1) 
			return 0 ;

		if (_OpenList.GetSize() > 0) 
			_OpenList.Empty() ;

		_CurrentContextValues = new int32_t[Problem.N()] ;
		if (NULL == _CurrentContextValues) 
			return 1 ;
		for (int32_t i = 0 ; i < Problem.N() ; i++) _CurrentContextValues[i] = -1 ;

		return 0 ;
	}

	virtual int32_t CollectAndNodes(std::vector<SearchAndNode*> & Nodes)
	{
		int32_t res = 1 ;

		int32_t DFS_height = 1 + 2*(1+MaxTreeHeight()) ;
		SearchNode **stack_node = new SearchNode*[DFS_height] ;
		SearchNode **stack_next_child2process = new SearchNode*[DFS_height] ;
		if (NULL == stack_node || NULL == stack_next_child2process) 
			goto done ;
		int32_t stack_depth ; stack_depth = 0 ;
		stack_node[0] = &_Root ;
		stack_next_child2process[0] = _Root.Children() ;
		
		while (stack_depth >= 0) {
			if (NULL == stack_next_child2process[stack_depth]) {
				SearchNode *n = stack_node[stack_depth] ;
				SearchAndNode *an = dynamic_cast<SearchAndNode*>(n) ;
				if (NULL != an) 
					Nodes.push_back(an) ;
				// backtrack
				--stack_depth ;
				}
			else {
				stack_node[stack_depth+1] = stack_next_child2process[stack_depth] ;
				stack_next_child2process[stack_depth] = stack_next_child2process[stack_depth]->Sibling() ;
				stack_next_child2process[++stack_depth] = stack_node[stack_depth]->Children() ;
				}
			}

done :
		if (NULL != stack_node) 
			delete stack_node ;
		if (NULL != stack_next_child2process) 
			delete stack_next_child2process ;
		return res ;
	}

	virtual int32_t DestroyTree(void)
	{
		int32_t res = 1 ;

		int32_t DFS_height = 1 + 2*(1+MaxTreeHeight()) ;
		SearchNode **stack_node = new SearchNode*[DFS_height] ;
		SearchNode **stack_next_child2process = new SearchNode*[DFS_height] ;
		if (NULL == stack_node || NULL == stack_next_child2process) 
			goto done ;
		int32_t stack_depth ; stack_depth = 0 ;
		stack_node[0] = &_Root ;
		stack_next_child2process[0] = _Root.Children() ;
		
		while (stack_depth >= 0) {
			if (NULL == stack_next_child2process[stack_depth]) {
				SearchNode *n = stack_node[stack_depth] ;
				if (&_Root != n) 
					delete n ;
				// backtrack
				--stack_depth ;
				}
			else {
				stack_node[stack_depth+1] = stack_next_child2process[stack_depth] ;
				stack_next_child2process[stack_depth] = stack_next_child2process[stack_depth]->Sibling() ;
				stack_next_child2process[++stack_depth] = stack_node[stack_depth]->Children() ;
				}
			}

		_Root.AttachChildren(NULL, 0) ;
		_nNodesInTree = 1 ;
		_nNodesCreated = 1 ;
		_nNodesISmerged = 0 ;
		_nNodesExactHComputed = 0 ;
		_nExpansionCalls = 0 ;
		_nTotalSumFrontierSize = 0 ;
		_nNumVarExpansions = 0 ;
		if (_OpenList.GetSize() > 0) 
			_OpenList.Empty() ;

done :
		if (NULL != stack_node) 
			delete stack_node ;
		if (NULL != stack_next_child2process) 
			delete stack_next_child2process ;
		return res ;
	}

	virtual int32_t Destroy(void)
	{
		DestroyTree() ;

		if (NULL != _CurrentContextValues) 
			{ delete [] _CurrentContextValues ; _CurrentContextValues = NULL ; }

/*		int N = NULL != _Problem ? _Problem->N() : _BFS_expansion_queue_size.size() ;
		for (int i = 0 ; i < N ; i++) {
			int64_t & size = _BFS_expansion_queue_size[i] ;
			SearchAndNode * & H = _BFS_expansion_queue_H[i] ;
			SearchAndNode * & T = _BFS_expansion_queue_T[i] ;
			while (NULL != H) {
				SearchNode *n = H ;
				H = (SearchAndNode*) H->NextInExpansionQueue() ;
				n->NextInExpansionQueue() = NULL ;
				delete n ;
				}
			H = T = NULL ; size = 0 ;
			}*/

		_nVarsPartitioned = -1 ;

		_nRandAbs = 0 ;
		_RandAbsFactors.clear() ;
		_RandAbsFactorPerVar.clear() ;

		BucketElimination::MBEworkspace::Destroy() ;
		return 0 ;
	}

	AbsSamplingWorkspace(AbsSamplingCompFn CompFn) :
		_Root(this, -1, -1), 
		_CompFn(CompFn), 
		_nNodesCreatedLimitBeforeLevellingOff(INT64_MAX), 
		_nNodesInTree(1), 
		_nNodesCreated(1), 
		_nNodesISmerged(0), 
		_nNodesExactHComputed(0), 
		_nExpansionCalls(-1), 
		_nTotalSumFrontierSize(-1), 
		_nNumVarExpansions(-1), 
		_nVarsPartitioned(-1), 
		_maxContextSize(-1),
		_max_abs_context_num_configs(-1), 
		_PTInducedWidth(-1),
		_nRandAbs(0), 
		_CurrentContextValues(NULL), 
		_OpenList(CompFn, 0, NULL, 100000, 10000) 
	{
	}

	virtual ~AbsSamplingWorkspace(void)
	{
		Destroy() ;
	}
} ;

}

#endif // AbsSamplingWorkspace_HXX_INCLUDED