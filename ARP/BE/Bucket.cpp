#include <stdlib.h>
#include <memory.h>

#include <Function.hxx>
#include <Bucket.hxx>
#include <MBEworkspace.hxx>

#define DO_PROD_SUM_MM

// 2018-01-13 KK : not recomputing Intermediate FNs is not working correctly, I think.
#define RECOMPUTE_INT_FNS_in_H_COMPUTATION

BucketElimination::Bucket::Bucket(void)
	:
	_Workspace(NULL),
	_IDX(-1), 
	_V(-1), 
	_PostComputeExtNote(NULL), 
	_Width(-1), 
	_Signature(NULL), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_ParentBucket(NULL), 
	_RootBucket(NULL), 
	_DistanceToRoot(-1), 
	_Height(-1), 
	_nSubtreeBuckets(-1), 
	_MaxDescendantNumVars(-1), 
	_MaxTreeHeight_BranchingVars(-1), 
	_ComputationNewFunctionSize(-1), 
	_MaxDescendantComputationNewFunctionSize(-1), 
	_nChildren(-1), 
	_ChildVars(NULL), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nAugmentedFunctions(0), 
	_AugmentedFunctions(NULL), 
	_AugmentedFunctionsArraySize(0), 
	_nIntermediateFunctions(0), 
	_IntermediateFunctions(NULL), 
	_IntermediateFunctionsArraySize(0), 
	_NextInOrder(NULL)
{
}


BucketElimination::Bucket::~Bucket(void)
{
	Destroy() ;
}


BucketElimination::Bucket::Bucket(MBEworkspace & WS, int32_t IDX, int32_t V)
	:
	_Workspace(&WS),
	_IDX(IDX), 
	_V(V),
	_PostComputeExtNote(NULL), 
	_Width(-1), 
	_Signature(NULL), 
	_nVars(0), 
	_Vars(NULL), 
	_VarsSpace(0), 
	_ParentBucket(NULL), 
	_RootBucket(NULL), 
	_DistanceToRoot(-1), 
	_Height(-1), 
	_nSubtreeBuckets(-1), 
	_MaxDescendantNumVars(-1), 
	_MaxTreeHeight_BranchingVars(-1), 
	_ComputationNewFunctionSize(-1), 
	_MaxDescendantComputationNewFunctionSize(-1), 
	_nChildren(-1), 
	_ChildVars(NULL), 
	_nOriginalFunctions(0), 
	_OriginalFunctions(NULL), 
	_OriginalWidth(-1), 
	_OriginalSignature(NULL), 
	_nAugmentedFunctions(0), 
	_AugmentedFunctions(NULL), 
	_AugmentedFunctionsArraySize(0), 
	_nIntermediateFunctions(0), 
	_IntermediateFunctions(NULL), 
	_IntermediateFunctionsArraySize(0), 
	_NextInOrder(NULL)
{
	// if Var is given, add it
	if (V >= 0) 
		AddVar(V) ;
}


void BucketElimination::Bucket::Destroy(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		mb->Destroy() ;
		delete mb ;
		}
	_MiniBuckets.clear() ;

	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	if (NULL != _AugmentedFunctions) {
		delete [] _AugmentedFunctions ;
		_AugmentedFunctions = NULL ;
		_AugmentedFunctionsArraySize = 0 ;
		}
	if (NULL != _IntermediateFunctions) {
		delete [] _IntermediateFunctions ;
		_IntermediateFunctions = NULL ;
		_IntermediateFunctionsArraySize = 0 ;
		}
	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_Width = -1 ;
	_nAugmentedFunctions = 0 ;
	_nIntermediateFunctions = 0 ;
	_OriginalWidth = -1 ;
	_nOriginalFunctions = 0 ;
	_nVars = 0 ;
	if (NULL != _Vars) {
		delete [] _Vars ;
		_Vars = NULL ;
		}
	_VarsSpace = 0 ;
	_MaxDescendantNumVars = -1 ;
	_MaxTreeHeight_BranchingVars = -1 ;
	_ComputationNewFunctionSize = -1 ;
	_MaxDescendantComputationNewFunctionSize = -1 ;
	_nChildren = -1 ;
	if (NULL != _ChildVars) { delete [] _ChildVars ; _ChildVars = NULL ; }
}


int32_t BucketElimination::Bucket::GetChildren(std::vector<BucketElimination::Bucket*> & Children)
{
	for (int32_t i = 0 ; i < _nChildren ; i++) {
		int32_t child = _ChildVars[i] ;
		Bucket *b = _Workspace->MapVar2Bucket(child) ;
		if (NULL != b) 
			Children.push_back(b) ;
		}
	return 0 ;
}


BucketElimination::Bucket *BucketElimination::Bucket::ChildBucket(int32_t idx) const
{
	return _Workspace->MapVar2Bucket(_ChildVars[idx]) ;
}


bool BucketElimination::Bucket::VerifyAugmentedFunctions(void)
{
	for (int32_t i = 0 ; i < _nAugmentedFunctions ; ++i) {
		ARE::Function *f = _AugmentedFunctions[i] ;
		int32_t j = 0 ;
		for (; j < f->N() ; ++j) {
			if (_V == f->Argument(j)) 
				break ;
			}
		if (j >= f->N()) 
			return 1 ;
		// check that this augmented fn is intermediate fn in all buckets between originating bucket and this bucket (exclusive)
		if (f->IDX() < 0) {
			int vOriginating = -(f->IDX() + 1) ;
			BucketElimination::Bucket *b = f->OriginatingBucket() ; // _Workspace->MapVar2Bucket(vOriginating) ;
			if (NULL == b) 
				return 1 ;
			for (BucketElimination::Bucket *b_ = b->ParentBucket() ; this != b_ ; b_ = b_->ParentBucket()) {
				if (NULL == b_) 
					return 1 ;
				int32_t res_idx = b_->ContainsIntermediateFunction(*f) ;
				if (res_idx < 0) 
					return 1 ;
				}
			}
		else {
			int32_t bug = 1 ;
			}
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::GetDFSorderofsubtree(
	// OUT
	std::vector<BucketElimination::Bucket*> & dfs_reverse_chain, 
	// IN : helper arrays
	std::vector<BucketElimination::Bucket*> & dfs_stack_b, std::vector<int32_t> & dfs_stack_i)
{
	dfs_stack_b.clear() ; dfs_stack_i.clear() ; dfs_stack_b.push_back(this) ; dfs_stack_i.push_back(0) ;
	while (dfs_stack_b.size() > 0) {
		// expand frontier
		int32_t depth = dfs_stack_b.size()-1 ;
		BucketElimination::Bucket *frontier = dfs_stack_b[depth] ;
		int32_t nProcessed = dfs_stack_i[depth] ;
		if (nProcessed >= frontier->nChildren()) {
			dfs_reverse_chain.push_back(frontier) ;
			dfs_stack_b.pop_back() ;
			dfs_stack_i.pop_back() ;
			if (--depth >= 0) dfs_stack_i[depth]++ ;
			}
		else {
			dfs_stack_b.push_back(_Workspace->MapVar2Bucket(frontier->ChildVar(nProcessed))) ;
			dfs_stack_i.push_back(0) ;
			}
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::GetDescendants(int32_t nLevelsDown, std::vector<BucketElimination::Bucket*> & Descendants)
{
	if (nLevelsDown <= 0) 
		return 0 ;
	int32_t j = Descendants.size() ;
	GetChildren(Descendants) ;
	int32_t n = 1 ; // 1 level is already done
	for (; n < nLevelsDown ; ++n) {
		int32_t m = Descendants.size() ;
		if (j >= m) 
			break ; // done; no more descendants.
		for (; j < m ; ++j) 
			Descendants[j]->GetChildren(Descendants) ;
		}
	return 0 ;
}

/*
int32_t BucketElimination::Bucket::SaveXMLString(const char *prefixspaces, const std::string & Dir, std::string & S)
{
	char s[1024] ;
	std::string temp ;
	int32_t i, j ;
	sprintf(s, "%s<bucket IDX=\"%d\" nVars=\"%d\" Vars=\"", prefixspaces, _IDX, _nVars) ;
	S += s ;
	for (i = 0 ; i < _nVars ; i++) {
		sprintf(s, "%d", _Vars[i]) ;
		if (i > 0) 
			S += ';' ;
		S += s ;
		}
	S += '"' ;
	if (NULL !=_ParentBucket) {
		sprintf(s, " parentbucket=\"%d\"", _ParentBucket->IDX()) ;
		S += s ;
		}
	S += ">" ;

	// save original functions
	if (_nOriginalFunctions > 0 && NULL != _OriginalFunctions) {
		sprintf(s, "\n%s <originalfunctions n=\"%d\" list=\"", prefixspaces, _nOriginalFunctions) ;
		S += s ;
		for (i = 0 ; i < _nOriginalFunctions ; i++) {
			ARE::Function *f = _OriginalFunctions[i] ;
			if (NULL == f) continue ;
			sprintf(s, "%d", f->IDX()) ;
			if (i > 0) 
				S += ';' ;
			S += s ;
			}
		S += "\"/>" ;
		}

	// save incoming child bucket functions
	if (_nChildBucketFunctions > 0 && NULL != _ChildBucketFunctions) {
		sprintf(s, "%s ", prefixspaces) ;
		for (i = 0 ; i < _nChildBucketFunctions ; i++) {
			ARE::Function *f = _ChildBucketFunctions[i] ;
			if (NULL == f) continue ;
			S += "\n" ;
			// bucketfunctions don't have IDX; set the IDX of the originating bucket as the IDX of this function
			int32_t idx = f->IDX() ;
			Bucket *b = f->OriginatingBucket() ;
			if (NULL != b) 
				f->SetIDX(b->IDX()) ;
			f->SaveXMLString(s, "childbucketfn", Dir, S) ;
			f->SetIDX(idx) ;
			}
		}

	// serialize _OutputFunction
	if (_OutputFunction.N() > 0) {
		sprintf(s, "%s ", prefixspaces) ;
		S += "\n" ;
		_OutputFunction.SaveXMLString(s, "ownbucketfn", Dir, S) ;
		}

	// serialize _OutputFunctionBlockComputationResult
	if (_OutputFunction.N() > 0 && NULL != _OutputFunctionBlockComputationResult) {
		__int64 nTB = _OutputFunction.nTableBlocks() ;
		int32_t size = (7 + nTB) >> 3 ;
		sprintf(s, "\n%s <ownbucketfncomputationresult nComputed=\"%I64d/%I64d\" bits=\"", prefixspaces, _nOutputFunctionBlocksComputed, nTB) ;
		S += s ;
		int32_t n = 0 ;
		for (i = 0 ; i < size && n < nTB ; i++) {
			for (j = 0 ; j < 8 && n < nTB ; j++, n++) {
				int32_t bit = _OutputFunctionBlockComputationResult[i] & (1 << j) ;
				if (0 != bit) 
					S += '1' ;
				else 
					S += '0' ;
				}
			}
		S += "\"/>" ;
		}

	sprintf(s, "\n%s</bucket>", prefixspaces) ;
	S += s ;

	return 0 ;
}
*/

int32_t BucketElimination::Bucket::SetOriginalFunctions(int32_t N, ARE::Function *FNs[]) 
{
	if (NULL != _OriginalFunctions) {
		delete [] _OriginalFunctions ;
		_OriginalFunctions = NULL ;
		}
	if (NULL != _OriginalSignature) {
		delete [] _OriginalSignature ;
		_OriginalSignature = NULL ;
		}
	_OriginalWidth = -1 ;
	_nOriginalFunctions = 0 ;
	return AddOriginalFunctions(N, FNs) ;
}


int32_t BucketElimination::Bucket::AddOriginalFunctions(int32_t N, ARE::Function *FNs[]) 
{
	if (N < 1) 
		return 0 ;

	int32_t i, j, k ;

	// check if functions in FNs are already in this bucket
	for (i = N-1 ; i >= 0 ; i--) {
		ARE::Function *f = FNs[i] ;
		for (j = 0 ; j < _nOriginalFunctions ; j++) {
			if (_OriginalFunctions[j] == f) {
				FNs[i] = FNs[--N] ;
				break ;
				}
			}
		}
	if (N < 1) 
		return 0 ;

	// reallocate original functions array
	int32_t space = _nOriginalFunctions + N ;
	ARE::Function **fnlist = new ARE::Function*[space] ;
	if (NULL == fnlist) 
		return 1 ;
	for (i = 0 ; i < _nOriginalFunctions ; i++) 
		fnlist[i] = _OriginalFunctions[i] ;
	for (; i < space ; i++) {
		fnlist[i] = FNs[i-_nOriginalFunctions] ;
		fnlist[i]->SetBucket(this) ;
		}
	delete [] _OriginalFunctions ;
	_OriginalFunctions = fnlist ;
	_nOriginalFunctions = space ;

	if (NULL != _OriginalSignature) { delete [] _OriginalSignature ; _OriginalSignature = NULL ; }
	_OriginalWidth = -1 ;

	if (0 != ARE::ComputeSignature(N, FNs, _OriginalWidth, _OriginalSignature)) 
		goto failed ;

	return 0 ;

failed :
	Destroy() ;
	return 1 ;
}


int32_t BucketElimination::Bucket::AddAugmentedFunction(ARE::Function & F)
{
	// check if we have enough space
	if (_nAugmentedFunctions+1 > _AugmentedFunctionsArraySize) {
		// 2016-05-06 KK : make it 16 for reallocations (used to be 8) so that there are fewer reallocations
		int32_t newsize = 0 == _AugmentedFunctionsArraySize ? 4 : _AugmentedFunctionsArraySize + 16 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nAugmentedFunctions > 0) 
			memcpy(newspace, _AugmentedFunctions, sizeof(ARE::Function *)*_nAugmentedFunctions) ;
		delete []_AugmentedFunctions ;
		_AugmentedFunctions = newspace ;
		_AugmentedFunctionsArraySize = newsize ;
		}

	_AugmentedFunctions[_nAugmentedFunctions++] = &F ;

/* 2016-12-06 KK : don't invalidate signature; sometimes we want to keep the previous signature (e.g. we run BE to compute BE signature, and when we run 
	MBE we don't want to change/delete BE signature; if we want to get rid of the signature, have calling fn do it explictly. */
/*
	// if fn has arguments, current signature (may be) is invalid
	if (F.N() > 0) {
		_Width = -1 ;
		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		}
*/

	return 0 ;
}


int32_t BucketElimination::Bucket::RemoveAugmentedFunction(ARE::Function & F, bool InvalidateSignature)
{
	int32_t i ;
	int32_t n = 0 ;
	for (i = _nAugmentedFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _AugmentedFunctions[i]) {
			_AugmentedFunctions[i] = _AugmentedFunctions[--_nAugmentedFunctions] ;
			F.SetBucket(NULL) ;
			n++ ;
			}
		}

	// if fn has arguments, current signature (may be) is invalid
	if (InvalidateSignature && n > 0 && F.N() > 0) {
		_Width = -1 ;
		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		}

	return 0 ;
}


int32_t BucketElimination::Bucket::AddIntermediateFunction(ARE::Function & F)
{
	// check if we have enough space
	if (_nIntermediateFunctions+1 > _IntermediateFunctionsArraySize) {
		// 2016-05-06 KK : make reallcation size larger (used to be 8) so that there are fewer reallocations
		int32_t newsize = 5 ; 
		if (_IntermediateFunctionsArraySize >= 65) 
			newsize = _IntermediateFunctionsArraySize + 60 ;
		else if (_IntermediateFunctionsArraySize >= 25) 
			newsize = _IntermediateFunctionsArraySize + 40 ;
		else if (_IntermediateFunctionsArraySize > 0) 
			newsize = _IntermediateFunctionsArraySize + 20 ;
//		int32_t newsize = (0 == _IntermediateFunctionsArraySize) ? 8 : _IntermediateFunctionsArraySize + 32 ;
		ARE::Function **newspace = new ARE::Function*[newsize] ;
		if (NULL == newspace) 
			return 1 ;
		if (_nIntermediateFunctions > 0) 
			memcpy(newspace, _IntermediateFunctions, sizeof(ARE::Function *)*_nIntermediateFunctions) ;
		delete []_IntermediateFunctions ;
		_IntermediateFunctions = newspace ;
		_IntermediateFunctionsArraySize = newsize ;
		}

	_IntermediateFunctions[_nIntermediateFunctions++] = &F ;

	return 0 ;
}


int32_t BucketElimination::Bucket::RemoveIntermediateFunction(ARE::Function & F)
{
	int32_t i ;
	for (i = _nIntermediateFunctions - 1 ; i >= 0 ; i--) {
		if (&F == _IntermediateFunctions[i]) {
			_IntermediateFunctions[i] = _IntermediateFunctions[--_nIntermediateFunctions] ;
			F.SetBucket(NULL) ;
			}
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeSignature(unsigned char instructions)
{
	bool IncludeIntermediateFunctions = instructions & 0x04 ;
	bool IncludeAugmentedFunctions = instructions & 0x02 ;
	bool IncludeOriginalFunctions = instructions & 0x01 ;

	if (NULL != _Signature) {
		delete [] _Signature ;
		_Signature = NULL ;
		}
	_Width = -1 ;

	if (_OriginalWidth < 0) {
		if (NULL != _OriginalSignature) {
			delete [] _OriginalSignature ;
			_OriginalSignature = NULL ;
			}
		if (_nOriginalFunctions > 0) {
			if (0 != ARE::ComputeSignature(_nOriginalFunctions, _OriginalFunctions, _OriginalWidth, _OriginalSignature)) 
				return 1 ;
			}
		else if (_nVars > 0) {
			_OriginalSignature = new int32_t[_nVars] ;
			if (NULL == _OriginalSignature) 
				return 1 ;
			for (int32_t i = 0 ; i < _nVars ; i++) 
				_OriginalSignature[i] = _Vars[i] ;
			_OriginalWidth = _nVars ;
			}
		if (_OriginalWidth < 0) 
			return 1 ;
		}

	int n_ = _OriginalWidth + _nAugmentedFunctions + _nIntermediateFunctions ;
	if (n_ <= 0) 
		{ _Width = 0 ; return 0 ; }

	int nF = 0 ;
	if (IncludeOriginalFunctions) nF += _nOriginalFunctions ;
	if (IncludeAugmentedFunctions) nF += _nAugmentedFunctions ;
	if (IncludeIntermediateFunctions) nF += _nIntermediateFunctions ;
	if (nF > 0) {
		std::vector<ARE::Function *> FNs ; 
		FNs.reserve(nF) ;
		if (FNs.capacity() != nF) 
			return 1 ;
		if (IncludeOriginalFunctions) {
			for (int32_t i = 0 ; i < _nOriginalFunctions ; i++) 
				FNs.push_back(_OriginalFunctions[i]) ;
			}
		if (IncludeAugmentedFunctions) {
			for (int32_t i = 0 ; i < _nAugmentedFunctions ; i++) 
				FNs.push_back(_AugmentedFunctions[i]) ;
			}
		if (IncludeIntermediateFunctions) {
			for (int32_t i = 0 ; i < _nIntermediateFunctions ; i++) 
				FNs.push_back(_IntermediateFunctions[i]) ;
			}

		if (NULL != _Signature) { delete [] _Signature ; _Signature = NULL ; }
		_Width = -1 ;

		if (0 != ARE::ComputeSignature(FNs.size(), FNs.data(), _Width, _Signature)) 
			return 1 ;
		}
	else if (_nVars > 0) {
		_Signature = new int32_t[_nVars] ;
		if (NULL == _Signature) 
			return 1 ;
		for (int32_t i = 0 ; i < _nVars ; i++) 
			_Signature[i] = _Vars[i] ;
		_Width = _nVars ;
		}

	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeOutputFunctionWithScopeWithoutTable(int32_t * & TempSpaceForArglist, int32_t TempSpaceForArglistSize, ARE::Function * & FN, int32_t & max_var)
{
	max_var = -1 ;
	if (NULL != FN) {
		FN->Destroy() ;
		FN->SetIDX(-(_V+1)) ;
		}
	if (NULL == _Workspace) 
		return 1 ;

	if (_Width < 0) {
		int32_t res = ComputeSignature() ;
		if (0 != res) {
//			if (createdFN && NULL != FN) 
//				{ delete FN ; FN = NULL ; }
			return 1 ;
			}
		}
	if (_Width <= 0) 
		return 0 ; // bucket has no variables; that should not happen
	// if _Width - _nVars <= 0, then all variables are eliminated; output is a const fn.

	bool createdFN = false ;
	if (NULL == FN) {
		FN = new ARE::Function(_Workspace, NULL != _Workspace ? _Workspace->Problem() : NULL, -(_V+1)) ;
		if (NULL == FN) 
			return 1 ;
		FN->SetOriginatingBucket(this) ;
		createdFN = true ;
		}

	// check memory is ok for sorting
	if (TempSpaceForArglistSize < _Width) {
		if (NULL != TempSpaceForArglist) 
			{ delete [] TempSpaceForArglist ; TempSpaceForArglist = NULL ; TempSpaceForArglistSize = 0 ; }
		TempSpaceForArglist = new int32_t[_Width] ;
		if (NULL == TempSpaceForArglist) 
			return 1 ;
		TempSpaceForArglistSize = _Width ;
		}

	// construct args list; find max argument.
	const int32_t *varpos = _Workspace->VarPos() ;
	int32_t nArgs = 0 ;
	for (int32_t i = 0 ; i < _Width ; i++) {
		int32_t v = _Signature[i], j ;
		for (j = 0 ; j < _nVars ; j++) 
			{ if (v == _Vars[j]) break ; }
		if (j < _nVars) 
			continue ;
		TempSpaceForArglist[nArgs++] = v ;
		if (max_var < 0) 
			max_var = v ;
		else if (varpos[max_var] < varpos[v]) 
			max_var = v ;
		}

	// set up fn; nArgs should be _Width - _nVars
	FN->SetArguments(nArgs, TempSpaceForArglist) ;

	return 0 ;
}


/*
int32_t BucketElimination::Bucket::ReorderFunctionScopesWrtParentBucket(bool IncludeOriginalFunctions, bool IncludeNewFunctions)
{
	// destroy table of all childbucketfunctions of this bucket;
	// reordering of scopes should be done in the beginning when these tables don't exist yet, 
	// so this should not be a problem.
	int32_t i ;
	if (IncludeNewFunctions) {
		for (i = 0 ; i < _nAugmentedFunctions ; i++) {
			ARE::Function *f = _AugmentedFunctions[i] ;
			f->DestroyTableData() ;
			}
		}

	// reorder scopes
	int32_t nOutputVars = _Width - _nVars ;
	if (nOutputVars < 1) {
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nAugmentedFunctions ; i++) {
				if (0 != _AugmentedFunctions[i]->ReOrderArguments(_nVars, _Vars, 0, NULL)) 
					return 1 ;
				}
			}
		if (IncludeOriginalFunctions) {
			for (i = 0 ; i < _nOriginalFunctions ; i++) {
				if (0 != _OriginalFunctions[i]->ReOrderArguments(_nVars, _Vars, 0, NULL)) 
					return 1 ;
				}
			}
		}
	else {
		for (MiniBucket & mb : _MiniBuckets) {
			}
		if (IncludeNewFunctions) {
			for (i = 0 ; i < _nAugmentedFunctions ; i++) {
				if (0 != _AugmentedFunctions[i]->ReOrderArguments(nOutputVars, _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		if (IncludeOriginalFunctions) {
			for (i = 0 ; i < _nOriginalFunctions ; i++) {
				if (0 != _OriginalFunctions[i]->ReOrderArguments(nOutputVars, _OutputFunction.Arguments(), _nVars, _Vars)) 
					return 1 ;
				}
			}
		}

	return 0 ;
}
*/

int32_t BucketElimination::Bucket::NoteOutputFunctionComputationCompletion(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		mb->NoteOutputFunctionComputationCompletion() ;
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeWMBEq(int32_t *ContextAssignment, double w_IS, double g_parent, double & fIntermediateValue, double & q)
{
	q = _Workspace->FnCombinationNeutralValue() ;

	// check if fInt is computed
	if (fIntermediateValue > 1.0e+32) {
		fIntermediateValue = _Workspace->FnCombinationNeutralValue() ;
		for (int32_t j = 0 ; j < nIntermediateFunctions() ; j++) {
			ARE::Function *f = IntermediateFunction(j) ;
			if (0 == f->N()) 
				_Workspace->ApplyFnCombinationOperator(fIntermediateValue, f->ConstValue()) ;
			else {
				__int64 adr = f->ComputeFnTableAdr(ContextAssignment, _Workspace->Problem()->K()) ;
				_Workspace->ApplyFnCombinationOperator(fIntermediateValue, f->TableEntry(adr)) ;
				}
			}
		}

	// compute h as the sum over all minibuckets : w_i * [\phi (x,Y) / m_i(x,Y)]^{1/w_i} 
	// where 
	//		w_i is WMBE weight of minibucket i
	//		/phi is the combination of all Original and Augmented functions in the minibucket
	//		m_i is the output function of the minibucket
	double h = _Workspace->VarEliminationDefaultValue() ;
	for (MiniBucket *mb : _MiniBuckets) {
		double WMBEweight = mb->WMBE_weight() ;
		double one_over_WMBEweight = WMBEweight < 1.0e+32 ? 1.0/WMBEweight : DBL_MAX ;
		double log_mbe_w = log10(WMBEweight) ;
		// compute combination of all MB functions
		double x = _Workspace->FnCombinationNeutralValue() ;
		for (int i = 0 ; i < mb->nFunctions() ; ++i) {
			ARE::Function *f = mb->Function(i) ;
			__int64 adr = f->ComputeFnTableAdr(ContextAssignment, _Workspace->Problem()->K()) ;
			_Workspace->ApplyFnCombinationOperator(x, f->TableEntry(adr)) ;
			}
		// divide the result with m_i; take the result to the power of 1/w; multiply with w_i
		ARE::Function & m = mb->OutputFunction() ; // m_i = outgoing msg of this minibucket
		double mValue ; // value of m_i corresponding to the given configuration
		if (m.N() > 0) {
			__int64 adr = m.ComputeFnTableAdr(ContextAssignment, _Workspace->Problem()->K()) ;
			mValue = m.TableEntry(adr) ;
			}
		else 
			mValue = m.ConstValue() ;
		if (_Workspace->Problem()->FunctionsAreConvertedToLogScale()) {
			x -= mValue ;
			x /= WMBEweight ;
			x += log_mbe_w ;
			}
		else {
			x /= mValue ;
			x = pow(x, one_over_WMBEweight) ;
			x *= WMBEweight ;
			}
		// add to h
		if (_Workspace->Problem()->FunctionsAreConvertedToLogScale()) {
			LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(h, h, x)
			}
		else 
			h += x ;
		}

	// q = w_IS * h * g_parent * fIntermediateValue
	q = w_IS ;
	_Workspace->ApplyFnCombinationOperator(q, h) ;
	_Workspace->ApplyFnCombinationOperator(q, g_parent) ;
	_Workspace->ApplyFnCombinationOperator(q, fIntermediateValue) ;

	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeCost(int32_t *FullAssignment, ARE_Function_TableType & Value)
{
	Value = _Workspace->FnCombinationNeutralValue() ;
	for (int32_t j = 0 ; j < nOriginalFunctions() ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
		_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeCostEx(int32_t *PathAssignment, ARE_Function_TableType & Value)
{
	Value = _Workspace->FnCombinationNeutralValue() ;
	for (int32_t j = 0 ; j < nOriginalFunctions() ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		__int64 adr = f->ComputeFnTableAdr_wrtLocalPermutation(PathAssignment, _Workspace->Problem()->K()) ;
		_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeHeuristic(int32_t *FullAssignment, ARE_Function_TableType & Value)
{
	Value = _Workspace->FnCombinationNeutralValue() ;
	for (int32_t j = 0 ; j < nAugmentedFunctions() ; j++) {
		ARE::Function *f = AugmentedFunction(j) ;
		if (f->N() <= 0) {
			_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
			}
		else {
			__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
			ARE_Function_TableType v = f->TableEntry(adr) ;
			f->CurrentValue() = v ;
			_Workspace->ApplyFnCombinationOperator(Value, v) ;
			}
		}
	// assume that for intermediate functions, CurrentValue was already computed (in an earlier bucket).
	for (int32_t j = 0 ; j < nIntermediateFunctions() ; j++) {
		ARE::Function *f = IntermediateFunction(j) ;
#ifdef RECOMPUTE_INT_FNS_in_H_COMPUTATION
		if (f->N() <= 0) {
			_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
			}
		else {
			__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
			_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
			}
#else
		_Workspace->ApplyFnCombinationOperator(Value, f->CurrentValue()) ;
#endif // RECOMPUTE_INT_FNS_in_H_COMPUTATION
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeHeuristicEx(int32_t *PathAssignment, ARE_Function_TableType & Value)
{
	Value = _Workspace->FnCombinationNeutralValue() ;
	for (int32_t j = 0 ; j < nAugmentedFunctions() ; j++) {
		ARE::Function *f = AugmentedFunction(j) ;
		if (f->N() <= 0) {
			_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
			}
		else {
			__int64 adr = f->ComputeFnTableAdr_wrtLocalPermutation(PathAssignment, _Workspace->Problem()->K()) ;
			ARE_Function_TableType v = f->TableEntry(adr) ;
			f->CurrentValue() = v ;
			_Workspace->ApplyFnCombinationOperator(Value, v) ;
			}
		}
	// assume that for intermediate functions, CurrentValue was already computed (in an earlier bucket).
	for (int32_t j = 0 ; j < nIntermediateFunctions() ; j++) {
		ARE::Function *f = IntermediateFunction(j) ;
#ifdef RECOMPUTE_INT_FNS_in_H_COMPUTATION
		if (f->N() <= 0) {
			_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
			}
		else {
			__int64 adr = f->ComputeFnTableAdr_wrtLocalPermutation(PathAssignment, _Workspace->Problem()->K()) ;
			_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
			}
#else
		_Workspace->ApplyFnCombinationOperator(Value, f->CurrentValue()) ;
#endif // RECOMPUTE_INT_FNS_in_H_COMPUTATION
		}
	return 0 ;
}


int32_t BucketElimination::Bucket::ComputeValue(int32_t *FullAssignment, int32_t Hints, ARE_Function_TableType & Value)
{
	Value = _Workspace->FnCombinationNeutralValue() ;
	if (0 != (Hints&1)) {
		for (int32_t j = 0 ; j < nOriginalFunctions() ; j++) {
			ARE::Function *f = OriginalFunction(j) ;
			__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
			_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
			}
		}

	if (0 != (Hints&2)) {
		for (int32_t j = 0 ; j < nAugmentedFunctions() ; j++) {
			ARE::Function *f = AugmentedFunction(j) ;
			if (f->N() <= 0) {
				_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
				}
			else {
				__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
				ARE_Function_TableType v = f->TableEntry(adr) ;
				f->CurrentValue() = v ;
				_Workspace->ApplyFnCombinationOperator(Value, v) ;
				}
			}
		}

	if (0 != (Hints&4)) {
		for (int32_t j = 0 ; j < nIntermediateFunctions() ; j++) {
			ARE::Function *f = IntermediateFunction(j) ;
#ifdef RECOMPUTE_INT_FNS_in_H_COMPUTATION
			if (f->N() <= 0) {
				_Workspace->ApplyFnCombinationOperator(Value, f->ConstValue()) ;
				}
			else {
				__int64 adr = f->ComputeFnTableAdr(FullAssignment, _Workspace->Problem()->K()) ;
				_Workspace->ApplyFnCombinationOperator(Value, f->TableEntry(adr)) ;
				}
#else
			_Workspace->ApplyFnCombinationOperator(Value, f->CurrentValue()) ;
#endif // RECOMPUTE_INT_FNS_in_H_COMPUTATION
			}
		}

	return 0 ;
}


__int64 BucketElimination::Bucket::ComputeProcessingComplexity(void)
{
	__int64 N = 0 ;
	for (MiniBucket *mb : _MiniBuckets) {
		__int64 n = mb->ComputeProcessingComplexity() ;
		if (n > 0) 
			N += n ;
		}
	return N ;
}


int32_t BucketElimination::Bucket::ComputeOutputFunctions(bool DoMomentMatching, signed char approx_bound, int64_t & TotalSumOutputFunctionsNumEntries)
{
	TotalSumOutputFunctionsNumEntries = 0 ;

	ARE::ARP *problem = _Workspace->Problem() ;
	int32_t varElimOperator = _Workspace->VarEliminationType() ;
	bool problem_is_optimization = problem->IsOptimizationProblem() ;
	bool problem_is_summation = problem->IsSummationProblem() ;

	// this is the avg max-marginals fn
	ARE::Function fAvgMM(_Workspace, _Workspace->Problem(), _MiniBuckets.size()) ;
	double *average_mm_table = NULL ;
	// this is max-marginals for each mb
	ARE::Function **max_marginals = NULL ;
	// this is joint-scope of all mb output functions
	int32_t *joint_scope = NULL, joint_scope_size = 0 ;

	// make sure width is computed
	int32_t min_width = INT_MAX, max_width = -INT_MAX ;
	for (MiniBucket *mb : _MiniBuckets) {
		mb->WMBE_weight() = 1.0 ;
		if (mb->Width() < 0) {
			if (0 != mb->ComputeSignature()) 
				return 1 ;
			if (mb->Width() < 0) 
				return 1 ;
			}
		if (min_width > mb->Width()) 
			min_width = mb->Width() ;
		if (max_width < mb->Width()) 
			max_width = mb->Width() ;
		}
	if (min_width <= 0 || max_width <= 0) {
		min_width = max_width = 0 ; // most likely nMiniBuckets is 0
		}

	int32_t res = 1 ;
	int32_t idx = 0 ;
	bool mbe_computed = false ;

	if (_MiniBuckets.size() > 1) {
		if (problem->IsMaximizationProblem() && approx_bound < 0) 
			// as of now (2017-03) cannot do lower bound of max problems
			goto done_MM ;
		else if (problem->IsMinimizationProblem() && approx_bound > 0) 
			// as of now (2017-03) cannot do upper bound of min problems
			goto done_MM ;
		}

	if (_MiniBuckets.size() > 1 && DoMomentMatching) {
		if (problem_is_optimization) {
			// compute joint scope of all MBs; next line, allocate some more space for other arrays
			joint_scope = new int32_t[min_width + max_width + max_width] ; if (NULL == joint_scope) goto done_MM ;
			int32_t *temp_scope = joint_scope + min_width, temp_scope_size ; 
			int32_t *temp2_scope = temp_scope + max_width ;
			bool first_mb = true ;
			for (MiniBucket *mb : _MiniBuckets) {
				const int32_t *sig = mb->SortedSignature() ;
				if (first_mb) 
					{ for (int32_t i = 0 ; i < mb->Width() ; i++) joint_scope[i] = sig[i] ; joint_scope_size = mb->Width() ; first_mb = false ; }
				else 
					{ ARE::SetIntersection(joint_scope, joint_scope_size, sig, mb->Width()) ; }
				if (0 == joint_scope_size) 
					goto done_MM ; // error; this is not supposed to happen
				}

			// compute max-marginals for each MB.
			max_marginals = new ARE::Function*[_MiniBuckets.size()] ; if (NULL == max_marginals) goto done_MM ;
			memset(max_marginals, 0, sizeof(ARE::Function*)*_MiniBuckets.size()) ;
			idx = 0 ;
			for (MiniBucket *mb : _MiniBuckets) {
				const int32_t *sig = mb->SortedSignature() ;
				for (int32_t i = 0 ; i < mb->Width() ; i++) temp_scope[i] = sig[i] ; temp_scope_size = mb->Width() ;
				ARE::SetMinus(temp_scope, temp_scope_size, joint_scope, joint_scope_size) ;
				max_marginals[idx] = new ARE::Function(_Workspace, _Workspace->Problem(), idx) ;
				if (NULL == max_marginals[idx]) 
					goto done_MM ;
				if (0 != mb->ComputeOutputFunction(varElimOperator, *(max_marginals[idx]), temp_scope, temp_scope_size, temp2_scope, DBL_MAX)) 
					goto done_MM ;
				++idx ;
				}
			// NOTE : all max_marginals[*] have same scope and the order of arguments is the same (i.e. argument lists are identical).

			// compute avg max-marginals
			INT64 table_size = max_marginals[0]->TableSize() ;
			if (table_size <= 0) 
				goto done_MM ;
			average_mm_table = new double[table_size] ;
			if (NULL == average_mm_table) 
				goto done_MM ;
			average_mm_table[0] = _Workspace->FnCombinationNeutralValue() ;
			for (int32_t i = 1 ; i < table_size ; i++) average_mm_table[i] = average_mm_table[0] ;
			for (int32_t j = _MiniBuckets.size() - 1 ; j >= 0 ; j--) {
				for (int32_t i = 0; i < table_size ; i++) 
					_Workspace->ApplyFnCombinationOperator(average_mm_table[i], max_marginals[j]->TableEntry(i)) ;
				}
			double N_ = (double) _MiniBuckets.size() ;
			for (int32_t i = 0; i < table_size ; i++) 
				average_mm_table[i] /= N_ ;

			if (0 != fAvgMM.SetArguments(max_marginals[0]->N(), max_marginals[0]->Arguments())) 
				goto done_MM ;
			if (0 != fAvgMM.SetTableData(table_size, average_mm_table)) 
				goto done_MM ;
			average_mm_table = NULL ;

			idx = 0 ;
			for (MiniBucket *mb : _MiniBuckets) {
				if (0 != mb->ComputeOutputFunction(varElimOperator, max_marginals[idx], &fAvgMM, DBL_MAX)) 
					goto done_MM ;
				++idx ;
				}
			mbe_computed = true ;
			}
		else if (problem_is_summation && approx_bound > 0) {
			// set w_i to 1/nMB; do weighted MBE, using the fU functions computed earlier
			double nb = _MiniBuckets.size() ;
			if (nb <= 0.0) goto do_default ;
			double w = 1.0/nb ;
#ifdef DO_PROD_SUM_MM
			// do local cost shifting; from paper "AND/OR Search for Marginal MAP" R.Marinescu, R.Dechter, A.Ihler.
			// in each minibucket, compute vector u_MB(X_B), where X_B is the variable of the bucket.
			std::vector<ARE::Function> fUs ; fUs.reserve(_MiniBuckets.size()) ; fUs.resize(_MiniBuckets.size()) ;
			for (int32_t i = 0 ; i < fUs.size() ; ++i) 
				{ ARE::Function & fU = fUs[i] ; fU.Initialize(_Workspace, _Workspace->Problem(), i) ; fU.SetArguments(1, &_V) ; }
			int32_t vars[MAX_NUM_VARIABLES_PER_BUCKET], temp_space_for_vars[MAX_NUM_VARIABLES_PER_BUCKET] ; int32_t nMBvars ;
			for (int32_t i = 0 ; i < fUs.size() ; ++i) {
				ARE::Function & fU = fUs[i] ;
				BucketElimination::MiniBucket *mb = _MiniBuckets[i] ;
				if (mb->Width() < 0) 
					{ if (0 != mb->ComputeSignature()) goto do_default ; }
				const int32_t *sig = mb->SortedSignature() ;
				// collect all vars of the MB, except X_B
				nMBvars = 0 ;
				for (int32_t j = 0 ; j < mb->Width() ; ++j) 
					{ if (sig[j] != _V) vars[nMBvars++] = sig[j] ; }
				int32_t res_mbe = mb->ComputeOutputFunction(VAR_ELIMINATION_TYPE_SUM, fU, vars, nMBvars, temp_space_for_vars, w) ;
				if (0 != res_mbe) 
					goto done_MM ;
				}
			// compute the combined u(X_B)
			ARE::Function FU ; FU.Initialize(_Workspace, _Workspace->Problem(), 0) ; FU.SetArguments(1, &_V) ; FU.AllocateTableData() ;
			int32_t k = _Workspace->Problem()->K(_V) ;
			double *FUtable = FU.TableData() ;
			for (int32_t j = 0 ; j < k ; ++j) {
				FUtable[j] = _Workspace->FnCombinationNeutralValue() ;
				for (int32_t i = 0 ; i < fUs.size() ; ++i) {
					ARE::Function & fU = fUs[i] ;
					double v = fU.TableEntry(j) ;
					if (_Workspace->ProblemIsLogScale()) v *= w ; else v = pow(v, w) ;
					_Workspace->ApplyFnCombinationOperator(FUtable[j], v) ;
					}
				}
			// do weighted MBE, using the fU functions computed earlier
			for (int32_t i = 0 ; i < _MiniBuckets.size() ; ++i) {
				BucketElimination::MiniBucket *mb = _MiniBuckets[i] ;
				ARE::Function & fU = fUs[i] ;
				mb->WMBE_weight() = w ;
				if (0 != mb->ComputeOutputFunction(varElimOperator, &FU, &fU, w)) 
					goto done_MM ;
				}
#else // DO_PROD_SUM_MM
			// do weighted MBE
			for (int32_t i = 0 ; i < _MiniBuckets.size() ; ++i) {
				BucketElimination::MiniBucket *mb = _MiniBuckets[i] ;
				mb->WMBE_weight() = w ;
				if (0 != mb->ComputeOutputFunction(varElimOperator, NULL, NULL, w)) 
					goto done_MM ;
				}
#endif // DO_PROD_SUM_MM
			mbe_computed = true ;
			}
		else if (problem_is_summation && approx_bound < 0) {
			// set w_0 to 2, others to 1/(nMB-1); do weighted MBE
			double nb = _MiniBuckets.size() ;
			if (nb <= 0.0) goto do_default ;
			double w_0 = 2.0 ;
			double w = -1.0/(nb-1) ;
			idx = 0 ;
			for (MiniBucket *mb : _MiniBuckets) {
				mb->WMBE_weight() = 0 == idx ? w_0 : w ;
				if (0 != mb->ComputeOutputFunction(varElimOperator, NULL, NULL, 0 == idx ? w_0 : w)) 
					goto done_MM ;
				++idx ;
				}
			mbe_computed = true ;
			}
		}

do_default :
	if (! mbe_computed) {
		idx = 0 ;
		for (MiniBucket *mb : _MiniBuckets) {
			if (0 != mb->ComputeOutputFunction(varElimOperator, NULL, NULL, DBL_MAX)) 
				goto done_MM ;
			++idx ;
			// if problem is summation, sum over first mini-bucket, max (or min) over other minibuckets.
			if (VAR_ELIMINATION_TYPE_SUM == varElimOperator) 
				varElimOperator = approx_bound >= 0 ? VAR_ELIMINATION_TYPE_MAX : VAR_ELIMINATION_TYPE_MIN ;
			}
		mbe_computed = true ;
		}

	res = 0 ;

done_MM :

	// if there are const functions, set their current_value equal to const_value
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::Function & f = mb->OutputFunction() ;
		if (f.N() <= 0) 
			f.CurrentValue() = f.ConstValue() ;
		else {
			if (mbe_computed) {
				TotalSumOutputFunctionsNumEntries += f.TableSize() ;
				}
			}
		}

	// done; delete stuff.
	if (NULL != average_mm_table) delete [] average_mm_table ;
	if (NULL != max_marginals) {
		for (int32_t j = 0 ; j < _MiniBuckets.size() ; j++) { if (NULL != max_marginals[j]) delete max_marginals[j] ; }
		delete [] max_marginals ;
		}
	if (NULL != joint_scope) 
		delete [] joint_scope ;

	if (NULL != _PostComputeExtNote) {
		for (MiniBucket *mb : _MiniBuckets) {
			(*_PostComputeExtNote)(*this,*mb) ;
			}
		}

	return res ;
}


int32_t MBpartitionFnSortCompFn(void *Obj1, void *Obj2)
{
	ARE::Function *f1 = (ARE::Function*) Obj1 ;
	ARE::Function *f2 = (ARE::Function*) Obj2 ;
	if (f1->N() > f2->N()) 
		return -1 ;
	else if (f1->N() < f2->N()) 
		return  1 ;
	int32_t *args1 = f1->SortedArgumentsList(true) ;
	int32_t *args2 = f2->SortedArgumentsList(true) ;
	for (int32_t i = 0 ; i < f1->N() ; ++i) {
		if (args1[i] < args2[i]) return -1 ;
		if (args1[i] > args2[i]) return  1 ;
		}
	return 0 ;
}

int32_t BucketElimination::Bucket::CreateMBPartitioning(int32_t iBound, bool CreateTables, bool doMomentMatching, signed char Bound, bool & AbandonIfActualPartitioning, std::vector<int32_t> & key, std::vector<int64_t> & data, std::vector<int32_t> & helperArray)
{
	bool abandonIfActualPartitioning = AbandonIfActualPartitioning ;
	AbandonIfActualPartitioning = false ;

//INT64 tB = ARE::GetTimeInMilliseconds() ;

	// check that all children have already been partitioned
	for (int32_t i = 0 ; i < _nChildren ; i++) {
		int32_t child = _ChildVars[i] ;
		Bucket *b = _Workspace->MapVar2Bucket(child) ;
		if (NULL == b) continue ;
		if (b->nMiniBuckets() <= 0) 
			return 1 ;
		}

	// sort functions in the order of decreasing scope size
	int32_t nF = nOriginalFunctions() + nAugmentedFunctions() ;
	if (key.capacity() < nF) { key.reserve(nF) ; if (key.capacity() < nF) goto failed ; }
	if (data.capacity() < nF) { data.reserve(nF) ; if (data.capacity() < nF) goto failed ; }
	key.clear() ; data.clear() ;
	for (int32_t j = _nOriginalFunctions - 1 ; j >= 0 ; j--) {
		ARE::Function *f = _OriginalFunctions[j] ;
		key.push_back(-f->N()) ;
		data.push_back((INT64) f) ;
		if (data.size() != key.size()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : data_size!=key_size(1); i=%d, v=%d, nMBs=%d", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size()) ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
		}
	for (int32_t j = _nAugmentedFunctions - 1 ; j >= 0 ; j--) {
		ARE::Function *f = _AugmentedFunctions[j] ;
		key.push_back(-f->N()) ;
		data.push_back((INT64) f) ;
		if (data.size() != key.size()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : data_size!=key_size(2); i=%d, v=%d, nMBs=%d", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size()) ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
		}
	int32_t left[32], right[32] ;
	QuickSortLong_i64(key.data(), key.size(), data.data(), left, right) ;

	// do secondary sort so that functions with equal scope size are also consistenty/persistently sorted
	int32_t j ; j = 0 ;
	for (int32_t i = 1 ; i <= data.size() ; ++i) {
		ARE::Function *f = i < data.size() ? (ARE::Function*) data[i] : NULL ;
		if (NULL != f ? f->N() != ((ARE::Function*) data[j])->N() : true) {
			int32_t n = i - j, k ;
			if (n > 1) {
				// sort data [j,i); use sorted args list
				for (k = j ; k < i ; ++k) {
					ARE::Function *f_ = (ARE::Function*) data[k] ;
					int32_t *args = f_->SortedArgumentsList(true) ;
					if (NULL == args) break ;
					}
				if (k < i) continue ;
				QuickSort((void**) (data.data()+j), n, left, right, MBpartitionFnSortCompFn) ;
				}
			if (NULL == f) 
				break ;
			j = i ;
			}
		}

/*
INT64 tS = ARE::GetTimeInMilliseconds() ;
if (_IDX < 2000) {
INT64 dt = tS - tB ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, iBound=%d, nFNs=%d, sorting time=%lld[msec]", (int32_t) _IDX, (int32_t) _V, (int32_t) iBound, (int32_t) data.size(), dt) ;
::fflush(ARE::fpLOG) ; }*/

#ifdef _DEBUG
		if (_MiniBuckets.size() > 0) {
			int32_t error = 1 ;
			}
		_MiniBuckets.clear() ;
#endif

//INT64 dtSUM1 = 0, dtSUM2 = 0, dtSUM3 = 0 ; ;

		int32_t MBmaxsize ; MBmaxsize = iBound+1 ;
		// process functions, one at a time; for each fn, try to add to an existing MB.
		for (int32_t j = 0 ; j < key.size() ; j++) {
			ARE::Function *f = (ARE::Function*) data[j] ;
			bool placed = false ;
			if (iBound > 0) {
				for (MiniBucket *mb : _MiniBuckets) {
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
					int32_t res = mb->AllowsFunction(*f, MBmaxsize, helperArray) ;
//INT64 t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
					if (res < 0) 
						goto failed ;
					if (res > 0) {
						placed = true ;
//INT64 t1_ = ARE::GetTimeInMilliseconds() ;
						if (0 != mb->AddFunction(*f, helperArray)) 
							goto failed ;
//INT64 t2_ = ARE::GetTimeInMilliseconds() ;
//dtSUM2 += t2_ - t1_ ;
						break ;
						}
					}
				}
			if (! placed) {
//int ss = _MiniBuckets.size() ;
				if (_MiniBuckets.size() > 0 && abandonIfActualPartitioning) {
					AbandonIfActualPartitioning = true ;
					goto failed ;
					}
//INT64 t3 = ARE::GetTimeInMilliseconds() ;
				MiniBucket *mb = new MiniBucket ;
				if (NULL == mb) 
					goto failed ;
				int32_t existing_size = _MiniBuckets.size() ;
				_MiniBuckets.push_back(mb) ;
				if (_MiniBuckets.size() != (1+existing_size)) 
					{ delete mb ; goto failed ; }
				mb->Initalize(*this, existing_size) ;
				for (int32_t idxV = 0 ; idxV < _nVars ; idxV++) 
					mb->AddVar(_Vars[idxV]) ;
				if (0 != mb->AddFunction(*f, helperArray)) 
					goto failed ;
//INT64 t4 = ARE::GetTimeInMilliseconds() ;
//dtSUM3 += t4 - t3 ;
				}
			}

/*INT64 tMB = ARE::GetTimeInMilliseconds() ;
int32_t div = _IDX%1000 ;
if (_IDX < 2000) {
	INT64 dt = tMB - tS ;
	fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, nMBs=%d, nFNs=%d, MB partitioning time=%lld(check=%lld, addfn=%lld, new=%lld)", (int32_t) _IDX, (int32_t) _V, (int32_t) _MiniBuckets.size(), (int32_t) data.size(), dt, dtSUM1, dtSUM2, dtSUM3) ;
	::fflush(ARE::fpLOG) ;
	}*/

	// create output fn for each mb...
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::FnConstructor f = MBoutputfncnstrctr(*mb) ;
		mb->CreateOutputFunction(*f) ;
		}

	// minibuckets for current bucket are now ready, process each and place resulting function
//dtSUM1 = dtSUM2 = dtSUM3 = 0 ; ;
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::Function & f = mb->OutputFunction() ;
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
		if (0 != mb->ComputeOutputFunctionWithScopeWithoutTable()) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb->ComputeOutputFunctionWithScopeWithoutTable() failed ...") ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
//INT64 t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
		Bucket *target_bucket = f.Bucket() ; // note target_bucket may be NULL (i.e. when 0==f.N())
		for (Bucket *mb_b = ParentBucket() ; mb_b != target_bucket && NULL != mb_b ; mb_b = mb_b->ParentBucket()) {
//t1 = ARE::GetTimeInMilliseconds() ;
//dtSUM1 += t2 - t1 ;
			if (0 != mb_b->AddIntermediateFunction(f)) {
				fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb_b->AddIntermediateFunction() failed ...") ;
				::fflush(ARE::fpLOG) ;
				goto failed ;
				}
//t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM2 += t2 - t1 ;
			}
		if (NULL != target_bucket) {
//t1 = ARE::GetTimeInMilliseconds() ;
			if (0 != target_bucket->AddAugmentedFunction(f)) {
				fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : target_bucket->AddAugmentedFunction() failed ...") ;
				::fflush(ARE::fpLOG) ;
				goto failed ;
				}
//t2 = ARE::GetTimeInMilliseconds() ;
//dtSUM3 += t2 - t1 ;
			}
		}

/*INT64 tE = ARE::GetTimeInMilliseconds() ;
if (_IDX < 2000) {
INT64 dt = tE - tMB ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, MB placing  time=%lld[msec] compOFN=%lld +IFN=%lld +AFN=%lld", (int32_t) _IDX, (int32_t) _V, dt, dtSUM1, dtSUM2, dtSUM3) ;
::fflush(ARE::fpLOG) ; }*/

	// if requested, do moment-matching
	if (CreateTables) {
//INT64 t1 = ARE::GetTimeInMilliseconds() ;
		int64_t totalSumOutputFunctionsNumEntries = 0 ;
		if (0 != ComputeOutputFunctions(doMomentMatching, Bound, totalSumOutputFunctionsNumEntries)) {
			fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : mb->ComputeOutputFunction() failed ...") ;
			::fflush(ARE::fpLOG) ;
			goto failed ;
			}
/*INT64 t2 = ARE::GetTimeInMilliseconds() ;
INT64 dt = t2 - t1 ;
if (_IDX < 2000) {
INT64 dt = t2 - t1 ;
fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning; i=%d, v=%d, ComputeOutputFunctions doMomentMatching=%c time=%lld[msec]", (int32_t) _IDX, (int32_t) _V, doMomentMatching ? 'Y' : 'N', dt) ;
::fflush(ARE::fpLOG) ; }*/
		}

	return 0 ;
failed :
	fprintf(ARE::fpLOG, "\n   BucketElimination::Bucket::CreateMBPartitioning ERROR : failed label reached (AbandonIfActualPartitioning=%c) ...", AbandonIfActualPartitioning ? 'Y' : 'N') ;
	::fflush(ARE::fpLOG) ;
	DestroyPartitioning() ;
	return 1 ;
}


int32_t BucketElimination::Bucket::ComputeFirstVariableDistribution(ARE_Function_TableType *dist)
{
	int32_t i, j, k, ret = 0 ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;
	int32_t w = Width() ;
	if (w < 0) {
		ComputeSignature() ;
		w = Width() ;
		if (w < 0) 
			return 0 ; // should not happen; means bucket has no functions/variables.
		}
	const int32_t *signature = Signature() ;
	const int32_t nOF = nOriginalFunctions() ;
	const int32_t nAF = nAugmentedFunctions() ;
	const int32_t nIF = nIntermediateFunctions() ;
	int32_t nFunctions_OA = nOF + nAF ;
	int32_t nTotalFunctions = nOF + nAF + nIF ;
	if (w < 0 || nTotalFunctions < 1) {
		return 0 ;
		}

	if (w > MAX_NUM_VARIABLES_PER_BUCKET) 
		return ERRORCODE_too_many_variables ;
	if (nTotalFunctions > MAX_NUM_FUNCTIONS_PER_BUCKET) 
		return ERRORCODE_too_many_functions ;

	int32_t values[MAX_NUM_VARIABLES_PER_BUCKET] ; // this is the current value combination of the arguments of this function (table block).
	ARE::Function *flist[MAX_NUM_FUNCTIONS_PER_BUCKET] ; // this is a list of input functions

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	int32_t nFNs = 0 ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nFunctions_OA ; j++) {
		ARE::Function *f = AugmentedFunction(j - nOF) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}
	for (; j < nTotalFunctions ; j++) {
		ARE::Function *f = IntermediateFunction(j - nFunctions_OA) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist[nFNs++] = f ; f->ComputeArgumentsPermutationList(w, signature) ; }
		}

	INT64 ElimSize = 1 ;
	for (j = 1 ; j < w ; j++) 
		ElimSize *= problem->K(signature[j]) ;

	for (j = 0 ; j < w ; j++) 
		values[j] = 0 ;
	int32_t v0 = _Vars[0] ;
	for (i = 0 ; i < problem->K(v0) ; i++) {
		ARE_Function_TableType V = bews->VarEliminationDefaultValue() ;
		for (INT64 ElimIDX = 0 ; ElimIDX < ElimSize ; ElimIDX++) {
			ARE_Function_TableType value = nv ;
			for (j = 0 ; j < nFNs ; j++) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
				INT64 adr = flist[j]->ComputeFnTableAdr_wrtLocalPermutation(values, problem->K()) ;
				bews->ApplyFnCombinationOperator(value, flist[j]->TableEntry(adr)) ;
				}
			bews->ApplyVarEliminationOperator(V, value) ;
			// go to next argument value combination
			ARE::EnumerateNextArgumentsValueCombination(w-1, signature+1, values+1, problem->K()) ;
			}
		bews->ApplyFnCombinationOperator(V, const_factor) ;
		dist[i] = V ;
		values[0]++ ;
		}
done :
	return ret ;
}


int32_t BucketElimination::Bucket::ComputeFirstVariableDistributionEx(int32_t *ContextValues, ARE_Function_TableType *dist)
{
	int32_t j ;

	if (1 != _nVars) 
		return ERRORCODE_generic ;

	MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
	if (NULL == bews) 
		return ERRORCODE_generic ;
	ARE::ARP *problem = bews->Problem() ;
	if (NULL == problem) 
		return ERRORCODE_generic ;
	ARE_Function_TableType nv = bews->FnCombinationNeutralValue() ;

	const int32_t nOF = nOriginalFunctions() ;
	const int32_t nAF = nAugmentedFunctions() ;
	const int32_t nIF = nIntermediateFunctions() ;
	int32_t nTotalFunctions = nOF + nAF + nIF ;
	if (nTotalFunctions < 1) {
		return 0 ;
		}

	std::vector<ARE::Function *> flist ; // this is a list of input functions
	flist.reserve(nTotalFunctions) ;
	if (flist.capacity() != nTotalFunctions) 
		return ERRORCODE_out_of_memory ;

	ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
	for (j = 0 ; j < nOF ; j++) {
		ARE::Function *f = OriginalFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist.push_back(f) ; }
		}
	for (j = 0 ; j < nAF ; j++) {
		ARE::Function *f = AugmentedFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else { flist.push_back(f) ; }
		}
	for (j = 0 ; j < nIF ; j++) {
		ARE::Function *f = IntermediateFunction(j) ;
		if (NULL == f) continue ;
		if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
		else {
			// all variables in the scope should be assigned
			INT64 adr = f->ComputeFnTableAdr(ContextValues, problem->K()) ;
			bews->ApplyFnCombinationOperator(const_factor, f->TableEntry(adr)) ;
			}
		}

	int32_t val_bucket_var = ContextValues[_V] ;
	for (int32_t i = 0 ; i < problem->K(_V) ; i++) {
		ContextValues[_V] = i ;
		ARE_Function_TableType V = const_factor ;
		for (j = flist.size() - 1 ; j >= 0 ; j--) { // note : it should be that 0 != flist[j]->N(). note : it is assumed that flist[j] has 1 block. note : it is assumed that order of flist[j] arguments is the same as signature.
			INT64 adr = flist[j]->ComputeFnTableAdr(ContextValues, problem->K()) ;
			bews->ApplyFnCombinationOperator(V, flist[j]->TableEntry(adr)) ;
			}
		dist[i] = V ;
		}
	ContextValues[_V] = val_bucket_var ;

	return 0 ;
}


int32_t BucketElimination::Bucket::DestroyPartitioning(void)
{
	for (MiniBucket *mb : _MiniBuckets) {
		ARE::Function & f = mb->OutputFunction() ;
		BucketElimination::Bucket *B = f.Bucket() ;
		if (NULL == B) 
			continue ; // only way this is possible is that this fn is const fn
		// remove output function from their buckets.
		for (BucketElimination::Bucket *b = ParentBucket() ; NULL != b && B != b ; b = b->ParentBucket()) {
			b->RemoveIntermediateFunction(f) ;
			}
		B->RemoveAugmentedFunction(f, true) ;
		// destroy output functio; so that its table can be released.
		f.Destroy() ;
		// destroy mb
		mb->Destroy() ;
		delete mb ;
		}
	_MiniBuckets.clear() ;

	return 0 ;
}

