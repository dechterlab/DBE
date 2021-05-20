// AbsSampling.cpp : Defines the entry point for the console application.

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <thread>
#include <time.h>
#include <cmath>
#include <limits>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#else
#include <sys/time.h>
#endif // WINDOWS

#include "CVO/VariableOrderComputation.hxx"
#include "BE/Bucket.hxx"
#include "BE/MBEworkspace.hxx"
#include "AbsSamplingWorkspace.hxx"

// 3000001 = first version with version info
// 3000002 = proper/non-proper both in same AbsSampling.cpp
// 3000003 = randomized abstractions
#define AS_VERSION 3000006

//#define PERFORM_SINGLTON_CONSISTENCY_CHECK

/*
-fUAI c:\uci\problems\uai08_test2.uai -fVO c:\uci\problems\uai08_test2.vo -fEV c:\uci\problems\uai08_test2.evi
-treeType OR -nR 1 -iB 15 -a proper -nCustomProper 0 -fUAI C:\UCI\problems\ObjectDetection.uai -fVO C:\UCI\problems\ObjectDetection.elim -fEV C:\UCI\problems\ObjectDetection.uai.evid -t 3600
-treeType OR -nR 1 -iB  4 -a customproper -nCustomProper 0 -fUAI C:\UCI\problems\ObjectDetection_11.uai -fVO C:\UCI\problems\ObjectDetection_11.elim -fEV C:\UCI\problems\ObjectDetection_11.uai.evid -t 3600
-treeType AO -nR 700000000 -iB 15 -a customproper -nCustomProper 1 -fUAI C:\UCI\problems\2bitcomp_5.cnf.uai -fVO C:\UCI\problems\2bitcomp_5.cnf.elim -fEV C:\UCI\problems\2bitcomp_5.cnf.uai.evid -t 3600 -printPeriod 100

simple test problem; AND/OR tree; unique abstraction
-treeType AO -nR 1 -iB 999 -a unique -fUAI C:\UCI\problems\handbuilt_burglary.uai -fVO C:\UCI\problems\handbuilt_burglary.vo -fEV C:\UCI\problems\0.evid -t 3600 -printPeriod 100

2017-06-07 DFS produces seg fault
-nR 1 -iB 10 -a customproper -nCustomProper 1 -fUAI C:\UCI\problems\Grids_13.uai -fVO C:\UCI\problems\Grids_13.elim  -fEV C:\UCI\problems\Grids_13.uai.evid -t 60 -treeType "OR" -nLevelsLimit -1

2007-06-20 test number of DFS generated nodes
-nR 1 -iB 10 -a proper -fUAI C:\UCI\problems\Grids_13.uai -fVO C:\UCI\problems\Grids_13.elim  -fEV C:\UCI\problems\Grids_13.uai.evid -t 60 -treeType "AO" -nLevelsLimit -1

-nR 1 -igEH 0 -iB 12 -a customproper -nCustomProper 2 -fUAI Grids_15.uai -fVO Grids_15.elim -fEV Grids_15.uai.evid -t 60 -treeType "AO" -nLevelsLimit -1
-treeType OR -nR 1000 -iB 10 -a customproper -nCustomProper 3 -fUAI BN_111.uai -fVO BN_111.elim -fEV BN_111.uai.evid -t 3600
*/

#if defined WINDOWS || _WINDOWS
#else
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif // WINDOWS

double logabsdiffexp10(double num1, double num2){
    if(num1 == num2){
        return -std::numeric_limits<double>::infinity();
    }
    num1 = num1*std::log(10.0);
    num2 = num2*std::log(10.0);
    double num_max = num1 < num2 ? num2 : num1 ; // std::max(num1, num2);
    double num_min = num1 > num2 ? num2 : num1 ; // std::min(num1, num2);
    double ans = (num_min + std::log(expm1(num_max - num_min)))/std::log(10.0) ;
    // std::cout << ans << std::endl;
    return ans;
}

int main(int argc, char* argv[])
{
#if defined WINDOWS || _WINDOWS
	time_t wall0 ;
	time(&wall0) ;
#else
    double wall0 = get_wall_time();
	double cpu0  = get_cpu_time();
    double wall1, cpu1;
#endif // WINDOWS

	int nParams = argc ;
	unsigned long randomGeneratorSeed = 0 ;
	std::string problem_filename, evidence_filename, vo_filename ;
	std::string tree_type = "OR";
	int32_t nrunstodo = 1 ;
	int32_t nArgs = (nParams-1)>>1 ;
	int32_t nContext = 0, nAbs = -1 ;
	int32_t time_limit = 60;
	int32_t print_period = 1;
	int64_t nLevelsLimit = INT64_MAX ;
	bool findPracticalVariableOrder = true ;
	bool properAlg = true ;
	bool verbose = false ;
	bool abs_indv = false ;
	bool exclude_prep_time = false;
	bool forORtreetype_use_DFS_order = false;
	int iB = -1 ;
	ARE::VarElimOrderComp::ObjectiveToMinimize objCodePrimary = ARE::VarElimOrderComp::Width, objCodeSecondary = ARE::VarElimOrderComp::None ;
	std::string algorithm ;
	double exactZ = 0.0;

	if (1 + 2*nArgs != nParams) {
		printf("\nBAD COMMAND LINE; will exit ...") ;
		return 1 ;
		}
	
	for (int i = 0 ; i < nArgs; i++) {
		std::string sArgID = (NULL != argv[1+2*i] ? argv[1+2*i] : "") ;
		std::string sArg   = (NULL != argv[2*(i+1)] ? argv[2*(i+1)] : "") ;
		if (0 == stricmp("-s", sArgID.c_str())) 
			randomGeneratorSeed = std::strtoul(sArg.c_str(), NULL, 0) ;
		else if (0 == stricmp("-fUAI", sArgID.c_str()))
			problem_filename = sArg ;
		else if (0 == stricmp("-fEV", sArgID.c_str()))
			evidence_filename = sArg ;
		else if (0 == stricmp("-fVO", sArgID.c_str()))
			vo_filename = sArg ;
		else if (0 == stricmp("-treeType", sArgID.c_str()))
			tree_type = sArg ;
		else if (0 == stricmp("-printPeriod", sArgID.c_str()))
			print_period = atoi(sArg.c_str());
		else if (0 == stricmp("-nR", sArgID.c_str()))
			nrunstodo = atoi(sArg.c_str());
		else if (0 == stricmp("-t", sArgID.c_str()))
			time_limit = atoi(sArg.c_str());
		else if (0 == stricmp("-nContext", sArgID.c_str()))
			nContext = atoi(sArg.c_str());
		else if (0 == stricmp("-nAbs", sArgID.c_str()))
			nAbs = atoi(sArg.c_str());
		else if (0 == stricmp("-iB", sArgID.c_str()))
			iB = atoi(sArg.c_str());
		else if (0 == stricmp("-Verbose", sArgID.c_str()))
			verbose = atoi(sArg.c_str()) > 0 ;
		else if (0 == stricmp("-igEH", sArgID.c_str())) 
			AndOrSearchSpace::EXPLOIT_EXACT_H = 0 == atoi(sArg.c_str()) ;
		else if (0 == stricmp("-a", sArgID.c_str()))
			algorithm = sArg ;
		else if (0 == stricmp("-fpvo", sArgID.c_str()))
			findPracticalVariableOrder = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0] ;
		else if (0 == stricmp("-proper", sArgID.c_str()))
			properAlg = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0] ;
		else if (0 == stricmp("-absindv", sArgID.c_str()))
			abs_indv = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0] ;
		else if (0 == stricmp("-xPrepTime", sArgID.c_str()))
			exclude_prep_time = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0];
		else if (0 == stricmp("-xORdfs", sArgID.c_str()))
			forORtreetype_use_DFS_order = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0];
		else if (0 == stricmp("-O1", sArgID.c_str()))
			objCodePrimary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
		else if (0 == stricmp("-O2", sArgID.c_str()))
			objCodeSecondary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
		else if (0 == stricmp("-nLevelsLimit", sArgID.c_str()))
			nLevelsLimit = atoll(sArg.c_str());
		else if (0 == stricmp("-exactZ", sArgID.c_str()))
			exactZ = atof(sArg.c_str());
	}
	
	if (nrunstodo < 0) 
		nrunstodo = 0 ;

	int maxNumProcessorThreads = std::thread::hardware_concurrency() ; // this requires c++11
	if (maxNumProcessorThreads <= 0) 
		maxNumProcessorThreads = 1 ;

	if (0 == problem_filename.length() || 0 == vo_filename.length()) 
		return 1 ;

	ARE::ARP p("AbsSamplingProblem") ;
	p.SetOperators(FN_COBINATION_TYPE_PROD, VAR_ELIMINATION_TYPE_SUM) ;
	int resLoad = p.LoadFromFile(problem_filename) ;
	if (0 != resLoad) 
		return 10 ;
	int resPostCA = p.PerformPostConstructionAnalysis() ;
	if (0 != resPostCA) 
		return 11 ;
	if (iB > p.N()) 
		iB = p.N() ;

	int32_t nEV = -1 ;
	if (evidence_filename.length() > 0) {
		int resLoadE = p.LoadFromFile_Evidence(evidence_filename, nEV) ;
		if (0 != resLoadE) 
			return 20 ;
		int resElimE = p.EliminateEvidence() ;
		if (0 != resElimE) 
			return 21 ;
	}

	int resELDV = p.EliminateSingletonDomainVariables() ;
	if (0 != resELDV) 
		return 30 ;

	// load vo
	{
		FILE *fpVO = fopen(vo_filename.c_str(), "rb") ;
		if (NULL == fpVO) 
			return 40 ;
		// get file size
		fseek(fpVO, 0, SEEK_END) ;
		int32_t filesize = ftell(fpVO) ;
		fseek(fpVO, 0, SEEK_SET) ;
		char *buf = new char[filesize+1] ;
		if (NULL == buf) 
			{ fclose(fpVO) ; return 41 ; }
		int32_t L = fread(buf, 1, filesize, fpVO) ;
		fclose(fpVO) ;
		if (filesize != L) 
			{ delete [] buf ; return 42 ; }
		buf[filesize] = 0 ;
		int32_t OrderType = 0 ;
		int resLoadVO = p.LoadVariableOrderingFromBuffer(OrderType, 1, true, buf) ;
		delete [] buf ;
		if (0 != resLoadVO) 
			return 43 ;
	}

#ifdef PERFORM_SINGLTON_CONSISTENCY_CHECK
	int32_t nNewSingletonDomainVariables = 0, nVarsWithReducedDomainSize = 0 ;
	int res_SC = p.ComputeSingletonConsistency(nNewSingletonDomainVariables, nVarsWithReducedDomainSize) ;
	if (res_SC < 0) {
		// problem inconsistent
		}
#endif 

//	AndOrSearchSpace::AbsSamplingWorkspace ws(AbsSamplingTwoAndNodeCompare_Unique) ;
	AbsSamplingCompFn *fn = AbsSamplingTwoAndNodeCompare_Knuth ;
	bool alg_is_DFS = true ;
	bool abs_is_randomized = 0 == stricmp("rand", algorithm.c_str()) ;
	if (abs_is_randomized) {
		if (nAbs <= 0) nAbs = 1 ;
		}
	if (algorithm.length() > 0) {
		if (0 == stricmp("unique", algorithm.c_str())) 
			fn = AbsSamplingTwoAndNodeCompare_Unique ;
		else {
			if (! properAlg) {
				if (abs_is_randomized) 
					fn = AbsSamplingTwoAndNodeCompare_RandCntxt ;
				else 
					fn = AbsSamplingTwoAndNodeCompare_ContextNonProper ;
				}
			else {
				if (abs_is_randomized) 
					fn = AbsSamplingTwoAndNodeCompare_RandCntxt ;
				else if (0 == stricmp("customproper", algorithm.c_str()) || 0 == stricmp("context", algorithm.c_str())) 
					fn = alg_is_DFS ? AbsSamplingTwoAndNodeCompare_CustomProper_DFS : AbsSamplingTwoAndNodeCompare_CustomProper ;
				}
			}
		}
	AndOrSearchSpace::AbsSamplingWorkspace ws(fn) ;
	int resWSinit = ws.Initialize(p, true, NULL, false) ;
	if (0 != resWSinit) 
		return 50 ;
	ws.ForORtreetype_use_DFS_order() = forORtreetype_use_DFS_order ;

	if (nLevelsLimit < 0) nLevelsLimit = 0 ;
	ws.nLevelsLimit() = nLevelsLimit ;

	bool isAO = false;
	if(tree_type == "AO"){
		isAO = true;
	}
	int resCB = ws.CreateBuckets(isAO, true /* keep original bucket signatures; later when do MB processing, don't want to overwrite orig signature */, true, false) ;
	if (0 != resCB) 
		return 51 ;
//	int resCB = ws.CreateBuckets(true, true, false, isAO) ;

	// NOTE : the MaxNumVarsInBucket() is BE (not MBE) based!!! i.e. as if i-bound=inf
	if (p.VarOrdering_InducedWidth() < 0 && ws.MaxNumVarsInBucket() >= 0) {
		ws.SetVarOrdering_InducedWidth(ws.MaxNumVarsInBucket()-1) ;
		p.SetVarOrdering_InducedWidth(ws.MaxNumVarsInBucket()-1) ;
		}

	ws.MaxTreeHeight_BranchingVars_Limit() = nLevelsLimit ;
	int32_t res_MNBV = 
		0 
//		ws.Enfoce_MaxTreeHeight_BranchingVars_Limit() 
		;
	if (0 != res_MNBV) {
		return 52 ;
		}

	ws.MaxSpaceAllowed_Log10() = 9.0 ; // 1GB
	int iBoundMin = 2 ;
	int ib, nBP, maxDPB ; double spaceused ;
	if (iB < 0) {
		INT64 tIBfindS = ARE::GetTimeInMilliseconds() ;
		int resFindMinIB = ws.FindIBoundForSpaceAllowed(iBoundMin, ib, spaceused, nBP, maxDPB) ;
		INT64 tIBfindE = ARE::GetTimeInMilliseconds() ;
		if (ib < 0 || ib > ws.Problem()->N()) {
			printf("\nFindIBoundForSpaceAllowed failure; res==%d ib=%d", resFindMinIB, (int) iB) ;
			return 61 ;
			}
		iB = ib ;
		}
	else {
		ws.iBound() = iB ;
		int32_t resMBE = ws.CreateMBPartitioning(false, false, 0, 0) ;
		if (0 != resMBE) {
			printf("\nCreateMBPartitioning failure; res==%d", resMBE) ;
			return 62 ;
			}
		}

	// compute MBE induced_width; this takes into account actual i-bound
	int resMBEindw = ws.ComputeMaxNumVarsInBucket(true) ;
	if (0 != resMBEindw) {
		printf("\nComputeMaxNumVarsInBucket failure; res==%d", resMBEindw) ;
		return 63 ;
		}
	ws.SetPseudoWidth(ws.MaxNumVarsInBucket()-1) ;

	// DEBUGGGGG
	if (_DEBUG) {
		int32_t N_vars = 0, N_parents = 0, N_children = 0, nAF=0, nOF=0, nIF=0, nSIG=0 ;
		ARE::ARP *p_ = ws.Problem() ;
		int32_t N_ = p_->N() ;
		for (int32_t i = N_-1 ; i >= 0 ; --i) {
			int32_t v = p_->VarOrdering_VarList()[i] ;
			BucketElimination::Bucket *b = ws.MapVar2Bucket(v) ;
			if (NULL == b) continue ; ++N_vars ;
			if (NULL != b->ParentBucket()) N_parents++ ;
			N_children += b->nChildren() ;
			nOF += b->nOriginalFunctions() ;
			nAF += b->nAugmentedFunctions() ;
			nIF += b->nIntermediateFunctions() ;
			if (b->Width() >= 0) 
				nSIG += b->Width() ;
			b->VerifyIntermediateFunctions() ;
			b->VerifyAugmentedFunctions() ;
			}
		}

	signed char approx_bound = 1 ; // 1=max
	bool do_moment_matching = false ;
/*
	int resCOFs = ws.ComputeOutputFunctions(do_moment_matching, approx_bound) ;
	int resPCP = ws.PostComputationProcessing() ;
	double BEvalue = ws.CompleteEliminationResult() ;
	printf("\n MBE done; result=%g", BEvalue) ;
*/
	do_moment_matching = true ;
	int resCOFs_w = -1 ;
	try {
		resCOFs_w = ws.ComputeOutputFunctions(do_moment_matching, approx_bound) ;
		}
	catch (...) {
		// if MBE failed, because out of memory, exception may be thrown; we will land here then.
		int MBE_memory_error = 1 ;
		}
	if (0 != resCOFs_w) {
		printf("\nws.ComputeOutputFunctions() error; res=%d", resCOFs_w) ;
		}
	int resPCP_w = ws.PostComputationProcessing() ;
	double BEvalue_w = ws.CompleteEliminationResult(approx_bound) ;
	//printf("\nwMBE done; result=%g", BEvalue_w) ;

#if defined WINDOWS || _WINDOWS
	time_t wall_preproc_done ;
	time(&wall_preproc_done) ;
	int64_t t_elapsed_preproc_done = wall_preproc_done - wall0 ;
#else
    double wall_preproc_done = get_wall_time();
	int64_t t_elapsed_preproc_done = (int64_t) (wall_preproc_done - wall0) ;
#endif // WINDOWS

	std::string output_filename;
	FILE* output_file;

	std::string snContext = std::to_string(nContext);
	std::string snAbs = std::to_string(nAbs);
	std::string sabs_indv = abs_indv ? "1" : "0" ;
	if (nLevelsLimit >= 0 && nLevelsLimit < 9000){
		snContext = std::to_string(nContext) + "_" + std::to_string(nLevelsLimit);
		}
    output_filename = abs_is_randomized ? 
		(problem_filename + "-" + tree_type + "-i-" + std::to_string(iB) + "-prop-" + (properAlg ? "1" : "0") + "-a-" + algorithm + "-nC-"  + snContext + "-nAbs-"  + snAbs + "-abs_indv-" + sabs_indv + "-nR-" +  std::to_string(nrunstodo) +".out") 
	    :
		(problem_filename + "-" + tree_type + "-i-" + std::to_string(iB) + "-prop-" + (properAlg ? "1" : "0") + "-a-" + algorithm + "-nC-"  + snContext + "-nR-" +  std::to_string(nrunstodo) +".out") ;
    output_file = fopen(output_filename.c_str(), "w");
    fprintf(output_file, "%s \t%d \t%s \t%s \t%d \t%d \t%d \t%d \t%d \t%g \t%g \t%lld \t%c \t%s\n", 
		problem_filename.c_str(), (int) iB, algorithm.c_str(), snContext.c_str(), (int) nrunstodo, (int) (AS_VERSION), 
		(int) ws.Problem()->N(), (int) (ws.MaxNumVarsInBucket()-1), (int) ws.MaxTreeHeight(), 
		(double) BEvalue_w, (double) exactZ, t_elapsed_preproc_done, exclude_prep_time ? 'Y' : 'N', tree_type.c_str()) ;
	fflush(output_file);
	if (exclude_prep_time)
		wall0 = wall_preproc_done;

 	// printf("\nBEWS init done; N=%d nBuckets=%d var_order width=%d mbe width=%d nVarsPartitioned=%d", (int) p.N(), (int) ws.nBuckets(), (int) (ws.VarOrdering_InducedWidth()), (int) (ws.PseudoWidth()), (int) ws.nBucketsWithPartitioning()) ;
	// printf("\nBEWS init done; N=%d nBuckets=%d width=%d nVarsPartitioned=%d", (int) p.N(), (int) ws.nBuckets(), (int) (ws.MaxNumVarsInBucket()-1), (int) ws.nBucketsWithPartitioning()) ;
	printf("%s \t%d \t%d \t%d \t%d \t%d \t%d", problem_filename.c_str(), (int) p.N(), (int) ws.nBuckets(), (int) (ws.VarOrdering_InducedWidth()), (int) ws.nBucketsWithPartitioning(), (int) ws.Height(), (int) (AS_VERSION)) ;

	// log stuff
	if (false) {
		ARE::ARP *p = ws.Problem() ;
		FILE *fp_log = fopen("OR-log.txt", "w");
		fprintf(fp_log, "problem=%s iB=%d", problem_filename.c_str(), (int) iB);
		for (int32_t i = 0 ; i < p->N() ; ++i) {
			BucketElimination::Bucket *b = ws.MapVar2Bucket(i) ;
			fprintf(fp_log, "\n var=%d nMBs=%d", i, (int) b->MiniBuckets().size());
			}
		fclose(fp_log) ;
		}

	// Store pseudotree derived from Bucket Tree
 	// std::string pseudotree_filename;
 	// FILE* pseudotree_file;
	// pseudotree_filename = "Pseudotree/" + problem_filename + ".pt";
 	// pseudotree_file = fopen(pseudotree_filename.c_str(), "w");
	// int32_t * ptArray;
	// ws.ExportPseudotree(ptArray);
	// for(int32_t i = 0; i < p.N(); i++){
	// 	fprintf(pseudotree_file, "%d, %d\n", i, ptArray[i]);
	// }
	// delete[] ptArray;	
	
#if defined WINDOWS || _WINDOWS
#else
    wall1 = get_wall_time();
    cpu1  = get_cpu_time();
    double wall_prep = wall1 - wall0;
    double cpu_prep= cpu1  - cpu0;
    //printf("\nBEresult=%g iB=%d", BEvalue, ws.iBound()) ;
    printf(" \t%g \t%d \t%s \t%g \t%g", BEvalue_w, ws.iBound(), hybrid_al.c_str(), wall_prep, cpu_prep);
#endif // WINDOWS

	// generate full AND/OR search tree
	bool ignoreclosestbranchingvariable = properAlg ? false : true ;
	int32_t resPREP = ws.PrepSampledSearchTree(nContext, ignoreclosestbranchingvariable) ;
	if (0 != resPREP) {
		printf("\nFAILED ABS SEARCH TREE PREP; res=%d ...", resPREP) ;
		}

	if (abs_is_randomized) {
		int res_RAsetup = ws.ComputeRandAbstractionFactors(nAbs) ;
		if (0 != res_RAsetup) {
			printf("\nFAILED ComputeRandAbstractionFactors(%d) ...", nAbs) ;
			}
		}

	int32_t nTriesGood = 0 ; int64_t sumTreeSize = 0, nNodesMerged = 0, nNodesCreated = 0;
	double sumOfTries;
	double varianceTerm;
	double varianceSum;
	double variance;
	double effective_variance;	

    double wall_inter, cpu_inter;

    int32_t runsdone = 0;
    double old_time = 0;
    double new_time = 0;
    bool do_print = false;
	int32_t dDeepestBucketExpanded = 0, nDFSBranchingPointsProcessed = 0 ;
	for (int32_t i = 0 ; i < nrunstodo ; i++) {
		runsdone++;
		if (ws.nNodesInTree() > 1) 
			ws.DestroyTree() ;
#ifdef doBFS
		int resGENtree = ws.GenerateSampledSearchTree() ;
//		int resGENtree = ws.GenerateFullSearchTree() ;
		if (verbose) {
			int32_t nN = ws.ComputeTreeSize() ;
			printf("\nAStree built res=%d nNodes=%d(%d) nNodesCreated=%d nNodesISmerged=%d nNodesPartitioned=%d", resGENtree, (int) ws.nNodesInTree(), nN, (int) ws.nNodesCreated(), (int) ws.nNodesISmerged(), (int) ws.nVarsPartitioned()) ;
			}

		// compute tree value
		int resTV = ws.ComputeSearchTreeValue() ;
#else
		if (abs_is_randomized && abs_indv) {
			int res_RAC_setup = ws.ComputeRandAbstractionFactors_Indv() ;
			}
		int32_t resTV = properAlg ? ws.RunDFS() : ws.RunDFS_nonProper(dDeepestBucketExpanded, nDFSBranchingPointsProcessed) ;
#endif // doBFS
		double SAMPLvalue = ws.Root().SubTreeValue() ;
		int isnan_res = isnan(SAMPLvalue) ;
		// if (0 != isnan_res) 
		// 	printf("\nAStree value isNaN; i=%d", i) ;
		if (verbose) 
			printf("\nAStree value res=%d SAMPLvalue=%g dDeepestBucketExpanded=%d nDFSBranchingPointsProcessed=%d", resTV, SAMPLvalue, dDeepestBucketExpanded, nDFSBranchingPointsProcessed) ;

		if (0 == resTV && 0 == isnan_res) { // ? -std::numeric_limits<float>::infinity() != SAMPLvalue : false) {
			// Compute variance term
			if (SAMPLvalue == -std::numeric_limits<double>::infinity()){
				varianceTerm = 2*exactZ;
			}
			else {
				varianceTerm = 2*logabsdiffexp10(SAMPLvalue, exactZ);
			}

			if (0 == nTriesGood++) {
		        sumOfTries = SAMPLvalue ;
		    	varianceSum = varianceTerm ;
		    }
		    else { 
		        LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(sumOfTries,sumOfTries,SAMPLvalue) ;
		        LOG_OF_SUM_OF_TWO_NUMBERS_GIVEN_AS_LOGS(varianceSum, varianceSum, varianceTerm) ;
		    }
		}
		else {
			bool is_inf = true ;
		}

		sumTreeSize += ws.nNodesInTree() ;
		nNodesMerged += ws.nNodesISmerged() ;
		nNodesCreated += ws.nNodesCreated() ;

#if defined WINDOWS || _WINDOWS
		time_t wall1 ;
		time(&wall1) ;
//		wall_inter = (wall1 - wall0)/1000.0 ;
		wall_inter = (wall1 - wall0) ; // time() returns elapsed time since 1970 in seconds
		cpu_inter = wall_inter ;
#else
        wall1 = get_wall_time();
        cpu1  = get_cpu_time();
        wall_inter = wall1 - wall0;
        cpu_inter = cpu1 - cpu0;
#endif // WINDOWS
        if(runsdone <=1){
        	old_time =  new_time = cpu_inter;
        }
        
        new_time = cpu_inter;

        if((new_time - old_time) > 1.0){
        	old_time = new_time;
        	do_print = true;
        }
        else{
        	do_print = false;
        }
        
		bool need_break = false ;
		if(cpu_inter >= time_limit) 
			need_break = true ;

        if (do_print || runsdone == 1 || need_break) { //runsdone%print_period) == 0 || print_period == 1  ||  runsdone == 1
			variance = varianceSum - log10((double) nTriesGood);
			effective_variance = variance + log10((double) nNodesCreated) - log10((double) runsdone);
        	//variance = nTriesGood > 0 ? varianceSum - log10((double)nTriesGood): -std::numeric_limits<double>::infinity();
			//effective_variance = nTriesGood > 0 ? variance + log10((double) nNodesCreated) - log10((double)runsdone) : -std::numeric_limits<double>::infinity();
		
        	fprintf(output_file, "%d \t%g \t%g \t%lld \t%lld \t%lld \t%lld \t%lld \t%g \t%g \t%lld \t%lld \t%g \t%g\n", 
				(int32_t) i+1, SAMPLvalue, (double) (nTriesGood > 0 ? sumOfTries - log10((double)nTriesGood): -std::numeric_limits<float>::infinity()),
        		(int32_t) nTriesGood, (int64_t) ws.nNodesInTree(), (int64_t) ws.nNodesISmerged(), (int64_t) (sumTreeSize/runsdone), (int64_t) (nNodesMerged/runsdone),
        		(double)wall_inter, (double)cpu_inter, (int64_t) ws.nNodesCreated(), (int64_t) (nNodesCreated/runsdone), (double)variance, (double)effective_variance);
			fflush(output_file);
    		}

    	/*if(do_print || runsdone == 1){ //runsdone%print_period) == 0 || print_period == 1  ||  runsdone == 1
		variance = nTriesGood > 0 ? varianceSum - log10((double)nTriesGood): -std::numeric_limits<float>::infinity();
		effective_variance = nTriesGood > 0 ? variance + log10((double) nNodesCreated) - log10((double)runsdone) : -std::numeric_limits<float>::infinity();
		fprintf(output_file, "%d \t%g \t%g \t%d \t%lld \t%lld \t%lld \t%lld \t%g \t%g \t%lld \t%lld t%g \t%g\n", (int) i+1, SAMPLvalue, nTriesGood > 0 ? sumOfTries - log10((double)nTriesGood): -std::numeric_limits<float>::infinity(),
        		(int) nTriesGood, ws.nNodesInTree(), ws.nNodesISmerged(), (int64_t) (sumTreeSize/runsdone), (int64_t) (nNodesMerged/runsdone),
        		wall_inter, cpu_inter, ws.nNodesCreated(), (int64_t) (nNodesCreated/runsdone), , variance, effective_variance);
    	}*/
		
        if (need_break) 
			break ;
		}

//  Stop timers
#if defined WINDOWS || _WINDOWS
	time_t wall1 ;
	time(&wall1) ;
    double wall_final = (wall1 - wall0)/1000.0 ;
	double cpu_final = wall_final ;
#else
    wall1 = get_wall_time();
    cpu1  = get_cpu_time();

    double wall_final = wall1 - wall0;
    double cpu_final = cpu1  - cpu0;
#endif // WINDOWS

 	// printf("\nMaxNBVC: %d\n", ws.MaxNumBranchingVariablesInChain());
 	// printf("NBV: %d\n", ws.NumBranchingVariables());
	// printf("Estimated STS: %lld\n", ws.EstimateSampleTreeSize());

	/*
	variance = nTriesGood > 0 ? varianceSum - log10((double)nTriesGood): -std::numeric_limits<float>::infinity();
    effective_variance = nTriesGood > 0 ? variance + log10((double) nNodesCreated) - log10((double)runsdone) : -std::numeric_limits<float>::infinity();																									   
	if (runsdone > 0) {
		//printf("\nASresultAvg=%g (%d tries; nNodes=%lld nNodesISmerged=%lld)", sumOfTries/((double)nTriesGood), (int) nTriesGood, (int64_t) (sumTreeSize/nTriesGood), (int64_t) (nNodesMerged/nTriesGood)) ;
   		printf(" \t%g \t%d \t%lld \t%lld \t%g \t%g \t%s \t%d \t%d \t%lld \t%d \t%d \t%lld \t%d \t%lf \t%lf \t%g \t%g \t%d\n", nTriesGood > 0 ? sumOfTries - log10((double)nTriesGood): -std::numeric_limits<float>::infinity(), (int) nTriesGood, 
   			   (int64_t) (sumTreeSize/runsdone), (int64_t) (nNodesMerged/runsdone), wall_final, cpu_final, output_filename.c_str(),
   			   ws.MaxNumBranchingVariablesInChain(), ws.NumBranchingVariables(), ws.EstimateSampleTreeSize(), (int) (ws.VarOrdering_InducedWidth()), (int) (ws.PseudoWidth()), (int64_t) (nNodesCreated/runsdone), runsdone, heuristicCoefficient, heuristicPower, 0.0, 0.0, ws.nNodesExactHComputed());																				
	}
    */

	//printf("%d\n", nTriesGood );

	//variance = nTriesGood > 0 ? varianceSum - log10((double)nTriesGood): -std::numeric_limits<double>::infinity();
    //effective_variance = nTriesGood > 0 ? variance + log10((double) nNodesCreated) - log10((double)runsdone) : -std::numeric_limits<double>::infinity();
	//printf("\n%g %g\n", variance, effective_variance);

	if (runsdone > 0) {
		//printf("\nASresultAvg=%g (%d tries; nNodes=%lld nNodesISmerged=%lld)", sumOfTries/((double)nTriesGood), (int) nTriesGood, (int64_t) (sumTreeSize/nTriesGood), (int64_t) (nNodesMerged/nTriesGood)) ;
   		printf(" \t%g \t%d \t%lld \t%lld \t%g \t%g \t%s \t%d \t%d \t%lld \t%d \t%d \t%lld \t%d \t%g \t%g \t%d\n", nTriesGood > 0 ? sumOfTries - log10((double)nTriesGood): -std::numeric_limits<float>::infinity(), (int) nTriesGood, 
   			   (int64_t) (sumTreeSize/runsdone), (int64_t) (nNodesMerged/runsdone), wall_final, cpu_final, output_filename.c_str(),
   			   ws.MaxNumBranchingVariablesInChain(), ws.NumBranchingVariables(), ws.EstimateSampleTreeSize(), (int) (ws.VarOrdering_InducedWidth()), (int) (ws.PseudoWidth()), (int64_t) (nNodesCreated/runsdone), runsdone, variance, effective_variance, ws.nNodesExactHComputed()) ;
	}

    fclose(output_file);
    // fclose(pseudotree_file);
    return 0 ;
}

