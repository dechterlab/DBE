#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <thread>
#include <time.h>
#include <cmath>
#include <limits>

#include <iostream>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#else
#include <sys/time.h>

#endif // WINDOWS

#include "CVO/VariableOrderComputation.hxx"
#include "BE/Bucket.hxx"
#include "BE/MBEworkspace.hxx"

#ifdef _MSC_VER
#ifndef _WIN64
// error - this project is meant as 64-bit...
#endif
#ifndef _M_X64
// error - this project is meant as 64-bit...
#endif
#endif

//#define PERFORM_SINGLTON_CONSISTENCY_CHECK

#if defined WINDOWS || _WINDOWS
#else
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
ARE::FnConstructor deepMBE_mboutputfncnstrctr(BucketElimination::MiniBucket & MB)
{
    // TODO : implement logic to decide what kind of output fn to use for the given MB...
//    if(MB->Width() == 20)
//            return ARE::FunctionNNConstructor;
    return ARE::FunctionConstructor ;
}
#include <torch/torch.h>
double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif // WINDOWS

#include <Config.h>

Config_NN global_config;

static MTRand RNG ;

double logabsdiffexp10(double num1, double num2) {
    if (num1 == num2) {
        return -std::numeric_limits<double>::infinity();
    }
    num1 = num1 * std::log(10.0);
    num2 = num2 * std::log(10.0);
    double num_max = num1 < num2 ? num2 : num1 ; // std::max(num1, num2);
    double num_min = num1 > num2 ? num2 : num1 ; // std::min(num1, num2);
    double ans = (num_min + std::log(expm1(num_max - num_min))) / std::log(10.0) ;
    // std::cout << ans << std::endl;
    return ans;
}

int32_t main(int32_t argc, char* argv[])
{
    torch::Tensor tensor2 = torch::eye(3);
    std::cout << tensor2 << std::endl;


#if defined WINDOWS || _WINDOWS
    time_t wall0 ;
	time(&wall0) ;
#else
    double wall0 = get_wall_time();
    double cpu0 = get_cpu_time();
    double wall1, cpu1;
#endif // WINDOWS

    int32_t nParams = argc ;
    unsigned long randomGeneratorSeed = 0 ;
    std::string problem_filename, evidence_filename, vo_filename ;
    char* p1, p2;
    int32_t nArgs = (nParams - 1) >> 1 ;
    int32_t time_limit = 60;
    int32_t print_period = 1;
    bool findPracticalVariableOrder = true ;
    bool verbose = false ;
    int32_t iB = -1 ;
    ARE::VarElimOrderComp::ObjectiveToMinimize objCodePrimary = ARE::VarElimOrderComp::Width, objCodeSecondary = ARE::VarElimOrderComp::None ;
    std::string algorithm ;
    double exactZ = 0.0;
    int32_t v2sample = -1, nsamples = -1 ;

    if (1 + 2 * nArgs != nParams) {
        printf("\nBAD COMMAND LINE; will exit ...") ;
        return 1 ;
    }

    for (int32_t i = 0 ; i < nArgs; i++) {
        std::string sArgID = (NULL != argv[1 + 2 * i] ? argv[1 + 2 * i] : "") ;
        std::string sArg = (NULL != argv[2 * (i + 1)] ? argv[2 * (i + 1)] : "") ;
        if (0 == stricmp("-s", sArgID.c_str()))
            randomGeneratorSeed = std::strtoul(sArg.c_str(), NULL, 0) ;
        else if (0 == stricmp("-fUAI", sArgID.c_str()))
            problem_filename = sArg ;
        else if (0 == stricmp("-fEV", sArgID.c_str()))
            evidence_filename = sArg ;
        else if (0 == stricmp("-fVO", sArgID.c_str()))
            vo_filename = sArg ;
        else if (0 == stricmp("-printPeriod", sArgID.c_str()))
            print_period = atoi(sArg.c_str());
        else if (0 == stricmp("-t", sArgID.c_str()))
            time_limit = atoi(sArg.c_str());
        else if (0 == stricmp("-iB", sArgID.c_str()))
            iB = atoi(sArg.c_str());
        else if (0 == stricmp("-Verbose", sArgID.c_str()))
            verbose = atoi(sArg.c_str()) > 0 ;
        else if (0 == stricmp("-v2sample", sArgID.c_str()))
            v2sample = atoi(sArg.c_str()) ;
        else if (0 == stricmp("-nsamples", sArgID.c_str()))
            nsamples = atoi(sArg.c_str()) ;
        else if (0 == stricmp("-a", sArgID.c_str()))
            algorithm = sArg ;
        else if (0 == stricmp("-fpvo", sArgID.c_str()))
            findPracticalVariableOrder = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0] ;
        else if (0 == stricmp("-O1", sArgID.c_str()))
            objCodePrimary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
        else if (0 == stricmp("-O2", sArgID.c_str()))
            objCodeSecondary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
        else if (0 == stricmp("-exactZ", sArgID.c_str()))
            exactZ = atof(sArg.c_str());
        else if (0 == stricmp("-h_dim",sArgID.c_str() ))
            global_config.h_dim = atoi(sArg.c_str());
        else if (0 == stricmp("-batch_size",sArgID.c_str() ))
            global_config.batch_size = atoi(sArg.c_str());
        else if (0 == stricmp("-n_epochs",sArgID.c_str() ))
            global_config.n_epochs = atoi(sArg.c_str());
        else if (0 == stricmp("-lr",sArgID.c_str() ))
            global_config.lr = atof(sArg.c_str());
        else if (0 == stricmp("-sample_size",sArgID.c_str() ))
            global_config.sample_size = atoi(sArg.c_str());
        else if (0 == stricmp("--out_file",sArgID.c_str() ))
            global_config.out_file = sArg;
        else if (0 == stricmp("--out_file2",sArgID.c_str() ))
            global_config.out_file2 = sArg;
        else if (0 == stricmp("--network",sArgID.c_str() ))
            global_config.network = sArg;
        else if (0 == stricmp("--sampling_method",sArgID.c_str() ))
            global_config.s_method = sArg;
        else if (0 == stricmp("--width_problem",sArgID.c_str() ))
            global_config.width_problem = atoi(sArg.c_str());
        else if (0 == stricmp("--loss_weight", sArgID.c_str()))
            global_config.loss_weight = atof(sArg.c_str());
        else if (0 == stricmp("--loss_weight_mse", sArgID.c_str()))
            global_config.loss_weight_mse = atof(sArg.c_str());
        else if (0 == stricmp("--train_stop", sArgID.c_str()))
            global_config.train_stop = sArg;
        else if (0 == stricmp("--stop_iter", sArgID.c_str()))
            global_config.stop_iter = atoi(sArg.c_str());
    }


    int32_t maxNumProcessorThreads = std::thread::hardware_concurrency() ; // this requires c++11
    if (maxNumProcessorThreads <= 0)
        maxNumProcessorThreads = 1 ;

    if (0 == problem_filename.length() || 0 == vo_filename.length())
        return 1 ;

    ARE::ARP p("mbe") ;
    p.SetOperators(FN_COBINATION_TYPE_PROD, VAR_ELIMINATION_TYPE_SUM) ;
    int32_t resLoad = p.LoadFromFile(problem_filename) ;
    if (0 != resLoad)
        return 10 ;
    int32_t resPostCA = p.PerformPostConstructionAnalysis() ;
    if (0 != resPostCA)
        return 11 ;
    if (iB > p.N())
        iB = p.N() ;

    int32_t nEV = -1 ;
    if (evidence_filename.length() > 0) {
        int32_t resLoadE = p.LoadFromFile_Evidence(evidence_filename, nEV) ;
        if (0 != resLoadE)
            return 20 ;
        int32_t resElimE = p.EliminateEvidence() ;
        if (0 != resElimE)
            return 21 ;
    }

    int32_t resELDV = p.EliminateSingletonDomainVariables() ;
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
        char *buf = new char[filesize + 1] ;
        if (NULL == buf)
        {
            fclose(fpVO) ; return 41 ;
        }
        int32_t L = fread(buf, 1, filesize, fpVO) ;
        fclose(fpVO) ;
        if (filesize != L)
        {
            delete[] buf ; return 42 ;
        }
        buf[filesize] = 0 ;
        int32_t OrderType = 0 ; // 0=OrderType means variable elimination order. 1=OrderType means variable bucket-tree order.
        int32_t resLoadVO = p.LoadVariableOrderingFromBuffer(OrderType, 1, true, buf) ;
        delete[] buf ;
        if (0 != resLoadVO)
            return 43 ;
    }
    int32_t k_max=0;

    for(int32_t i=0; i<p.N();i++)
    {
        int32_t k = p.K(i);
        if (k>k_max)
            k_max=k;
    }

#ifdef PERFORM_SINGLTON_CONSISTENCY_CHECK
    int32_t nNewSingletonDomainVariables = 0, nVarsWithReducedDomainSize = 0 ;
	int32_t res_SC = p.ComputeSingletonConsistency(nNewSingletonDomainVariables, nVarsWithReducedDomainSize) ;
	if (res_SC < 0) {
		// problem inconsistent
	}
#endif
    BucketElimination::MBEworkspace ws ;
    int32_t delete_used_tables = 0 ; // keep tables in  memory so that we can sample from them...
    int32_t resWSinit = ws.Initialize(p, true, NULL, delete_used_tables) ;
    if (0 != resWSinit)
        return 50 ;
    bool isAO = true ;
    int32_t resCB = ws.CreateBuckets(isAO, true /* keep original bucket signatures; later when do MB processing, don't want to overwrite orig signature */, true, false) ;
    if (0 != resCB)
        return 51 ;

    // NOTE : the MaxNumVarsInBucket() is BE (not MBE) based!!! i.e. as if i-bound=inf
    if (p.VarOrdering_InducedWidth() < 0 && ws.MaxNumVarsInBucket() >= 0) {
        ws.SetVarOrdering_InducedWidth(ws.MaxNumVarsInBucket() - 1) ;
        p.SetVarOrdering_InducedWidth(ws.MaxNumVarsInBucket() - 1) ;
    }

    ws.MaxSpaceAllowed_Log10() = 9.0 ; // 1GB
    int32_t iBoundMin = 2 ;
    int32_t ib, nBP, maxDPB ; double spaceused ;
    if (iB < 0) {
        INT64 tIBfindS = ARE::GetTimeInMilliseconds() ;
        int32_t resFindMinIB = ws.FindIBoundForSpaceAllowed(iBoundMin, ib, spaceused, nBP, maxDPB) ;
        INT64 tIBfindE = ARE::GetTimeInMilliseconds() ;
        if (ib < 0 || ib > ws.Problem()->N()) {
            printf("\nFindIBoundForSpaceAllowed failure; res==%d ib=%d", resFindMinIB, (int)iB) ;
            return 61 ;
        } else
            printf("\n FindIboundForSpaceAllowed %d", (int)ib);
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
    int32_t resMBEindw = ws.ComputeMaxNumVarsInBucket(true) ;
    if (0 != resMBEindw) {
        printf("\nComputeMaxNumVarsInBucket failure; res==%d", resMBEindw) ;
        return 63 ;
    }
    ws.SetPseudoWidth(ws.MaxNumVarsInBucket() - 1) ;

    // print some stats
    printf("*************************");
    printf("\nSTATS : nBuckets=%d MaxNumVarsInBucket=%d nPartitionedBuckets=%d", (int) ws.nBuckets(), (int) ws.MaxNumVarsInBucket(), (int) ws.nBucketsWithPartitioning()) ;

    signed char approx_bound = 1 ; // 1=max
    bool do_moment_matching = true ;
    int32_t resCOFs_w = -1 ;
    time_t tStart = 0 ; time(&tStart) ;
    int64_t totalSumOutputFunctionsNumEntries = 0 ;
    try {
        resCOFs_w = ws.ComputeOutputFunctions(do_moment_matching, approx_bound, totalSumOutputFunctionsNumEntries) ;
    }
    catch (...) {
        // if MBE failed, because out of memory, exception may be thrown; we will land here then.
        int32_t MBE_memory_error = 1 ;
    }
    time_t tEnd = 0 ; time(&tEnd) ;
    int64_t tElapsed = tEnd - tStart ;
    if (0 != resCOFs_w) {
        printf("\nws.ComputeOutputFunctions() error; res=%d", resCOFs_w) ;
    }
    int32_t resPCP_w = ws.PostComputationProcessing() ;
    double BEvalue_w = ws.CompleteEliminationResult(approx_bound) ;
    printf("\nwMBE done; result=%g runtime=%lldsec tablesmemory=%lld bytes", BEvalue_w, tElapsed, (int64_t) sizeof(double)*totalSumOutputFunctionsNumEntries) ;


    std::string o_file = global_config.out_file + ".txt";
    printf("Writing to -- %s", o_file.c_str());

    std::ofstream to_write;
    to_write.open(o_file,std::ios_base::app);

    if (to_write.is_open())
    {
        to_write <<  BEvalue_w  << std::endl ;
        to_write.close();
    }

    std::string o_file2 = global_config.out_file2 + ".txt";

    std::ofstream f;
    f.open(o_file2,std::ios_base::app);

    if (f.is_open())
    {
        f << (int)iB << '\t' << k_max << '\t' << global_config.width_problem << '\t' << (int) ws.nBuckets()  << '\t' << (int) ws.count_ComputeOutputFunction_NN << '\t' << global_config.avg_val_mse << '\t' << global_config.avg_test_mse << '\t' << global_config.avg_val_w_mse << '\t' << global_config.avg_test_w_mse << '\t'<< BEvalue_w << '\t' << (float) tElapsed/3600 << std::endl ;
        f.close();
    }

    BucketElimination::Bucket *B = ws.GetBucketWithMostVariables() ;
    int32_t V = B->V() ;
    printf("\nlargest bucket v=%d", (int32_t) V) ;
    if (v2sample < 0)
        v2sample = V ;
    if (nsamples <= 0)
        nsamples = 1 ;
    BucketElimination::Bucket *B_ = ws.MapVar2Bucket(v2sample) ;
    if (NULL == B_) {
        printf("\nERROR : B_ == NULL") ;
        exit(1) ;
    }
    if (B_->V() != v2sample) {
        printf("\nERROR : B_->V() != v2sample") ;
        exit(1) ;
    }
    std::vector<BucketElimination::MiniBucket *> & MBs = B_->MiniBuckets() ;
    printf("\nwill sample v=%d nsamples=%d bSize=%d vars nMBs=%d", (int32_t) v2sample, (int32_t) nsamples, (int32_t) B_->Width(), (int32_t) MBs.size()) ;
    if (MBs.size() <= 0) {
        printf("\nERROR : nMBS <= 0") ;
        exit(1) ;
    }
    ARE::Function & output_fn = MBs[0]->OutputFunction() ;
    printf("\noutput fn size=%lld entries\n", (int64_t) output_fn.TableSize()) ;
    printf(" ComputeOutputFunction_NN Time: %f Count: %d Average: %f\n", ws.time_ComputeOutputFunction_NN, ws.count_ComputeOutputFunction_NN, ws.time_ComputeOutputFunction_NN/ws.count_ComputeOutputFunction_NN);
    printf(" TableEntryEx Time: %f Count: %d Average: %f\n", ws.time_TableEntryEx, ws.count_TableEntryEx, ws.time_TableEntryEx / ws.count_TableEntryEx);
    printf(" Train Time: %f Count: %d Average: %f\n", ws.time_Train, ws.count_Train, ws.time_Train / ws.count_Train);
    // samples randomly
//	std::vector<int32_t> signature ; signature.reserve(output_fn.N()) ; signature.resize(output_fn.N(), -1) ;
//	for (int32_t i = 0 ; i < output_fn.N() ; ++i) signature[i] = output_fn.Argument(i) ;
//	output_fn.ComputeArgumentsPermutationList(signature.size(), signature.data()) ;
    std::vector<int32_t> signature_assignment ; signature_assignment.reserve(output_fn.N()) ; signature_assignment.resize(output_fn.N(), -1) ;
    printf("\nfn signature : ") ;
    for (int32_t j = 0 ; j < output_fn.N() ; ++j) {
        if (0 == j) printf("%d", (int32_t) output_fn.Argument(j)) ;
        else printf(" %d", (int32_t) output_fn.Argument(j)) ;
    }
    for (int32_t i = 0 ; i < nsamples ; ++i) {
        // generate signature assignment
        for (int32_t j = 0 ; j < output_fn.N() ; ++j) {
            int32_t v = output_fn.Argument(j) ;
            int32_t domain_size_of_v = p.K(v) ;
            int32_t value = RNG.randInt(domain_size_of_v-1) ;
            signature_assignment[j] = value ;
        }
        // fetch value
        int64_t fn_idx = output_fn.ComputeFnTableAdr_wrtSignatureAssignment(signature_assignment.data(), p.K()) ;
        if (fn_idx < 0 || fn_idx >= output_fn.TableSize()) {
            printf("\nERROR : fn_idx out of bounds...") ;
            continue ;
        }
        double fn_value = output_fn.TableEntry(fn_idx) ;
        // DEBUG : print sample
        printf("\nsample : ") ;
        for (int32_t j = 0 ; j < output_fn.N() ; ++j) {
            if (0 == j) printf("%d", (int32_t) signature_assignment[j]) ;
            else printf(" %d", (int32_t) signature_assignment[j]) ;
        }
        printf(" = %g logscale (idx=%lld)", fn_value, fn_idx) ;
    }

    return 0 ;
}
