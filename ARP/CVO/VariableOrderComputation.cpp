// CVO.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <sstream>
#include <thread>
#include <time.h>

#if defined WINDOWS || _WINDOWS
#include <process.h>    /* _beginthread, _endthread */
#endif // WINDOWS

#include "Utils/AVLtreeSimple.hxx"
#include "Utils/RandomProblemGenerator.hxx"
#include "Utils/MiscUtils.hxx"
#include "CVO/VariableOrderComputation.hxx"
#include "BE/Bucket.hxx"
#include "BE/MBEworkspace.hxx"

#define SAVE_ORDER
//#define VERBOSE_CVO
//#define RUN_FOLDER_PROBLEMS
//#define PRINT_BEST_TD_AT_END

#ifdef LINUX
#include <signal.h>
#include <pthread.h>
#include <dirent.h>
#include <sys/stat.h>
pthread_mutex_t nRunsSumMutex = PTHREAD_MUTEX_INITIALIZER ;
pthread_mutex_t stopSignalMutex = PTHREAD_MUTEX_INITIALIZER ;
pthread_mutex_t printTDMutex = PTHREAD_MUTEX_INITIALIZER ;
#endif

#ifdef WINDOWS
std::string wstringToString(const std::wstring& wString)
{
	unsigned int size = WideCharToMultiByte(CP_UTF8, 0, wString.c_str(), wString.length(), NULL, 0, NULL, NULL) ;
	if (0 == size)
		return ""; // error
	char* utf8 = new char[size + 1];
	::WideCharToMultiByte(CP_UTF8, 0, wString.c_str(), wString.length(), utf8, size+1, NULL, NULL) ;
	utf8[size] = '\0';
	std::string result(utf8, size);
	delete[] utf8;
	return result ;
}
#endif

static int ReadDirectory(const char *Dir, const char *Pattern, const char *prefix, const char *ext, std::vector<std::string> & Filenames)
{
#if defined _WINDOWS || WINDOWS
    HANDLE hFile ;
	WIN32_FIND_DATAA FileInformation ;
	std::string strPattern(Dir) ;
	if (strPattern.length() > 0) {
		if ('\\' != strPattern[strPattern.length()-1] && L'/' != strPattern[strPattern.length()-1])
			strPattern += '\\' ;
		}
	if (NULL != Pattern) {
		strPattern += Pattern ;
		}
	else {
		if (NULL != prefix)
			strPattern += prefix ;
		strPattern += "*." ;
		if (NULL != ext)
			strPattern += ext ;
		}
	hFile = ::FindFirstFileA(strPattern.c_str(), &FileInformation) ;
	if (INVALID_HANDLE_VALUE != hFile) {
		do {
			bool isdot = '.' == FileInformation.cFileName[0] ;
			bool isdir = 0 != (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ;
			if (! isdot && ! isdir) {
				Filenames.push_back(std::string(FileInformation.cFileName)) ;
				}
			} while(TRUE == ::FindNextFileA(hFile, &FileInformation)) ;
		::FindClose(hFile) ;
		}
#else

    std::string extStr ; if (ext != NULL) extStr = ext ;
    if(extStr.length() > 0) {
        if(extStr[0] != '.')
            extStr.insert(0, ".");
        std::transform(extStr.begin(), extStr.end(), extStr.begin(), ::tolower);  // assume ascii
	    }

    DIR *pDir = NULL;
    dirent *pEnt = NULL;
		std::string PathSeparator = "/";
		char wPathSeparator = '/';
    struct stat st;

    std::string sDir(Dir);
    if(sDir[sDir.length() - 1] != wPathSeparator)
        sDir += wPathSeparator;

    pDir = opendir(sDir.c_str());
    if(NULL == pDir)
        return 1;

    while ((pEnt = readdir(pDir)) != NULL) {
        std::string fName = pEnt->d_name ;
        if (fName[0] == '.')
            continue ;
        std::string fullFileName = sDir + PathSeparator + fName ;
        if (stat(fullFileName.c_str(), &st) == -1)
            continue ;

        const bool isDirectory = (st.st_mode & S_IFDIR) != 0 ;
        if (isDirectory)
            continue ;

        // TODO: Implement pattern matching, prefix, etc. For now, check extension only.
        if(extStr.length() > 0) {
            if(fName.length() <= extStr.length())
                continue;
            std::string fNameEnd = fName.substr((int)fName.size() - (int)extStr.size(), (int)extStr.size());
            std::transform(fNameEnd.begin(), fNameEnd.end(), fNameEnd.begin(), ::tolower);  // assume ascii
            if(fNameEnd != extStr)
                continue;
			}
        // Filenames.push_back(stringToWstring(fName)) ;
				Filenames.push_back(fName) ;
		}
    closedir(pDir) ;
#endif

    return 0 ;
}


static void GetCurrentDTsec(char *strDT, time_t & ttNow)
{
	// get current date/time string so that we can log when stuff is happening
	if (0 == ttNow) time(&ttNow) ;
	struct tm *pTime = localtime(&ttNow) ;
	char *sDT = asctime(pTime) ; // from manual : The string result produced by asctime contains exactly 26 characters and has the form Wed Jan 02 02:03:55 1980\n\0.
	memcpy(strDT, sDT, 26) ;
	strDT[24] = 0 ;
}

static void GetCurrentDTmsec(char *strDT, int64_t  & tNow)
{
	// get current date/time string so that we can log when stuff is happening
	if (0 == tNow) tNow = ARE::GetTimeInMilliseconds() ;
	sprintf(strDT, "%lld", tNow) ;
}

int ARE::VarElimOrderComp::CVOcontext::NoteVarOrderComputationCompletion(int w_IDX, ARE::Graph & G)
{
	if (G._OrderLength != _Problem->N()) {
		int32_t error = 1 ;
		return 1 ;
		}
	++_nRunsCompleted ;

	int32_t bt_height = -1 ; // BT height of order of G
	int32_t nBTpartitionings = -1 ;
	int32_t ib = -1 ;
	int32_t bt_maxDepthPartitionedNode = -1 ;
	BucketElimination::MBEworkspace *bews = NULL ;
	// check if we need BT
	if (HasBTObjective()) {
		if (NULL == bews) {
			bews = new BucketElimination::MBEworkspace ;
			if (NULL == bews)
				goto done ;
			}
		int res_BT = ComputeTreeDecomposition(*_Problem, G._VarElimOrder, G._VarElimOrderWidth, *bews) ;
		if (0 != res_BT)
			goto done ;
		bt_height = bews->MaxTreeHeight() ;
		// check if has BT partitioning objective
		if (HasBTPartitioningObjective()) {
			double spaceused = -1.0 ;
			int32_t iBoundMin = 2 ;
			bews->MaxSpaceAllowed_Log10() = 9.0 ; // 1GB
			int res_MBE = bews->FindIBoundForSpaceAllowed(iBoundMin, ib, spaceused, nBTpartitionings, bt_maxDepthPartitionedNode) ;
			if (0 != res_MBE) {
				nBTpartitionings = -2 ;
				ib = -2 ;
				bt_maxDepthPartitionedNode = -2 ;
				}
			}
		}

	// check if the new order is better than the old : primary objective
	{
	switch (_ObjCodePrimary) {
		case Width :
			if (G._VarElimOrderWidth > _BestOrder->_Width)
				goto done ;
			if (G._VarElimOrderWidth < _BestOrder->_Width)
				goto is_better ;
			break ;
		case StateSpaceSize :
			if (G._TotalVarElimComplexity_Log10 > _BestOrder->_TotalVarElimComplexity_Log10)
				goto done ;
			if (G._TotalVarElimComplexity_Log10 < _BestOrder->_TotalVarElimComplexity_Log10)
				goto is_better ;
			break ;
		case BTheight :
			if (_BestOrder->_BTheight < 0)
				goto is_better ;
			if (bt_height > _BestOrder->_BTheight)
				goto done ;
			if (bt_height < _BestOrder->_BTheight)
				goto is_better ;
			break ;
		case BTib :
			if (_BestOrder->_BTib < 0)
				goto is_better ;
			if (ib < _BestOrder->_BTib)
				goto done ;
			if (ib > _BestOrder->_BTib)
				goto is_better ;
			break ;
		case nBTnodesPartitioned :
			if (_BestOrder->_nBTpartitionings < 0)
				goto is_better ;
			if (nBTpartitionings > _BestOrder->_nBTpartitionings)
				goto done ;
			if (nBTpartitionings < _BestOrder->_nBTpartitionings)
				goto is_better ;
			break ;
		case BTmaxDepthPartitionedNode :
			if (_BestOrder->_BTmaxDepthPartitionedNode < 0)
				goto is_better ;
			if (bt_maxDepthPartitionedNode > _BestOrder->_BTmaxDepthPartitionedNode)
				goto done ;
			if (bt_maxDepthPartitionedNode < _BestOrder->_BTmaxDepthPartitionedNode)
				goto is_better ;
			break ;
		default :
			goto done ;
		}
	}

	// check if the new order is better than the old : secondary objective
	{
	switch (_ObjCodeSecondary) {
		case Width :
			if (G._VarElimOrderWidth > _BestOrder->_Width)
				goto done ;
			if (G._VarElimOrderWidth < _BestOrder->_Width)
				goto is_better ;
			break ;
		case StateSpaceSize :
			if (G._TotalVarElimComplexity_Log10 > _BestOrder->_TotalVarElimComplexity_Log10)
				goto done ;
			if (G._TotalVarElimComplexity_Log10 < _BestOrder->_TotalVarElimComplexity_Log10)
				goto is_better ;
			break ;
		case BTheight :
			if (_BestOrder->_BTheight < 0)
				goto is_better ;
			if (bt_height > _BestOrder->_BTheight)
				goto done ;
			if (bt_height < _BestOrder->_BTheight)
				goto is_better ;
			break ;
		case BTib :
			if (_BestOrder->_BTib < 0)
				goto is_better ;
			if (ib < _BestOrder->_BTib)
				goto done ;
			if (ib > _BestOrder->_BTib)
				goto is_better ;
			break ;
		case nBTnodesPartitioned :
			if (_BestOrder->_nBTpartitionings < 0)
				goto is_better ;
			if (nBTpartitionings > _BestOrder->_nBTpartitionings)
				goto done ;
			if (nBTpartitionings < _BestOrder->_nBTpartitionings)
				goto is_better ;
			break ;
		case BTmaxDepthPartitionedNode :
			if (_BestOrder->_BTmaxDepthPartitionedNode < 0)
				goto is_better ;
			if (bt_maxDepthPartitionedNode > _BestOrder->_BTmaxDepthPartitionedNode)
				goto done ;
			if (bt_maxDepthPartitionedNode < _BestOrder->_BTmaxDepthPartitionedNode)
				goto is_better ;
			break ;
		default :
			goto done ;
		}
	}

	// two orders are equal; keep old.
	goto done ;

is_better :

	// for info, compute BT stuff always
	if (NULL == bews && _ComputeBTstats) {
		bews = new BucketElimination::MBEworkspace ;
		if (NULL == bews)
			goto done ;
		int res_BT = ComputeTreeDecomposition(*_Problem, G._VarElimOrder, G._VarElimOrderWidth, *bews) ;
		if (0 != res_BT)
			goto done ;
		bt_height = bews->MaxTreeHeight() ;
		}
	if (ib < 0 && NULL != bews) {
		double spaceused = -1.0 ;
		int32_t iBoundMin = 2 ;
		bews->MaxSpaceAllowed_Log10() = 9.0 ; // 1GB
		int res_MBE = bews->FindIBoundForSpaceAllowed(iBoundMin, ib, spaceused, nBTpartitionings, bt_maxDepthPartitionedNode) ;
		if (0 != res_MBE) {
			nBTpartitionings = -2 ;
			ib = -2 ;
			bt_maxDepthPartitionedNode = -2 ;
			}
		}

	int64_t tNow ; tNow = ARE::GetTimeInMilliseconds() ;
	if (NULL != _fpLOG) {
		fprintf(_fpLOG, "\n%I64d worker %2d found better solution : width=%d complexity=%g bt_height=%d space(#elements)=%g", tNow, (int) w_IDX, (int) G._VarElimOrderWidth, (double) G._TotalVarElimComplexity_Log10, (int) bt_height, (double) G._TotalNewFunctionStorageAsNumOfElements_Log10) ;
		fflush(_fpLOG) ;
		}
#ifdef VERBOSE_CVO
	printf("%I64d worker %2d found better solution : width=%d complexity=%g BTh=%d BTib=%d nBTp=%d maxBTPd=%d space(#elements)=%g\n", tNow, (int) w_IDX, (int) G._VarElimOrderWidth, (double) G._TotalVarElimComplexity_Log10, (int) bt_height, (int) ib, (int) nBTpartitionings, (int) bt_maxDepthPartitionedNode, (double) G._TotalNewFunctionStorageAsNumOfElements_Log10) ;
#endif

	{
	ARE::VarElimOrderComp::ResultSnapShot & result_record = _Improvements[_nImprovements++] ;
	result_record._dt = tNow - _tStart ;
	result_record._width = G._VarElimOrderWidth ;
	result_record._complexity = G._TotalVarElimComplexity_Log10 ;
	result_record._BTheight = bt_height ;
	result_record._nBTpartitionings = nBTpartitionings ;
	result_record._BTib = ib ;
	result_record._BTmaxDepthPartitionedNode = bt_maxDepthPartitionedNode ;

	_BestOrder->_Width = G._VarElimOrderWidth ;
	_BestOrder->_nFillEdges = G._nFillEdges ;
	_BestOrder->_MaxSingleVarElimComplexity = G._MaxVarElimComplexity_Log10 ;
	_BestOrder->_BTheight = bt_height ;
	_BestOrder->_nBTpartitionings = nBTpartitionings ;
	_BestOrder->_BTib = ib ;
	_BestOrder->_BTmaxDepthPartitionedNode = bt_maxDepthPartitionedNode ;
	_BestOrder->_TotalVarElimComplexity_Log10 = G._TotalVarElimComplexity_Log10 ;
	_BestOrder->_TotalNewFunctionStorageAsNumOfElements_Log10 = G._TotalNewFunctionStorageAsNumOfElements_Log10 ;
	for (int i = 0 ; i < _Problem->N() ; i++)
		_BestOrder->_VarListInElimOrder[i] = (G._VarElimOrder)[i] ;

#ifdef VERBOSE_CVO
	cout << "c status " << (1+_BestOrder->_Width) << ' ' << tNow << std::endl ;
	cout << flush ;
#endif
	}

done :

	if (NULL != bews)
		delete bews ;

	if (G._VarElimOrderWidth >= 0 && G._VarElimOrderWidth < 1024) {
		_Width2CountMap[G._VarElimOrderWidth]++ ;
		if (_Width2MaxComplexityMap[G._VarElimOrderWidth] > G._TotalVarElimComplexity_Log10)
			_Width2MaxComplexityMap[G._VarElimOrderWidth] = G._TotalVarElimComplexity_Log10 ;
		if (_Width2MaxComplexityMap[G._VarElimOrderWidth] < G._TotalVarElimComplexity_Log10)
			_Width2MaxComplexityMap[G._VarElimOrderWidth] = G._TotalVarElimComplexity_Log10 ;
		}

	return 0 ;
}


#if defined WINDOWS || _WINDOWS
typedef unsigned int (__stdcall *pWorkerThreadFn)(void *X) ;
static unsigned int __stdcall WorkerThreadFn(void *X)
#elif defined (LINUX)
typedef void *(*pWorkerThreadFn)(void *X) ;
static void *WorkerThreadFn(void *X)
#endif
{
	ARE::VarElimOrderComp::Worker *w = (ARE::VarElimOrderComp::Worker *) X;
	ARE::VarElimOrderComp::CVOcontext & CVOcontext = *(w->_CVOcontext) ;
	ARE::ARP & Problem = *(CVOcontext._Problem) ;
	ARE::VarElimOrderComp::Order & best_order = *(CVOcontext._BestOrder) ;

	// initialize best order to something bad; we don't want anything worse than that.
	int bestWidth = 255 ;
	double bestComplexity = DBL_MAX ;
	{
	ARE::utils::AutoLock lock(CVOcontext._BestOrderMutex) ;
	bestWidth = best_order._Width ;
	bestComplexity = best_order._TotalVarElimComplexity_Log10 ;
	}

	int64_t tNow = 0 ;
	int res ;

	if (NULL != CVOcontext._fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d start; objCode1=%d objCode2=%d alg=%d ...", tNow, (int) w->_IDX, (int) CVOcontext._ObjCodePrimary, (int) CVOcontext._ObjCodeSecondary, (int) CVOcontext._AlgCode) ;
		fflush(CVOcontext._fpLOG) ;
		}

	try {
		*(w->_G) = CVOcontext._MasterGraph ;
		if (! w->_G->_IsValid)
			goto done ;
		}
	catch (...) {
		if (NULL != CVOcontext._fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d prep0 exception ...", tNow, (int) w->_IDX) ;
			fflush(CVOcontext._fpLOG) ;
			}
		}

	int nCompleteRunsTodo ; nCompleteRunsTodo = 3 ;
	while (! w->_ThreadStop) {
		tNow = 0 ;
		++(w->_nRunsDone) ;
#if defined WINDOWS || _WINDOWS
		long v = InterlockedIncrement(&(CVOcontext._nRunsStarted)) ;
#else
		pthread_mutex_lock(&nRunsSumMutex) ;
		long v = ++CVOcontext._nRunsStarted ;
		pthread_mutex_unlock(&nRunsSumMutex) ;
#endif
		if (0 == (v % CVOcontext._LogIncrement)) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d will do run %d ...", tNow, (int) w->_IDX, (int) v) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}
// DEBUGGG
//printf("\nworker %d starting; nRunsSum=%d ...", (int) w->_IDX, (int) v) ;
		try {
			bool earlyTerminationOk = nCompleteRunsTodo-- > 0 ? false : true ;
			int32_t widthLimit = bestWidth ;
			double spaceLimit = bestComplexity ;
			if (CVOcontext._FindPracticalVariableOrder) {
				if (widthLimit > CVOcontext._PracticalOrderLimit_W)
					widthLimit = CVOcontext._PracticalOrderLimit_W ;
				if (spaceLimit > CVOcontext._PracticalOrderLimit_C)
					spaceLimit = CVOcontext._PracticalOrderLimit_C ;
				}
			else {
				widthLimit = CVOcontext.HasNonWidthObjective() ? 10 + bestWidth : bestWidth ; // set width cutoff at slightly above best known, so that orderings close to the best (in terms of width) get a chance too
				spaceLimit = DBL_MAX ;
				}
// 2014-03-21 KK : even if MinFill algorithm is used, run Simple(), not Simple_wMinFillOnly(), because Simple() will compute complexity also, as so we can try to minimize complexity too.
//			if (ARE::VarElimOrderComp::MinFill == CVOcontext._AlgCode)
//				res = w->_G->ComputeVariableEliminationOrder_Simple_wMinFillOnly(widthLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_W, false, 1, CVOcontext._nRandomPick, CVOcontext._eRandomPick, w->_TempAdjVarSpace, TempAdjVarSpaceSize) ;
//			else
				res = w->_G->ComputeVariableEliminationOrder_Simple(CVOcontext._AlgCode, widthLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_W, spaceLimit, earlyTerminationOk && CVOcontext._EarlyTerminationOfBasic_C, false, 1, CVOcontext._nRandomPick, CVOcontext._eRandomPick, w->_TempAdjVarSpaceSizeExtraArrayN, w->_TempAdjVarSpaceSizeExtraArray) ;
// DEBUGGG_AAA
//ARE::VarElimOrderComp::DeleteNewAdjVarList(w->_TempAdjVarSpaceSizeExtraArrayN, w->_TempAdjVarSpaceSizeExtraArray) ;
			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %2d exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				goto done ;
				}
			}
// DEBUGGG
//		printf("\nworker %d finished res=%d; width=%d complexity=%I64d", (int) w->_IDX, (int) res, (int) w->_G->_VarElimOrderWidth, (int64_t) w->_G->_TotalVarElimComplexity) ;
//		if (NULL != CVOcontext._fpLOG) {
//			fprintf(CVOcontext._fpLOG, "\nworker %d finished res=%d; width=%d complexity=%I64d", (int) w->_IDX, (int) res, (int) w->_G->_VarElimOrderWidth, (int64_t) w->_G->_TotalVarElimComplexity) ;
//			fflush(CVOcontext._fpLOG) ;
//			}
		if (w->_ThreadStop)
			goto done ; // if stop requested, abandon
		try {
			ARE::utils::AutoLock lock(CVOcontext._BestOrderMutex) ;
			if (w->_ThreadStop)
				goto done ; // if stop requested, abandon
//GetCurrentDTmsec(strDT, tNow) ;
//printf("\n%s worker %2d found width=%d complexity=%I64d space(#elements)=%I64d res=%d", strDT, (int) w->_IDX, (int) w->_G->_VarElimOrderWidth, (int64_t) w->_G->_TotalVarElimComplexity, (int64_t) w->_G->_TotalNewFunctionStorageAsNumOfElements, (int) res) ;
			if (0 == res) {
				CVOcontext.NoteVarOrderComputationCompletion(w->_IDX, *(w->_G)) ;
				}
			else {
				// DEBUGGG
//				printf("\nThread %d res=%d", w->_IDX, res) ;
				}
			bestWidth = best_order._Width ;
			bestComplexity = best_order._TotalVarElimComplexity_Log10 ;
			if (CVOcontext._nRunsStarted >= CVOcontext._nRunsToDoMax) {
// DEBUGGG
//				printf("\nworker %d out of runs (%d >= %d) ...", (int) w->_IDX, (int) nRunsSum, (int) CVOcontext._nRunsToDoMax) ;
				goto done ;
				}
			if (CVOcontext._tToStop > 0) {
				if (0 == tNow) tNow = ARE::GetTimeInMilliseconds() ;
				if (tNow >= CVOcontext._tToStop) {
// DEBUGGG
//					printf("\nworker %d out of time ...", (int) w->_IDX) ;
					goto done ;
					}
				}

			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %d summary exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}

		if (w->_ThreadStop)
			goto done ; // if stop requested, abandon
		try {
			*(w->_G) = CVOcontext._MasterGraph ;
			}
		catch (...) {
			if (NULL != CVOcontext._fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(CVOcontext._fpLOG, "\n%I64d worker %d prep exception ...", tNow, (int) w->_IDX) ;
				fflush(CVOcontext._fpLOG) ;
				}
			}
		}

done :
	if (NULL != CVOcontext._fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(CVOcontext._fpLOG, "\n%I64d worker %d done ...", tNow, (int) w->_IDX) ;
		fflush(CVOcontext._fpLOG) ;
		}
	w->_ThreadHandle = 0 ;
	w->_ThreadIsRunning = false ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}


#if defined WINDOWS || _WINDOWS
typedef unsigned int (__stdcall *pCVOThreadFn)(void *X) ;
static unsigned int __stdcall CVOThreadFn(void *X)
#elif defined (LINUX)
typedef void *(*pCVOThreadFn)(void *X) ;
static void *CVOThreadFn(void *X)
#endif
{
	ARE::VarElimOrderComp::CVOcontext *context = (ARE::VarElimOrderComp::CVOcontext *)(X) ;
	context->Reset() ;
	int nWorkers = context->_nThreads ;
	int nRunsToDoMax = context->_nRunsToDoMax ;
	ARE::ARP & p = *(context->_Problem) ;
	ARE::VarElimOrderComp::Order & BestOrder = *(context->_BestOrder) ;
	int & ret = context->_ret, i ;
	ARE::Graph & OriginalGraph = context->_OriginalGraph ;
	ARE::Graph & MasterGraph = context->_MasterGraph ;
	ARE::VarElimOrderComp::Order & best_order = *(context->_BestOrder) ;

	ARE::VarElimOrderComp::Worker *Workers = NULL ;

	int nRunning, nWrunning ;
	long stop_signalled = 0 ;

	char strDT[64] ;
	int64_t tNow = 0 ; // ARE::GetTimeInMilliseconds() ;
	GetCurrentDTmsec(strDT, tNow) ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%s CVO control thread; start preprocessing ...", strDT) ;
		fflush(context->_fpLOG) ;
		}

	if (NULL == BestOrder._VarListInElimOrder) {
		if (0 != BestOrder.Initialize(p.N())) {
			ret = 1000 ;
			if (NULL != context->_fpLOG) {
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; bestorder init failed ...", tNow) ;
				fflush(context->_fpLOG) ;
				}
			goto done ;
			}
		}

	// create problem graph
	if (context->_RandomGeneratorSeed > 0)
		OriginalGraph.RNG().seed(context->_RandomGeneratorSeed) ; // set seed so that starting point can be duplicated
	OriginalGraph.Create(p) ;
	if (! OriginalGraph._IsValid) {
		ret = 1001 ;
		goto done ;
		}

	// do all the easy eliminations; this will give us a starting point for large-scale randomized searches later.
	tNow = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; eliminate easy vars from graph ...", tNow) ;
		fflush(context->_fpLOG) ;
		}
	MasterGraph = OriginalGraph ;
	if (context->_RandomGeneratorSeed > 0)
		MasterGraph.RNG().seed(context->_RandomGeneratorSeed) ; // set seed so that starting point can be duplicated
	if (! MasterGraph._IsValid) {
		ret = 1002 ;
		goto done ;
		}
// 2014-03-21 KK : even if MinFill algorithm is used, run Simple(), not Simple_wMinFillOnly(), because Simple() will compute complexity also, as so we can try to minimize complexity too.
//	if (ARE::VarElimOrderComp::MinFill == context->_AlgCode)
//		i = MasterGraph.ComputeVariableEliminationOrder_Simple_wMinFillOnly(INT_MAX, false, false, 1, 1, 0.0, context->_TempAdjVarSpace, TempAdjVarSpaceSize) ;
//	else
		i = MasterGraph.ComputeVariableEliminationOrder_Simple(0, INT_MAX, false, DBL_MAX, false, true, 1, 1, 0.0, context->_TempAdjVarSpaceSizeExtraArrayN, context->_TempAdjVarSpaceSizeExtraArray) ;
	// check if the problem was solved completely
	if (MasterGraph._OrderLength >= MasterGraph._nNodes) {
		context->NoteVarOrderComputationCompletion(-1, MasterGraph) ;
		ret = 0 ;
		goto done ;
		}
// DEBUGGG_AAA
//ARE::VarElimOrderComp::DeleteNewAdjVarList(context->_TempAdjVarSpaceSizeExtraArrayN, context->_TempAdjVarSpaceSizeExtraArray) ;
	MasterGraph.ReAllocateEdges() ;
	tNow = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; %d vars eliminated, %d remaining ...", tNow, (int) MasterGraph._OrderLength, (int) MasterGraph._nRemainingNodes) ;
		fflush(context->_fpLOG) ;
		}

	context->_tStart = tNow ;
	if (context->_TimeLimitInMilliSeconds > 0)
		context->_tToStop = context->_tStart + context->_TimeLimitInMilliSeconds ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; preprocessing done; ready to run : nRunsToDo=%d/%d, tLimitmsec=%I64d (tAutoStop=%I64d) ...", tNow, context->_nRunsToDoMin, context->_nRunsToDoMax, context->_TimeLimitInMilliSeconds, context->_tToStop) ;
		fflush(context->_fpLOG) ;
		}

	{
	ARE::Graph *g = new ARE::Graph ;
	if (NULL == g)
		{ ret = 1003 ; goto done ; }
	*g = MasterGraph ;
	if (g->_IsValid) {
		if (context->_RandomGeneratorSeed > 0)
			g->RNG().seed(1+context->_RandomGeneratorSeed) ; // set seed so that starting point can be duplicated; move by 1.
		i = g->ComputeVariableEliminationOrder_LowerBound() ;
		best_order._WidthLowerBound = g->_VarElimOrderWidth ;
		if (NULL != context->_fpLOG) {
			fprintf(context->_fpLOG, "\n%s lower_bound=%d", strDT, (int) best_order._WidthLowerBound) ;
			fflush(context->_fpLOG) ;
			}
		}
	delete g ;
	}

	// if strictly best order (whatever the width/complexity) is required, execute one run here to get some real bound on width/complexity.
	// when we launch multi-threaded search for good order, knowing a decent bound will help the threads right away.
	// when we want a practical order, width/complexity limits are (should be) quite low.
	best_order._Width = p.N() ;
	best_order._TotalVarElimComplexity_Log10 = DBL_MAX ;
	// 2015-12-14 KK : do this always, so that we have some result (order); we may quite quickly, and if we did  not run this, we may have nothing.
//	if (! context->_FindPracticalVariableOrder) {
		{
		ARE::Graph g ;
		g = MasterGraph ;
		if (g._IsValid) {
			if (context->_RandomGeneratorSeed > 0)
				g.RNG().seed(2+context->_RandomGeneratorSeed) ; // set seed so that starting point can be duplicated; move by 2.
			int widthLimit = INT_MAX ;
			double spaceLimit = DBL_MAX ;
			i = g.ComputeVariableEliminationOrder_Simple(0, widthLimit, false, spaceLimit, false, false, 10, 1, 0.0, context->_TempAdjVarSpaceSizeExtraArrayN, context->_TempAdjVarSpaceSizeExtraArray) ;
			++context->_nRunsStarted ;
			tNow = 0 ;
			GetCurrentDTmsec(strDT, tNow) ;
			if (0 == i) {
				if (NULL != context->_fpLOG) {
					fprintf(context->_fpLOG, "\n%s Initial computation width=%d; MaxSingleVarElimComplexity=%g, TotalVarElimComplexity=%g, TotalNewFunctionStorageAsNumOfElements=%g", strDT, (int) g._VarElimOrderWidth, (double) g._MaxVarElimComplexity_Log10, (double) g._TotalVarElimComplexity_Log10, (double) g._TotalNewFunctionStorageAsNumOfElements_Log10) ;
					fflush(context->_fpLOG) ;
					}
				context->NoteVarOrderComputationCompletion(-1, g) ;
				}
			else {
				if (NULL != context->_fpLOG) {
					fprintf(context->_fpLOG, "\n%s Initial computation failed; res=%d; will continue ...", strDT, i) ;
					fflush(context->_fpLOG) ;
					}
				}
			}
		}

	if (context->_nRunsStarted >= context->_nRunsToDoMax)
		goto done ;
#if defined WINDOWS || _WINDOWS
	stop_signalled = InterlockedCompareExchange(&(context->_StopAndExit), 1, 1) ;
#else
	pthread_mutex_lock(&stopSignalMutex);
	stop_signalled = context->_StopAndExit;
	pthread_mutex_unlock(&stopSignalMutex);
#endif
	if (0 != stop_signalled) {
		if (NULL != context->_fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(context->_fpLOG, "\n%I64d CVO control thread; ready to create workers, but stop signalled ...", tNow) ;
			fflush(context->_fpLOG) ;
			}
		goto done ;
		}

	if (NULL != context->_fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; nWorkers=%d", tNow, nWorkers) ;
		fflush(context->_fpLOG) ;
		}

	if (nWorkers < 1) {
		if (NULL != context->_fpLOG) {
			fprintf(context->_fpLOG, "\nNo threads; will quit ...") ;
			fflush(context->_fpLOG) ;
			goto done ;
			}
		}
//goto done ;

	// create workers; don't let them run yet
	Workers = new ARE::VarElimOrderComp::Worker[nWorkers] ;
	if (NULL == Workers) {
		ret = 1004 ;
		if (NULL != context->_fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(context->_fpLOG, "\n%I64d CVO control thread; failed to create workers, will quit ...", tNow) ;
			fflush(context->_fpLOG) ;
			}
		goto done ;
		}
	for (i = 0 ; i < nWorkers ; i++) {
		Workers[i]._IDX = i ;
		Workers[i]._G = new ARE::Graph ;
		Workers[i]._CVOcontext = context ;
		if (context->_RandomGeneratorSeed > 0)
			Workers[i]._G->RNG().seed(100+i+context->_RandomGeneratorSeed) ; // set random seed; move by 100+idx.
		if (NULL == Workers[i]._G) {
			ret = 1005 ;
			goto done ;
			}
		}

	context->_LogIncrement = nRunsToDoMax/20 ;
	if (context->_LogIncrement < 1) context->_LogIncrement = 1 ;
	nWrunning = 0 ;
	for (i = 0 ; i < nWorkers ; i++) {
		Workers[i]._ThreadIsRunning = true ;
		Workers[i]._ThreadStop = false ;
#if defined WINDOWS || _WINDOWS
		Workers[i]._ThreadHandle = _beginthreadex(NULL, 0, WorkerThreadFn, &(Workers[i]), 0, NULL) ;
#else
		pthread_create(&(Workers[i]._ThreadHandle), NULL, WorkerThreadFn, &(Workers[i])) ; // TODO third argument
#endif
		if (0 == Workers[i]._ThreadHandle) {
			Workers[i]._ThreadIsRunning = false ;
			Workers[i]._ThreadStop = true ;
			if (NULL != context->_fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; FAILED to create thread %d", tNow, i) ;
				fflush(context->_fpLOG) ;
				}
			}
		else
			nWrunning++ ;
		}
	if (0 == nWrunning) {
		if (NULL != context->_fpLOG) {
			tNow = ARE::GetTimeInMilliseconds() ;
			fprintf(context->_fpLOG, "\n%I64d CVO control thread; no worker threads, will quite ...", tNow) ;
			fflush(context->_fpLOG) ;
			}
		ret = 1006 ;
		goto done ;
		}

	// loop here waiting for stop signal of until all workers threads have exited
	nRunning = -1 ;
	stop_signalled = 0 ;
	while (true) {
		SLEEP(50) ;
#if defined WINDOWS || _WINDOWS
		stop_signalled = InterlockedCompareExchange(&(context->_StopAndExit), 1, 1) ;
#else
		pthread_mutex_lock(&stopSignalMutex);
		stop_signalled = context->_StopAndExit;
		pthread_mutex_unlock(&stopSignalMutex);
#endif
		if (0 != stop_signalled) {
			if (NULL != context->_fpLOG) {
				tNow = ARE::GetTimeInMilliseconds() ;
				fprintf(context->_fpLOG, "\n%I64d CVO control thread; stop_signalled ...", tNow) ;
				fflush(context->_fpLOG) ;
				}
			}
		nRunning = 0 ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			if (0 != stop_signalled) {
				if (NULL != context->_fpLOG) {
					tNow = ARE::GetTimeInMilliseconds() ;
					fprintf(context->_fpLOG, "\n%I64d CVO control thread; stopping worker %d ...", tNow, i) ;
					fflush(context->_fpLOG) ;
					}
				Workers[i]._ThreadStop = true ;
				}
			++nRunning ;
			}
		if (nRunning <= 0 || 0 != stop_signalled)
			break ;
		}
	// wait for threads to stop
	while (nRunning > 0) {
		SLEEP(50) ;
		nRunning = 0 ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			++nRunning ;
			}
		}
	if (NULL != context->_fpLOG) {
		tNow = ARE::GetTimeInMilliseconds() ;
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; all worker threads have closed ...", tNow) ;
		fflush(context->_fpLOG) ;
		}

/*
	if (TimeLimitInSeconds > 0) {
		Sleep(1000*TimeLimitInSeconds) ;
		for (i = 0 ; i < nWorkers ; i++)
			Workers[i]._ThreadStop = true ;
		for (i = 0 ; i < nWorkers ; i++) {
			if (0 == Workers[i]._ThreadHandle) continue ;
			DWORD r = WaitForSingleObject((HANDLE) Workers[i]._ThreadHandle, 10) ;
			if (WAIT_TIMEOUT == r)
				TerminateThread((HANDLE) Workers[i]._ThreadHandle, 0) ;
			CloseHandle((HANDLE) Workers[i]._ThreadHandle) ;
			Workers[i]._ThreadHandle = 0 ;
			}
		}
*/

	ret = 0 ;
done :
	tNow = context->_tEnd = ARE::GetTimeInMilliseconds() ;
	if (NULL != context->_fpLOG) {
		fprintf(context->_fpLOG, "\n%I64d CVO control thread; done; time=%I64dmsec ...", tNow, (int64_t) (context->_tEnd - context->_tStart)) ;
		fflush(context->_fpLOG) ;
		}
//	for (i = 0 ; i < nWorkers ; i++)
//		nRuns += Workers[i]._nRunsDone ;
	if (NULL != Workers)
		delete [] Workers ;
	if (best_order._Width < p.N()) {
		// some ordering was found
		}

	// remove redundant fill edges
	if (false && best_order._Width > 0) {
/* 2015-12-01 KK : commented out this block of code.
		ARE::Graph g ;
		g = OriginalGraph ;
		g._OrderLength = p.N() ;
		for (i = 0 ; i < p.N() ; i++) {
			int v_i = best_order._VarListInElimOrder[i] ;
			g._VarType[i] = 0 ;
			g._PosOfVarInList[v_i] = i ;
			g._VarElimOrder[i] = v_i ;
			}
		g._nFillEdges = best_order._nFillEdges ;
		g.RemoveRedundantFillEdges() ;
		int nEdgesRemoved = (OriginalGraph._nEdges + best_order._nFillEdges) - g._nEdges ;
		fprintf(context->_fpLOG, "\n# of redundant fill edges = %d", nEdgesRemoved) ;
		fflush(context->_fpLOG) ;
*/
		}

	context->_ThreadHandle = 0 ;
#if defined WINDOWS || _WINDOWS
	_endthreadex(0) ;
	return 0  ;
#else
	return NULL ;
#endif
}


int ARE::VarElimOrderComp::CVOcontext::CreateCVOthread(void)
{
#if defined WINDOWS || _WINDOWS
	_ThreadHandle = _beginthreadex(NULL, 0, CVOThreadFn, this, 0, NULL) ;
#else
	pthread_create(&_ThreadHandle, NULL, CVOThreadFn, this) ; // TODO third argument
#endif
	return 0 != _ThreadHandle ? 0 : 1 ;
}


int ARE::VarElimOrderComp::CVOcontext::RequestStopCVOthread(void)
{
	if (0 == _ThreadHandle) {
		if (NULL != _fpLOG) {
			int64_t tNowLog = ARE::GetTimeInMilliseconds() ;
			fprintf(_fpLOG, "\n%I64d CVO_th : request stop variable order computation; already stopped ...", tNowLog) ;
			fflush(_fpLOG) ;
			}
		return 0 ;
		}
#if defined WINDOWS || _WINDOWS
	long previous_value = InterlockedCompareExchange(&_StopAndExit, 1, 0) ;
	if (previous_value > 0)
		return -1 ;
#else
	long previous_value = 0 ;
	pthread_mutex_lock(&stopSignalMutex) ;
	if (0 == _StopAndExit) {
		_StopAndExit = 1 ;
		}
	else
		previous_value = 1 ;
	pthread_mutex_unlock(&stopSignalMutex) ;
	if (previous_value > 0)
		return -1 ;
#endif
	return 1 ;
}


int ARE::VarElimOrderComp::CVOcontext::StopCVOthread(int64_t TimeoutInMilliseconds)
{
	if (0 == _ThreadHandle) {
		if (NULL != _fpLOG) {
			int64_t tNowLog = ARE::GetTimeInMilliseconds() ;
			fprintf(_fpLOG, "\n%I64d CVO_th : stop variable order computation; already stopped ...", tNowLog) ;
			fflush(_fpLOG) ;
			}
		return 0 ;
		}
#if defined WINDOWS || _WINDOWS
	InterlockedCompareExchange(&_StopAndExit, 1, 0) ;
#else
	pthread_mutex_lock(&stopSignalMutex);
	if (_StopAndExit == 0) {
		_StopAndExit = 1 ;
		}
	pthread_mutex_unlock(&stopSignalMutex);
#endif
	int64_t tStart = ARE::GetTimeInMilliseconds() ;
	if (NULL != _fpLOG) {
		fprintf(_fpLOG, "\n%I64d    CVO_th : stop variable order computation; stop signalled, will wait ...", tStart) ;
		fflush(_fpLOG) ;
		}
	if (TimeoutInMilliseconds < 1)
		TimeoutInMilliseconds = 1 ;
	else if (TimeoutInMilliseconds > 86400000)
		TimeoutInMilliseconds = 86400000 ;
	while (true) {
		SLEEP(50) ;
		if (0 == _ThreadHandle)
			break ;
		int64_t tNow = ARE::GetTimeInMilliseconds() ;
		int64_t dt = tNow - tStart ;
		if (dt > TimeoutInMilliseconds) {
			// we asked the thread to stop and waited for it to stop, but it won't stop, so kill the thread.
			if (NULL != _fpLOG) {
				int64_t tNowLog = ARE::GetTimeInMilliseconds() ;
				fprintf(_fpLOG, "\n%I64d CVO_th : stop variable order computation, hard kill ...", tNowLog) ;
				fflush(_fpLOG) ;
				}
#if defined WINDOWS || _WINDOWS
			TerminateThread((HANDLE) _ThreadHandle, 0) ;
			CloseHandle((HANDLE) _ThreadHandle) ;
			_ThreadHandle = 0 ;
#else
			pthread_exit(&_ThreadHandle);
#endif
			break ;
			}
		}
	return 0 ;
}


int ARE::VarElimOrderComp::Compute(
	// IN
	const std::string & ProblemInputFile,
	ARE::VarElimOrderComp::ObjectiveToMinimize objCodePrimary,
	ARE::VarElimOrderComp::ObjectiveToMinimize objCodeSecondary,
	ARE::VarElimOrderComp::NextVarPickCriteria algcode,
	int nthreads,
	int nrunstodo,
	int64_t TimeLimitInMilliSeconds,
	int nRP,
	double eRP,
	bool PerformSingletonConsistencyChecking,
	bool EliminateSingletonDomainVariables,
	bool earlyterminationofbasic_W,
	bool earlyterminationofbasic_C,
	bool FindPracticalVariableOrder,
	unsigned long random_seed,
	// OUT
	ARE::VarElimOrderComp::Order & BestOrder,
	ARE::VarElimOrderComp::CVOcontext * & Context
	)
{
	if (ProblemInputFile.length() < 1)
		return 1 ;

	char filetype ;

	int ret = 1 ;

	char strDT[64] ;
	int64_t tNow = 0 ;

	int64_t tStart;
	int64_t tStopSignalled;
	int i;

	ARE::VarElimOrderComp::CVOcontext *cvocontext = Context ;
	ARE::VarElimOrderComp::CVOcontext *localCVOcontext = NULL ;
	if (NULL == cvocontext) {
		localCVOcontext = new ARE::VarElimOrderComp::CVOcontext ;
		if (NULL == localCVOcontext)
			goto done ;
		Context = cvocontext = localCVOcontext ;
		}
	cvocontext->_FindPracticalVariableOrder = FindPracticalVariableOrder ;
	GetCurrentDTmsec(strDT, tNow) ;
#ifdef VERBOSE_CVO
	printf("\n%s CVO : Start ...", strDT) ;
#endif
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s CVO : Start ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	if (NULL == cvocontext->_Problem) {
		cvocontext->_Problem = new ARE::ARP("DefaultCVOproblem") ;
		if (NULL == cvocontext->_Problem)
			goto done ;
		}
	ARE::ARP *p ; p = cvocontext->_Problem ;
	cvocontext->_AlgCode = algcode ;
	cvocontext->_ObjCodePrimary = objCodePrimary ;
	cvocontext->_ObjCodeSecondary = objCodeSecondary ;
	cvocontext->_RandomGeneratorSeed = random_seed ;

	cvocontext->_nRunsToDoMax = nrunstodo ;
	if (cvocontext->_nRunsToDoMax < 1)
		cvocontext->_nRunsToDoMax = 1 ;
	else if (cvocontext->_nRunsToDoMax > 1000000000)
		cvocontext->_nRunsToDoMax = 1000000000 ;
	cvocontext->_nRunsToDoMin = cvocontext->_nRunsToDoMax ;

	if (TimeLimitInMilliSeconds < 0)
		TimeLimitInMilliSeconds = 0 ;
	else if (TimeLimitInMilliSeconds > 86400000)
		TimeLimitInMilliSeconds = 86400000 ;
	cvocontext->_TimeLimitInMilliSeconds = TimeLimitInMilliSeconds ;

	cvocontext->_eRandomPick = eRP ;
	cvocontext->_nRandomPick = nRP ;
	if (cvocontext->_nRandomPick < 1)
		cvocontext->_nRandomPick = 1 ;

	cvocontext->_EarlyTerminationOfBasic_W = earlyterminationofbasic_W ;
	cvocontext->_EarlyTerminationOfBasic_C = earlyterminationofbasic_C ;
	cvocontext->_BestOrder = &BestOrder ;

	cvocontext->_nThreads = nthreads ;
	if (cvocontext->_nThreads > cvocontext->_nRunsToDoMax-1)
		cvocontext->_nThreads = cvocontext->_nRunsToDoMax-1 ;
	else if (cvocontext->_nThreads < 1)
		cvocontext->_nThreads = 1 ;
// DEBUGGG
//nThreads = 1 ;

	// at least one of nRunsToDo/TimeLimitInSeconds has to be given
	if (cvocontext->_nRunsToDoMin < 1 && cvocontext->_TimeLimitInMilliSeconds < 1)
		cvocontext->_nRunsToDoMin = 1 ;

	// fn may have dir in it; extract filename.
	i = ProblemInputFile.length() - 1 ;
	for (; i >= 0 ; i--) {
#ifdef LINUX
		if ('\\' == ProblemInputFile[i] || '/' == ProblemInputFile[i]) break ;
#else
		if ('\\' == ProblemInputFile[i] || '//' == ProblemInputFile[i]) break ;
#endif
		}
	filetype = 0 ; // 1=uai, 2=gr
	{
	std::string fn(ProblemInputFile.substr(i+1)) ;
	p->SetName(fn) ;
	}

	if (0 != p->LoadFromFile(ProblemInputFile)) {
		ret = 2 ;
#ifdef VERBOSE_CVO
		printf("\nload failed ...") ;
#endif
		goto done ;
		}
	if (0 != p->PerformPostConstructionAnalysis()) {
		ret = 3 ;
#ifdef VERBOSE_CVO
		printf("\nPerformPostConstructionAnalysis failed ...") ;
#endif
		goto done ;
		}
	if (p->N() < 1) {
		ret = 4 ;
#ifdef VERBOSE_CVO
		printf("\nN=%d; will exit ...", p->N()) ;
#endif
		goto done ;
		}
/*
	if (p.N() > (1 << 16)) {
//		printf("\nN=%d; currently nVars limit is 64K since AVL tree (edgesadded) with key type long is used to store 2 variable indeces; will exit ...", p.N()) ;
		printf("\nN=%d; currently nVars limit is 64K since heap uses 2 bytes for cost and 2 bytes for var index, for total 4 bytes; change heap key to __in64; will exit ...", p.N()) ;
		goto done ;
		}
*/

	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
#ifdef VERBOSE_CVO
	printf("\n%s File loaded; start preprocessing; N=%d ...", strDT, p->N()) ;
#endif
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s File loaded; start preprocessing ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	// eliminate singleton-domain variables; do this before ordering is computed; this is easy and should be done by any algorithm processing the network
	if (EliminateSingletonDomainVariables) {
#ifdef VERBOSE_CVO
		printf("\n%s EliminateSingletonDomainVariables(): start, N=%d ...", strDT, (int) p->N()) ;
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s EliminateSingletonDomainVariables(): start, N=%d ...", strDT, (int) p->N()) ;
			fflush(cvocontext->_fpLOG) ;
			}
		p->EliminateSingletonDomainVariables() ;
		tNow = 0 ;
		GetCurrentDTmsec(strDT, tNow) ;
#ifdef VERBOSE_CVO
		printf("\n%s EliminateSingletonDomainVariables(): done; nSingletonDomainVariables=%d ...", strDT, (int) p->nSingletonDomainVariables()) ;
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s EliminateSingletonDomainVariables(): done; nSingletonDomainVariables=%d ...", strDT, (int) p->nSingletonDomainVariables()) ;
			fflush(cvocontext->_fpLOG) ;
			}
		}
	else {
#ifdef VERBOSE_CVO
		printf("\n%s EliminateSingletonDomainVariables not requested ...", strDT);
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s EliminateSingletonDomainVariables not requested ...", strDT);
			fflush(cvocontext->_fpLOG);
			}
		}

	// prune domains of variables by checking which values participate in no complete assignment with probability>0.
	// basically, we compute a minimal domain for each variable.
	if (PerformSingletonConsistencyChecking) {
		int32_t nNewSingletonDomainVariablesAfterSingletonConsistency = 0, nVarsWithReducedDomainSize = 0 ;
		int32_t res_sc = p->ComputeSingletonConsistency(nNewSingletonDomainVariablesAfterSingletonConsistency, nVarsWithReducedDomainSize) ;
		tNow = 0 ;
		GetCurrentDTmsec(strDT, tNow) ;
		if (-1 == res_sc) {
			// domain of some variable is empty
#ifdef VERBOSE_CVO
			printf("\n%s ComputeSingletonConsistency(): domain of some variable is empty; will quit ...", strDT) ;
#endif
			if (NULL != cvocontext->_fpLOG) {
				fprintf(cvocontext->_fpLOG, "\n%s ComputeSingletonConsistency(): domain of some variable is empty; will quit ...", strDT) ;
				fflush(cvocontext->_fpLOG) ;
				}
			goto done ;
			}
#ifdef VERBOSE_CVO
		printf("\n%s ComputeSingletonConsistency(): nNewSingletonDomainVariablesAfterSingletonConsistency = %d", strDT, nNewSingletonDomainVariablesAfterSingletonConsistency) ;
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s ComputeSingletonConsistency(): nNewSingletonDomainVariablesAfterSingletonConsistency = %d", strDT, nNewSingletonDomainVariablesAfterSingletonConsistency) ;
			fflush(cvocontext->_fpLOG) ;
			}
		if (nNewSingletonDomainVariablesAfterSingletonConsistency > 0)
			p->EliminateSingletonDomainVariables() ;
		}
	else {
#ifdef VERBOSE_CVO
		printf("\n%s Singleton-Consistency check not requested ...", strDT) ;
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\n%s Singleton-Consistency check not requested ...", strDT) ;
			fflush(cvocontext->_fpLOG) ;
			}
		}

	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
#ifdef VERBOSE_CVO
	printf("\n%s Launching CVO thread ...\n", strDT) ;
#endif
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s Launching CVO thread ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	cvocontext->CreateCVOthread() ;
	if (0 == cvocontext->_ThreadHandle) {
#ifdef VERBOSE_CVO
		printf("\nFAILED to create cvo thread ...") ;
#endif
		if (NULL != cvocontext->_fpLOG) {
			fprintf(cvocontext->_fpLOG, "\nFAILED to create cvo thread ...") ;
			fflush(cvocontext->_fpLOG) ;
			}
		ret = 8 ;
		goto done ;
		}

	tStart = ARE::GetTimeInMilliseconds();
	tStopSignalled = 0 ;
	while (true) {
		SLEEP(100) ;
		if (0 == cvocontext->_ThreadHandle)
			break ;
		int64_t tNow = ARE::GetTimeInMilliseconds() ;
		int64_t dt = tNow - tStart ;
		if (dt < cvocontext->_TimeLimitInMilliSeconds)
			continue ;
		if (0 == tStopSignalled) {
			tStopSignalled = tNow ;
#if defined WINDOWS || _WINDOWS
			InterlockedCompareExchange(&(cvocontext->_StopAndExit), 1, 0) ;
#else
			pthread_mutex_lock(&stopSignalMutex);
			if (cvocontext->_StopAndExit == 0) {
				cvocontext->_StopAndExit = 1;
				}
			pthread_mutex_unlock(&stopSignalMutex);
#endif
			continue ;
			}
		dt = tNow - tStopSignalled ;
		if (dt > 10000) {
			// we asked the thread to stop and waited for it to stop, but it won't stop, so kill the thread.
#if defined WINDOWS || _WINDOWS
			TerminateThread((HANDLE) cvocontext->_ThreadHandle, 0) ;
			CloseHandle((HANDLE) cvocontext->_ThreadHandle) ;
			cvocontext->_ThreadHandle = 0 ;
#else
			// TODO
			pthread_exit(&cvocontext->_ThreadHandle);
#endif
			break ;
			}
		}
	tNow = 0 ;
	GetCurrentDTmsec(strDT, tNow) ;
#ifdef VERBOSE_CVO
	printf("\n%s CVO : thread has closed ...", strDT) ;
#endif
	if (NULL != cvocontext->_fpLOG) {
		fprintf(cvocontext->_fpLOG, "\n%s CVO : thread has closed ...", strDT) ;
		fflush(cvocontext->_fpLOG) ;
		}

	ret = 0 ;
done :

	return ret ;
}


int32_t ARE::VarElimOrderComp::ComputeTreeDecomposition(ARE::ARP & P, int32_t *Order, int32_t Width, BucketElimination::MBEworkspace & bews)
{
	if (Width < 0 || Width >= INT_MAX) {
		return -1 ;
		}

	int res_copyorder = P.SetVarElimOrdering(Order, Width) ;

	int res_BEwsInit = bews.Initialize(P, true, NULL, 0) ;
	if (0 != res_BEwsInit)
		return 1 ;
	int res_BEwsCB = bews.CreateBuckets(true, true, true, false) ;
	if (0 != res_BEwsCB)
		return 1 ;
	for (int i = 0 ; i < bews.nBuckets() ; i++) {
		BucketElimination::Bucket *b = bews.getBucket(i) ;
		if (NULL == b) continue ;
		if (b->Width() < 0)
			b->ComputeSignature() ;
		}
	if (bews.MaxNumVarsInBucket() < 0)
		bews.ComputeMaxNumVarsInBucket(false) ;
	if (bews.MaxTreeHeight() < 0)
		bews.ComputeBucketTreeStatistics() ;

	return 0 ;
}


int32_t ARE::VarElimOrderComp::Order::SerializeTreeDecomposition(ARE::ARP & P, BucketElimination::MBEworkspace & bews, bool one_based_indexing, bool ConnectedComponents, std::string & sOutput)
{
	char s[128] ;
	if (P.N() <= 0) {
		sOutput = "s td 0 0 0" ;
		return 0 ;
		}
	if (_Width < 0 || _Width >= INT_MAX) {
		sprintf(s, "s td 0 0 %d", (int) P.N()) ; // s td <nBags> <tw+1> <N>
		sOutput = s ;
		return 0 ;
		}

	int res_BT = ComputeTreeDecomposition(P, _VarListInElimOrder, _Width, bews) ;
	if (0 != res_BT)
		return res_BT ;

	int extra_idx = one_based_indexing? 1 : 0 ;

/*
// DEBUGGG
	int idxv4 = bews.VarPos()[3] ;
	BucketElimination::Bucket *b4 = bews.MapVar2Bucket(3) ;
	if (idxv4 < 0 || idxv4 >= bews.N())
		{ sprintf(s, "c ERORR v4 has no ordeer idx\n") ; sOutput += s ; }
	if (NULL == b4)
		{ sprintf(s, "c ERORR v4 has no bucket\n") ; sOutput += s ; }
	int idxBv4 = b4->IDX() ; { sprintf(s, "c INFO v4 bidx=%d (bV=%d) idxv4inOrder=%d\n", extra_idx + idxBv4, extra_idx + b4->V(), extra_idx + idxv4) ; sOutput += s ; }
// DEBUGGG
	BucketElimination::Bucket *b7070 = bews.getBucket(7070) ;
	int vB7070 = b7070->V() ;
		{ sprintf(s, "c INFO b7071 v=%d width=%d nC=%d P=%c\n", extra_idx + vB7070, b7070->Width(), b7070->nChildren(), NULL != b7070->ParentBucket() ? 'y': 'n') ; sOutput += s ; }
// DEBUGGG
	int c4 = 0 ;
	for (int i = 0 ; i < bews.N() ; i++) {
		if (3 == bews.VarOrder()[i]) c4++ ;
		}
		{ sprintf(s, "c INFO c4=%d\n", c4) ; sOutput += s ; }
*/
	// if there are connected components, we may need to connect them to dummy root
	int32_t idxDummyRoot = -1 ;
	if (ConnectedComponents && bews.Problem()->nConnectedComponents() > 1)
		idxDummyRoot = extra_idx + bews.nBuckets() ;

	sprintf(s, "s td %d %d %d", (int) (bews.nBuckets() + (idxDummyRoot >= 0 ? 1 : 0)), bews.MaxNumVarsInBucket(), (int) P.N()) ; // s td <nBags> <tw+1> <N>
	sOutput += s ;

	// nBags X vars of the bag
	for (int i = 0 ; i < bews.nBuckets() ; i++) {
		BucketElimination::Bucket *b = bews.getBucket(i) ;
		if (NULL == b) continue ;
		sprintf(s, "\nb %d", extra_idx + b->IDX()) ; sOutput += s ;
		for (int j = 0 ; j < b->Width() ; j++)
			{ sprintf(s, " %d", extra_idx + b->Signature()[j]) ; sOutput += s ; }
		}
	// if there are connected components, we may need to connect them to dummy root
	std::vector<int> roots ;
	if (ConnectedComponents && idxDummyRoot >= 0) {
		idxDummyRoot = extra_idx + bews.nBuckets() ;
		sprintf(s, "\nb %d", idxDummyRoot) ; sOutput += s ;
		roots.reserve(bews.Problem()->nConnectedComponents()) ;
		for (int i = 0 ; i < bews.nBuckets() ; i++) {
			BucketElimination::Bucket *b = bews.getBucket(i) ;
			if (NULL == b) continue ;
			if (NULL == b->ParentBucket())
				roots.push_back(extra_idx + b->IDX()) ;
			}
		}

	// edges of the bucket tree
	for (int i = 0 ; i < bews.nBuckets() ; i++) {
		BucketElimination::Bucket *b = bews.getBucket(i) ;
		if (NULL == b) continue ;
		BucketElimination::Bucket *p = b->ParentBucket() ;
		if (NULL == p) continue ;
		sprintf(s, "\n%d %d", extra_idx + b->IDX(), extra_idx + p->IDX()) ; sOutput += s ;
		}
	// extra edges to connect BT roots to dummy root
	if (ConnectedComponents && idxDummyRoot >= 0) {
		for (int i = 0 ; i < roots.size() ; i++)
			{ sprintf(s, "\n%d %d", idxDummyRoot, roots[i]) ; sOutput += s ; }
		}

	return 0 ;
}

#if defined DEFINE_PACE16_MAIN_FN

static ARE::VarElimOrderComp::Order BestOrder ;
static ARE::VarElimOrderComp::CVOcontext Context ;
static int64_t tStart = 0 ;
static int nThreads2Use = -1 ;
static long nTDprintsAttempted = 0 ;
static long nTDprintsDone = 0 ;

static int SerializeBestOrderTD(std::string *fn)
{
	// count how many times it has been printed
	long v = -1 ;
	{
#if defined WINDOWS || _WINDOWS
		v = InterlockedIncrement(&nTDprintsAttempted) ;
#else
		pthread_mutex_lock(&printTDMutex) ;
		v = ++nTDprintsAttempted ;
		pthread_mutex_unlock(&printTDMutex) ;
#endif
	}

	// just one print
	if (v <= 1) {
		FILE *fp = NULL != fn ? fopen(fn->c_str(), "w") : NULL ;
		ARE::utils::AutoLock lock(Context._BestOrderMutex) ;
//		int64_t tNow = ARE::GetTimeInMilliseconds() ;
//		cout << "c BEST ORDER so far : N=" << Context._Problem->N() << " width=" << BestOrder._Width << " lowerbound=" << BestOrder._WidthLowerBound << " varElimComplexity(log10)=" << BestOrder._TotalVarElimComplexity_Log10 << " nRunsStarted=" <<  Context._nRunsStarted << " nRunsCompleted=" << Context._nRunsCompleted << " nImprovements=" << Context._nImprovements << " nTrivialVars=" << Context._MasterGraph._OrderLength << " runtime=" << (tNow - tStart) << "msec" << " nThreads=" << nThreads2Use << std::endl ;
		std::string sTD ;
		BucketElimination::MBEworkspace bews ;
		BestOrder.SerializeTreeDecomposition(*(Context._Problem), bews, true, true, sTD) ;
		if (NULL != fp) {
			fwrite(sTD.c_str(), 1, sTD.length(), fp) ;
			fflush(fp) ;
			}
		else {
//			cout << "c mbews stats : N=" << bews.N() << " nCC=" << bews.Problem()->nConnectedComponents() << " nB=" << bews.nBuckets() << " nVwoB=" << bews.nVarsWithoutBucket() << " maxDTheight=" << bews.MaxTreeHeight() << " nBwoC=" << bews.nBucketsWithNoChildren() << " nRoots=" << bews.nRoots() << std::endl ;
			cout << sTD ;
			cout << flush ;
			}
		if (NULL != fp)
			fclose(fp) ;
		++nTDprintsDone ;
		}

	return 0 ;
}

static int SerializeBestOrderW(void)
{
	ARE::utils::AutoLock lock(Context._BestOrderMutex) ;
	cout << (1 + BestOrder._Width) << "\n" ;
}

#ifdef LINUX
static void handle_signal(int signal)
{
    sigset_t pending;

    // Find out which signal we're handling
    switch (signal) {
        case SIGUSR1:
            SerializeBestOrderW() ;
            break;
        case SIGINT:
        case SIGTERM:
			Context.RequestStopCVOthread() ;
            SerializeBestOrderTD(NULL) ;
			// wait until print is finished before quitting
			while (nTDprintsDone <= 0) {
				SLEEP(1) ;
				}
            exit(0);
        default:
            return;
    }
}
#endif

int main(int argc, char* argv[])
{
	int nParams = argc ;

	unsigned long randomGeneratorSeed = 999999937 ;
	std::string problem_filename ;
	int nrunstodo = 1000000000 ;
	int nArgs = (nParams-1)>>1 ;
	bool findPracticalVariableOrder = true ;
	ARE::VarElimOrderComp::ObjectiveToMinimize objCodePrimary = ARE::VarElimOrderComp::Width, objCodeSecondary = ARE::VarElimOrderComp::None ;
	if (1 + 2*nArgs != nParams) {
		printf("\nBAD COMMAND LINE; will exit ...") ;
		return 1 ;
		}
	for (int i = 0 ; i < nArgs; i++) {
		std::string sArgID = (NULL != argv[1+2*i] ? argv[1+2*i] : "") ;
		std::string sArg   = (NULL != argv[2*(i+1)] ? argv[2*(i+1)] : "") ;
		if (0 == stricmp("-s", sArgID.c_str()))
			randomGeneratorSeed = std::strtoul(sArg.c_str(), NULL, 0) ;
		else if (0 == stricmp("-f", sArgID.c_str()))
			problem_filename = sArg ;
		else if (0 == stricmp("-nR", sArgID.c_str()))
			nrunstodo = atoi(sArg.c_str());
		else if (0 == stricmp("-fpvo", sArgID.c_str()))
			findPracticalVariableOrder = '1' == sArg[0] || 'y' == sArg[0] || 'Y' == sArg[0] ;
		else if (0 == stricmp("-O1", sArgID.c_str()))
			objCodePrimary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
		else if (0 == stricmp("-O2", sArgID.c_str()))
			objCodeSecondary = (ARE::VarElimOrderComp::ObjectiveToMinimize) atoi(sArg.c_str()) ;
		}
	if (nrunstodo < 1)
		nrunstodo = 1 ;

	int maxNumProcessorThreads = std::thread::hardware_concurrency() ; // this requires c++11
	if (maxNumProcessorThreads <= 0)
		maxNumProcessorThreads = 1 ;

	// if filename is given, figure out file type
	bool file_is_uai = false ;
	std::string problem_filename_wo_dir ;
	std::vector<std::string> files2process ;
	unsigned int inBufSize = 0 ;
	if (problem_filename.length() > 0 ? '\\' == problem_filename[problem_filename.length()-1] || '/' == problem_filename[problem_filename.length()-1] : false) {
		int res = ReadDirectory(problem_filename.c_str(), NULL, NULL, "uai", files2process) ;
		for (int32_t ii = 0 ; ii < files2process.size() ; ++ii) {
			std::string & fn = files2process[ii] ;
			fn = problem_filename + fn ;
			}
		problem_filename.clear() ;
		}
	else if (problem_filename.length() > 0) {
		files2process.push_back(problem_filename) ;
/*		int i = problem_filename.length() - 1;
		for (; i >= 0; --i) {
#ifdef LINUX
			if ('\\' == problem_filename[i] || '/' == problem_filename[i])
#else
			if ('\\' == problem_filename[i] || '//' == problem_filename[i])
#endif
				break;
			}
		problem_filename_wo_dir = problem_filename.substr(i + 1);
		std::size_t pos = problem_filename_wo_dir.find(".uai") ;
		if (std::string::npos != pos)
			file_is_uai = true ;*/
		}
	else { // if (0 == problem_filename.length()) {
		problem_filename = "cin.gr" ;
		files2process.push_back(problem_filename) ;
		}

#ifdef LINUX
	// set up signal handling
	struct sigaction sa ;
	sa.sa_handler = handle_signal ;
//	sigemptyset(&sa.sa_mask) ;
	sigfillset(&sa.sa_mask) ; // block every signal during handling
//	sa.sa_flags = 0 ;
	sa.sa_flags = SA_RESTART ; // Restart functions if interrupted by handler
	sigaction(SIGTERM, &sa, NULL) ;
	sigaction(SIGUSR1, &sa, NULL) ;
#endif

#ifdef RUN_FOLDER_PROBLEMS
	cout << "RUN_FOLDER_PROBLEMS\n" ;
	std::wstring wfolder(L"D:\\UCI\\PACE-testbed\\problems\\") ;
	std::string folder = wstringToString(wfolder) ;

	std::wstring wfiles[1024] ; std::string files[1024] ; std::string output[1024] ; int nFiles = 0 ;
	{
		HANDLE hFile ;
		WIN32_FIND_DATA FileInformation ;
		std::wstring strPattern(wfolder) ;
		strPattern += L"*" ;
//		char fn[260] ;
		hFile = ::FindFirstFile(strPattern.c_str(), &FileInformation) ;
		if (INVALID_HANDLE_VALUE != hFile) {
			do {
				bool isdot = '.' == FileInformation.cFileName[0] ;
				bool isdir = 0 != (FileInformation.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) ;
				if (isdir) {
					}
				else if (isdot) {
					}
				else {
					wfiles[nFiles] = FileInformation.cFileName ;
					if (wfiles[nFiles].length() < 4) continue ;
					std::wstring wEXT(wfiles[nFiles].substr(wfiles[nFiles].length()-3, 3)) ;
					if (L".gr" != wEXT) continue ;
					files[nFiles] = folder ;
					files[nFiles] += wstringToString(wfiles[nFiles]) ;
					output[nFiles] = files[nFiles] ;
					output[nFiles].erase(output[nFiles].length()-3) ;
					output[nFiles] += "-out.td" ;
					nFiles++ ;
					}
				} while(TRUE == ::FindNextFile(hFile, &FileInformation)) ;
			::FindClose(hFile) ;
			}
	}
	cout << "Processing " << nFiles << " files ..." << "\n" ;

	for (int ii = 0 ; ii < nFiles ; ii++) {
		std::string & fn = files[ii] ;
		std::string & of = output[ii] ;
		nTDprintsAttempted = nTDprintsDone = 0 ;
		problem_filename = fn ;
		Context.Destroy() ;
		cout << "Processing " << fn << "\n" ;
#else
	std::string sResultsSummaryWidthList ;
	for (int32_t ii = 0 ; ii < files2process.size() ; ++ii) {
		std::string & fn = files2process[ii] ;
		std::size_t pos = fn.find(".uai") ;
		if (std::string::npos != pos)
			file_is_uai = true ;
#endif

		int32_t lenFN = fn.length() - 1, jj = lenFN-1 ;
		for (; jj >= 0 ; --jj) {
#ifdef LINUX
			if ('\\' == fn[jj] || '/' == fn[jj])
#else
			if ('\\' == fn[jj] || '//' == fn[jj])
#endif
				break ;
			}
		problem_filename_wo_dir = jj < 0 ? fn : fn.substr(jj + 1) ;

		int64_t TimeLimitInMilliSeconds = 86400000 ;
		int nRP = 8 ; double eRP = 0.5 ;
		bool PerformSingletonConsistencyChecking = false, EliminateSingletonDomainVariables = file_is_uai ? true : false, earlyterminationofbasic_W = true, earlyterminationofbasic_C = false ;
		nThreads2Use = maxNumProcessorThreads > 1 ? maxNumProcessorThreads - 1 : 1 ;
#ifdef VERBOSE_CVO
		printf("\nnThreads2Use=%d nrunstodo=%d TimeLimitInMilliSeconds=%lld O1=%d O2=%d fpvo=%c", (int)nThreads2Use, (int)nrunstodo, (int64_t)TimeLimitInMilliSeconds, (int) objCodePrimary, (int) objCodeSecondary, findPracticalVariableOrder ? 'Y' : 'N');
#endif
		BestOrder.Destroy() ;
		Context.Destroy() ;
		Context._BestOrder = &BestOrder ;
		ARE::VarElimOrderComp::CVOcontext *context = &Context ;
		tStart = ARE::GetTimeInMilliseconds();
		int res = ARE::VarElimOrderComp::Compute(fn,
			objCodePrimary, objCodeSecondary, ARE::VarElimOrderComp::MinFill,
			nThreads2Use, nrunstodo, TimeLimitInMilliSeconds, nRP, eRP,
			PerformSingletonConsistencyChecking, EliminateSingletonDomainVariables, earlyterminationofbasic_W, earlyterminationofbasic_C, findPracticalVariableOrder, randomGeneratorSeed,
			BestOrder, context) ;
		int64_t tEnd = ARE::GetTimeInMilliseconds();

#ifdef VERBOSE_CVO
		printf("\nBEST ORDER width = %d, varElimComplexity = %f, lower bound = %d, nRunsDone = %d/%d, nImprovements = %d, nTrivialVars = %d, runtime = %lldmsec\n", (int) BestOrder._Width, (double)BestOrder._TotalVarElimComplexity_Log10, (int) BestOrder._WidthLowerBound, (int)Context._nRunsStarted, (int)Context._nRunsCompleted, (int) Context._nImprovements, (int) Context._MasterGraph._OrderLength, (int64_t) (tEnd - tStart)) ;
#endif

		// add to summary list
		{
		if (sResultsSummaryWidthList.length() > 0) sResultsSummaryWidthList += '\n' ;
		char s[64] ;
		sprintf(s, "width=%d ", (int) BestOrder._Width) ;
		sResultsSummaryWidthList += s ; sResultsSummaryWidthList += problem_filename_wo_dir ;
		}

		// save order
#ifdef SAVE_ORDER
		std::string vofn(fn);
		std::string::size_type ext_pos = vofn.rfind('.');
		vofn.erase(ext_pos);
		vofn += "_var_elim_order";
		if (objCodePrimary > ARE::VarElimOrderComp::None) {
			vofn += "-" ;
			vofn += ARE::VarElimOrderComp::MapObjectiveToMinimize2Str(objCodePrimary) ;
			}
		if (objCodeSecondary > ARE::VarElimOrderComp::None) {
			vofn += "-" ;
			vofn += ARE::VarElimOrderComp::MapObjectiveToMinimize2Str(objCodeSecondary) ;
			}
		vofn += ".vo";
		BestOrder.SerializeAsElimOrder(vofn.c_str());
#endif

		// save best tree decomposition
#ifdef RUN_FOLDER_PROBLEMS
		SerializeBestOrderTD(&of) ;
#elif PRINT_BEST_TD_AT_END
		SerializeBestOrderTD(NULL) ;
#else
		nTDprintsDone = 999 ;
#endif
//	std::string sTD ;
//	BestOrder.SerializeTreeDecomposition(*(Context._Problem), true, sTD) ;
//	cout << "\n\n" << sTD ;

//	if (NULL != Context)
//		{ delete Context ; Context = NULL ; }

		// wait this printout is done; SIGTERM may have caused it too.
		while (nTDprintsDone <= 0) {
			SLEEP(1) ;
			}
#ifdef RUN_FOLDER_PROBLEMS
		}
	FILE *fp_validate = fopen("validate.bat", "w") ;
	for (int ii = 0 ; ii < nFiles ; ii++) {
		std::string & fn = files[ii] ;
		std::string & of = output[ii] ;
		fprintf(fp_validate, "td-validate %s %s\n", fn.c_str(), of.c_str()) ;
		}
	fflush(fp_validate) ; fclose(fp_validate) ;
#else
		}
#endif

	if (sResultsSummaryWidthList.length() > 0) {
		FILE *fpSummary = fopen("summary.txt", "w") ;
		if (NULL != fpSummary) {
			fwrite(sResultsSummaryWidthList.c_str(), 1, sResultsSummaryWidthList.length(), fpSummary) ;
			fclose(fpSummary) ;
			}
		}

	return 0 ;
}
#endif
