#include <stdlib.h>
#include <memory.h>

#include <Function.hxx>
#include <Function-NN.hxx>
#include <Bucket.hxx>
#include <MBEworkspace.hxx>
#include <Bucket.hxx>
#include <MiniBucket.hxx>
#include <Sort.hxx>
#include "DATA_SAMPLES.h"
#include "Utils/MersenneTwister.h"
#include <exception>
#include "Config.h"
#include <chrono>

//#include <iostream> //delete
//#include <fstream> //delete
static MTRand RNG ;

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int32_t BucketElimination::MiniBucket::ComputeOutputFunction_NN(int32_t varElimOperator, ARE::Function *FU, ARE::Function *fU, double WMBEweight)
{
    auto start = std::chrono::high_resolution_clock::now();
//    int32_t i, j, k ;
     int32_t i, k;
    bool convert_exp=false;
    if (0 != global_config.network.compare("net")){
        printf("Masked net ----");
        convert_exp = true;
    }
    else
        printf("Net ...");
    // TODO : if _OutputFunction is of type FunctionNN, then do .....
    ARE::FunctionNN *fNN = dynamic_cast<ARE::FunctionNN *>(_OutputFunction) ;
    if (NULL == fNN)
        return 1 ;

    MBEworkspace *bews = dynamic_cast<MBEworkspace*>(_Workspace) ;
    if (NULL == bews)
        return ERRORCODE_generic ;
    ARE::ARP *problem = bews->Problem() ;
    if (NULL == problem)
        return ERRORCODE_generic ;
    if (_Width < 0) {
        if (0 != ComputeSignature())
            return ERRORCODE_generic ;
        if (_Width < 0)
            return ERRORCODE_generic ;
    }
    if (nVars() < 1)
        // nothing to do; should not happen; however, allow this to pass as ok; calling fn should be able to handle this.
        return 0 ;

    const int32_t w = Width() ;
    if (w < 0 || _nFunctions < 1)
        return 0 ;
    const int32_t *signature = Signature() ;

    // generate some number of random samples...
    Config config;
    int32_t nSamples = global_config.sample_size;

    std::vector<int32_t> values_, vars_ ;
    values_.resize(problem->N(), 0) ;
    vars_.resize(problem->N(), 0) ;

    if (values_.size() != problem->N() || vars_.size() != problem->N())
        return 1 ;
    int32_t *vals = values_.data();
    int32_t *vars = vars_.data();

    k = problem->K(_V);
    for ( i = 0 ; i < fNN->N() ; ++i)
        vars[i] = fNN->Argument(i) ;

    printf("\n In Minibucket NN --");
    vars[i] = _V ;

   // int32_t ** samples_signiture;
   // samples_signiture = new int32_t *[nSamples];

    double s = 0.8;
    double v_s = 0.1;
    int n_val_samples = int(v_s*nSamples);
    int n_train_samples = int(s*nSamples);
    int n_test_samples = int(v_s*nSamples);

    int32_t ** train_samples;
    train_samples = new int32_t *[n_train_samples];

    int32_t ** val_samples;
    val_samples = new int32_t *[n_val_samples];

    int32_t ** test_samples;
    test_samples = new int32_t *[n_test_samples];

    for (int nn=0; nn<int(s*nSamples); nn++){
       // samples_signiture[nn] = new int32_t[fNN->N()];
        train_samples[nn] = new int32_t[fNN->N()];

    }
    for (int nn=0; nn<int(v_s*nSamples); nn++){
        // samples_signiture[nn] = new int32_t[fNN->N()];
        test_samples[nn] = new int32_t[fNN->N()];
    }
    for (int nn=0; nn<int(v_s*nSamples); nn++){
        // samples_signiture[nn] = new int32_t[fNN->N()];
        val_samples[nn] = new int32_t[fNN->N()];
    }
   // float* sample_values = new float[nSamples];
    float* train_sample_values = new float[n_train_samples];
    float* val_sample_values = new float[n_val_samples];
    float* test_sample_values = new float[n_test_samples];

   //printf("Declared train,val,test arrays. Size : %d,%d,%d", n_train_samples,n_val_samples,n_test_samples);
    ARE_Function_TableType const_factor = bews->FnCombinationNeutralValue() ;
    int32_t nFNs = 0 ;
    std::vector<ARE::Function *> flist ; flist.reserve(_nFunctions) ; if (flist.capacity() != _nFunctions) return 1 ;
    for (int32_t j = 0 ; j < _nFunctions ; j++) {
        ARE::Function *f = _Functions[j] ;
        if (NULL == f) return ERRORCODE_fn_ptr_is_NULL ;
        if (0 == f->N()) bews->ApplyFnCombinationOperator(const_factor, f->ConstValue()) ;
        else {
            flist.push_back(f) ;
            f->ComputeArgumentsPermutationList(w, vars); }
    }
    float zero = 0.1*pow(10,-34);
    int count_non_zero =0, t_non_zero=0, v_non_zero=0, test_non_zero=0;
    int c1=0;
    printf("Uniform Sampling -----");

    //if save samples for a bucket
    //ofstream myfile ("bucket_samples_level_" + std::to_string(this->Workspace()->count_Train) + ".txt",std::ios_base::app);

    for (i = 0; i < nSamples; ++i) {
        // generate assignment to fNN arguments
        ARE_Function_TableType V = bews->VarEliminationDefaultValue();
        for (int32_t j = 0; j < fNN->N(); ++j) {
            int32_t v = fNN->Argument(j);
            int32_t domain_size_of_v = problem->K(v);
            //            printf("%d ", domain_size_of_v);
            int32_t value = RNG.randInt(domain_size_of_v - 1);
            vals[j] = value;
        }
        // enumerate all current variable values; compute bucket value for each configuration and combine them using elimination operator...

        for (int32_t j = 0; j < k; ++j) {
            vals[fNN->N()] = j;
            ARE_Function_TableType v = bews->FnCombinationNeutralValue();
            // compute value for this configuration : fNN argument assignment + _V=j
            for (int32_t l = 0; l < _nFunctions; ++l) {
                ARE::Function *f = _Functions[l];
                if (NULL == f) continue;
                double fn_v = f->TableEntryEx(vals,
                                              problem->K()); //This would make us problem specifially when privious ones be NN
                bews->ApplyFnCombinationOperator(v, fn_v);
            }
            ApplyVarEliminationOperator(varElimOperator, problem->FunctionsAreConvertedToLogScale(), V, v);
        }

        bews->ApplyFnCombinationOperator(V, const_factor);
        //save the samples in arrays to make them tensors
        c1++;
        //uncomment to save samples of a bucket
        /*
        if (myfile.is_open()) {
            for (int32_t m = 0; m < fNN->N(); ++m) {
                myfile << vals[m] << "\t";
            }
            if (convert_exp) {
                myfile << exp(V) << std::endl;
            }
            else {
                myfile << V << std::endl;
            }

            //myfile.close();
        }
        */
        if (i < n_train_samples) {
            for (int32_t m = 0; m < fNN->N(); ++m) {
                train_samples[i][m] = vals[m];
            }
            if (convert_exp) {
                V = exp(V);
            }
            if(V>=zero){
                t_non_zero++;
            }
            train_sample_values[i] = V; //should be V
        } else if (i >= n_train_samples && i < n_train_samples + n_val_samples) {
            for (int32_t m = 0; m < fNN->N(); ++m) {
                val_samples[i - n_train_samples][m] = vals[m];
            }
            if (convert_exp) {
                V = exp(V);
            }
            if(V>=zero){
                v_non_zero++;
            }
            val_sample_values[i - n_train_samples] = V; //should be V
        } else if (i >= n_train_samples + n_val_samples && i < n_train_samples + n_val_samples + n_test_samples) {
            for (int32_t m = 0; m < fNN->N(); ++m) {
                test_samples[i - n_train_samples - n_val_samples][m] = vals[m]; //vals[fNN->ArgumentsPermutationList()[m]]; // vals[fnn->_Argument]
            }
            if (convert_exp) {
                V = exp(V);
            }
            if(V>=zero){
                test_non_zero++;
                // printf("sample value - %f, count -- %d",V,count_non_zero);
            }
           // printf("Non-zero sample test ---- %f, zero -- %f",V, zero);
            test_sample_values[i - n_train_samples - n_val_samples] = V; //should be V
        }
        if(V>=zero){
            count_non_zero++;
           // printf("sample value - %f, count -- %d",V,count_non_zero);
        }
        if (V > fNN->max_value) {
            fNN->max_value = V;
        }
    }
    std::string o_file = global_config.out_file + ".txt";
    std::ofstream f;
    f.open(o_file,std::ios_base::app);
    //sampling finished
    if(count_non_zero<0.001*nSamples)
    {
        printf("exiting the program because number of non-zero values is not enough; %d out of %d samples",count_non_zero,nSamples);
        if (f.is_open()) {
            f <<  "program exited" << std::endl;
        }
        f.close();
        exit(0);
    } else
    {
       // if (f.is_open()) {
       //     f << count_non_zero << '\t' << t_non_zero << '\t' << v_non_zero << '\t' << test_non_zero << '\t';
       // }
        f.close();
    }
    DATA_SAMPLES *DS_train, *DS_val, *DS_test;
    t_non_zero=0, v_non_zero=0, test_non_zero=0;
    DS_train = fNN->samples_to_data(train_samples, train_sample_values, fNN->N(), n_train_samples, t_non_zero);
    DS_val = fNN->samples_to_data(val_samples, val_sample_values, fNN->N(), n_val_samples,v_non_zero);
    printf("samples created  ---%d %d", t_non_zero, v_non_zero);
    fNN->Train(DS_train, DS_val);
    DS_test = fNN->samples_to_data(test_samples, test_sample_values, fNN->N(), n_test_samples,test_non_zero);
    fNN->test(DS_test);

    auto elapsed = std::chrono::high_resolution_clock::now() - start;
    long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
    double duration = microseconds*1.0 / 1000000;
    //printf(" Table Duration: %f\n", duration);
    this->Workspace()->time_ComputeOutputFunction_NN += duration;
    this->Workspace()->count_ComputeOutputFunction_NN++;

    //write stats to file
    /*
    f.open(o_file,std::ios_base::app);
    if (f.is_open())
    {
        printf("writing to the file %s", o_file.c_str());
        f << c1 << '\t' << this->Workspace()->time_ComputeOutputFunction_NN << '\t'  << this->Workspace()->time_ComputeOutputFunction_NN/this->Workspace()->count_ComputeOutputFunction_NN << '\t' << this->Workspace()->time_TableEntryEx << '\t' << this->Workspace()->time_TableEntryEx / this->Workspace()->count_TableEntryEx << '\t' << this->Workspace()->time_Train << '\t' << this->Workspace()->time_Train / this->Workspace()->count_Train << std::endl ;
        f << c1 << std::endl ;
        f.close();
    }
    */

    printf(" ComputeOutputFunction_NN Time: %f Count: %d Average: %f\n", this->Workspace()->time_ComputeOutputFunction_NN, this->Workspace()->count_ComputeOutputFunction_NN, this->Workspace()->time_ComputeOutputFunction_NN/this->Workspace()->count_ComputeOutputFunction_NN);
    printf(" TableEntryEx Time: %f Count: %d Average: %f\n", this->Workspace()->time_TableEntryEx, this->Workspace()->count_TableEntryEx, this->Workspace()->time_TableEntryEx / this->Workspace()->count_TableEntryEx);
    printf(" Train Time: %f Count: %d Average: %f\n", this->Workspace()->time_Train, this->Workspace()->count_Train, this->Workspace()->time_Train / this->Workspace()->count_Train);

    printf("Bucket number ---- %d", _V);
    return 0 ;
}




