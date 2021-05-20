#ifndef FunctionNN_HXX_INCLUDED
#define FunctionNN_HXX_INCLUDED

#include <climits>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>

#include "Utils/Mutex.h"
#include "Utils/Sort.hxx"
#include "Problem/Globals.hxx"
#include "Problem/Function.hxx"
#include "Problem/Workspace.hxx"
#include "Net.h"
#include "DATA_SAMPLES.h"
#include <torch/torch.h>
#include "Config.h"
#include <chrono>
#include <fstream>

namespace BucketElimination { class Bucket ; class MiniBucket ; }

namespace ARE
{

class ARP ;

class FunctionNN : public Function
{

private:
   //Net* model = NULL;
   Masked_Net * model = NULL;
//   DATA_SAMPLES * DS = NULL;
//   Config config;
//   torch::Tensor empty_tensor = torch::empty(_nArgs);
//   torch::DeviceType device_type_inf = torch::kCPU;
//   torch::Device device_inf(device_type_inf);

public :
    float max_value = std::numeric_limits<float>::min();
	virtual int32_t AllocateTableData(void)
	{
		// NO TABLE!!!
		return 0 ;
	}

	virtual ARE_Function_TableType TableEntryEx(int32_t *BEPathAssignment, const int32_t *K) const 
	/*
		desc = return fn value corresponding to given input configuration...
		BEPathAssignment = assignment to all variables on the path from the bucket to the root of the bucket tree...
	    this fn is a wrapper for these two lines... when functions are table-based...
		adr = fn->ComputeFnTableAdr_wrtLocalPermutation(BEPathAssignment, K) ;
		double v = fn->TableEntry(adr) ;
		which is 
		return fn->TableEntry(fn->ComputeFnTableAdr_wrtLocalPermutation(BEPathAssignment, K)) ;
		OVERWRITE in a NN-based fn...
	*/
	{
        bool isNet,exp_converted=true;
        isNet = 0 == global_config.network.compare("net");
        if(isNet)
            exp_converted=false;
        auto start = std::chrono::high_resolution_clock::now();
        torch::DeviceType device_type = torch::kCPU;
        torch::Device device(device_type);
        model->to(device);
        model->eval();
        auto empty_tensor = torch::empty(_nArgs);
        float* myData = empty_tensor.data_ptr<float>();

//        printf("******************************\n");
//        printf("BEPathAssignment:\n");
//        for (int i=0; i<_nArgs; i++){
//            printf("%d ", BEPathAssignment[i]);
//        }
//        printf("\n_ArgumentsPermutationList:\n");
//        for (int i=0; i<_nArgs; i++){
//            printf("%d ", _ArgumentsPermutationList[i]);
//        }
////        printf("\nBEPathAssignment[_ArgumentsPermutationList[i]:\n");
//        for (int i=0; i<_nArgs; i++){
//            printf("%d ", BEPathAssignment[_ArgumentsPermutationList[i]]);
//        }
//        printf("\n");

        for (int i=0; i <_nArgs; i++)
            *myData++ = (float)BEPathAssignment[_ArgumentsPermutationList[i]];//(float)BEPathAssignment[_ArgumentsPermutationList[i]]
        torch::Tensor input = empty_tensor.resize_(_nArgs).clone();
        input = input.to(device);
        auto output = model->forward(input,true,isNet);
        double out_value = (double)output.x.item<double>()*max_value;
        // torch::Tensor output = model->forward(input);
        //calculating log_value here
        if(exp_converted)
            out_value = log(out_value);
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        double duration = microseconds*1.0 / 1000000;
        //printf(" Table Duration: %f\n", duration);
        this->WS()->time_TableEntryEx += duration;
        this->WS()->count_TableEntryEx++;
		return out_value;
	}

public :

	virtual void Initialize(Workspace *WS, ARP *Problem, int32_t IDX)
	{

		Function::Initialize(WS, Problem, IDX) ;
		//here we assume the structure is known just initialize the most simple FF neural network with no training.

	}

	void Destroy(void)
	{
		// TODO own stuff here...
		Function::Destroy() ;

	}
	void After_Train(void)
    {
        torch::DeviceType device_type = torch::kCPU;
        torch::Device device(device_type);
        model->to(device);
        model->eval();
    }

  //  torch::Tensor weighted_mse_loss(auto input,auto target,auto weight);
    void confusion_matrix(torch::Tensor prediction,torch::Tensor truth, int &true_positives_batch, int &false_positives_batch, int &true_negatives_batch, int &false_negatives_batch){
            torch::Tensor confusion_vector = prediction / truth;
            true_positives_batch = (confusion_vector==1).sum().template item<int>();
            false_positives_batch = (confusion_vector== float('inf')).sum().template item<int>();
            true_negatives_batch = isnan(confusion_vector).sum().template item<int>();
            false_negatives_batch = (confusion_vector==0).sum().template item<int>();
        }

    void Train(DATA_SAMPLES *DS_train, DATA_SAMPLES *DS_val){
        auto start = std::chrono::high_resolution_clock::now();
        _TableData = NULL;
        _TableSize = 0;
        bool isNet=false;
        if (0 == global_config.network.compare("net")){
            isNet=true;
        }
        if (model == NULL)
            model = new Masked_Net(_nArgs);
        Masked_Net * w_model = new Masked_Net(_nArgs);
        torch::DeviceType device_type;
        if (torch::cuda::is_available()) {
            std::cout << "CUDA available! Training on GPU." << std::endl;
            device_type =  torch::kCUDA;//kCUDA
        } else {
            std::cout << "Training on CPU." << std::endl;
            device_type = torch::kCPU;
        }
        torch::Device device(device_type);
        w_model->to(device);
        w_model->train();
        // Making data  loaders
        auto dataset_train = DS_train->map(torch::data::transforms::Stack<>());
        int64_t batch_size = global_config.batch_size;
        auto data_loader_train = torch::data::make_data_loader<torch::data::samplers::RandomSampler>(dataset_train,batch_size);

        auto dataset_val = DS_val->map(torch::data::transforms::Stack<>());
        auto data_loader_val = torch::data::make_data_loader<torch::data::samplers::RandomSampler>(dataset_val,batch_size);
        //Optimizer for NN
        torch::optim::Adam optimizer(w_model->parameters(), torch::optim::AdamOptions(global_config.lr));

        int64_t n_epochs = global_config.n_epochs;
        float best_mse = std::numeric_limits<float>::max();
        float mse=0, val_mse=0, prev_val_mse=std::numeric_limits<float>::max();
        float w_mse=0, w_val_mse=0,mse_to_compare=0, c_loss=0, val_c_loss=0;
        int count=0, epoch,epoch_t=0;
        bool isTest;
        float zero = 0.1*pow(10,-34);;
        float non_zeros =0.0,train_non_zeros=0.0, val_non_zeros=0.0;
        int total_batch =0, neg_batch =0, neg=0;
        int true_positives=0, false_positives=0, true_negatives=0, false_negatives=0;
        int true_positives_batch=0, false_positives_batch=0, true_negatives_batch=0, false_negatives_batch=0;
        int true_positives_t=0, false_positives_t=0, true_negatives_t=0, false_negatives_t=0;
        int prev_false_negatives=10000,best_false_negatives=10000, fn_to_compare=10000;
        float val_l_t=0. ,  w_val_l_t =0., loss_to_compare=0. ;
        torch::Tensor val_loss_total, val_w_loss_total, w_loss_total, loss_total;

        for (epoch = 1; epoch <= n_epochs; epoch++) {
            isTest=false;
            size_t batch_idx = 0;
            size_t val_batch_idx = 0;
            mse = 0.; // mean squared error
            val_mse = 0.;
            c_loss = 0.;
            val_c_loss = 0.;
            w_mse = 0.; // weighted mean squared error
            w_val_mse = 0.;
            for (auto &batch : *data_loader_train) {
                auto imgs = batch.data;
                torch::Tensor loss, w_loss, loss_n;
                auto labels = batch.target.squeeze();
                imgs = imgs.to(device);
                labels = labels.to(device);
                optimizer.zero_grad();
                auto output = w_model->forward(imgs,isTest,isNet);
                if (isNet){
                    auto w = labels;
                    w = w/w.sum();
                    w_loss = (w*((output.x - labels).pow(2))).mean();
                    loss = torch::nn::functional::mse_loss(labels, output.x);
                    w_loss_total = w_loss;
                    loss_total = loss;
                }
                else{
                    auto w = labels;
                    float w1 = w.sum().template item<float>();
                    //auto output = w_model->forward(imgs, isTest,isNet);
                    if(w1!=0){
                        w = w/w.sum();
                    }
                    float w_max = w.max().template item<float>();
                    auto label_binary = torch::where(labels>0, torch::ones_like(labels), torch::zeros_like(labels));
                    loss_n =  torch::binary_cross_entropy(output.masked,label_binary);
                    non_zeros = (label_binary==1).sum().template item<float>();

                    w_loss_total = global_config.loss_weight_mse*(w*((output.x - labels).pow(2))).mean() + global_config.loss_weight *loss_n;
                    loss_total = global_config.loss_weight_mse*torch::nn::functional::mse_loss(labels, output.x) + global_config.loss_weight *loss_n;
                    loss = torch::nn::functional::mse_loss(labels, output.x);
                    w_loss = (w*((output.x - labels).pow(2))).mean();
                    neg_batch = (label_binary==0).sum().template item<float>();
                    neg = neg + neg_batch;
                    confusion_matrix(output.masked, label_binary, true_positives_batch, false_positives_batch, true_negatives_batch, false_negatives_batch);
                    true_positives_t += true_positives_batch;
                    false_positives_t += false_positives_batch;
                    true_negatives_t += true_negatives_batch;
                    false_negatives_t += false_negatives_batch;
                }
                if(0 == global_config.s_method.compare("is")){
                    w_loss_total.backward();
                }
                else
                    loss_total.backward();
                optimizer.step();
                float l = loss.template item<float>();
                float w_l = w_loss.template item<float>();
                mse += l;
                w_mse += w_l;
                if (!isNet){
                    float c_l = loss_n.template item<float>();
                    c_loss += c_l;
                    train_non_zeros += non_zeros;
                }
                batch_idx++;
            }
            printf("\t  %d %d %d %d \t ",true_positives_t, false_positives_t, true_negatives_t, false_negatives_t);
            mse /= (float) batch_idx;
            w_mse /= (float) batch_idx;
            if (!isNet){
                c_loss /= (float) batch_idx;
                printf("Non zeros in training set -- %d, cross-entropy loss -- %f", train_non_zeros, c_loss);
            }
            printf("Epoch number : %d, train mse : %f ", epoch, mse );
            torch::Tensor val_loss, w_val_loss, val_loss_n;

            // evaluating on validation set
            true_positives=0, false_positives=0, true_negatives=0, false_negatives=0;
            printf("Validation set error to calculate at epoch %d",epoch);
            isTest=true;
            batch_idx = 0;
            for (auto &batch : *data_loader_val) {
                auto imgs = batch.data;
                auto labels = batch.target.squeeze();
                imgs = imgs.to(device);
                labels = labels.to(device);
                auto output = w_model->forward(imgs,isTest,isNet);
                auto w = labels;
                float w1 = w.sum().template item<float>();
                if(w1!=0){
                    w = w/w.sum();
                }
                if (!isNet){
                    auto label_binary = torch::where(labels>0, torch::ones_like(labels), torch::zeros_like(labels));
                    val_loss_n =  torch::binary_cross_entropy(output.masked,label_binary);
                    val_c_loss += val_loss_n.template item<float>();
                    non_zeros = (label_binary==1).sum().template item<float>();
                    val_non_zeros += non_zeros;
                    total_batch = label_binary.sizes()[0];
                    neg_batch = (label_binary==0).sum().template item<float>();
                    neg_batch = total_batch - non_zeros;
                    neg = neg + neg_batch;
                    confusion_matrix(output.masked, label_binary, true_positives_batch, false_positives_batch, true_negatives_batch, false_negatives_batch);
                    val_w_loss_total = global_config.loss_weight_mse*(w*((output.x - labels).pow(2))).mean() + global_config.loss_weight *val_loss_n;
                    val_loss_total = global_config.loss_weight_mse*torch::nn::functional::mse_loss(labels, output.x) + global_config.loss_weight *val_loss_n;
                    val_loss = torch::nn::functional::mse_loss(labels, output.x);
                    true_positives += true_positives_batch;
                    false_positives += false_positives_batch;
                    true_negatives += true_negatives_batch;
                    false_negatives += false_negatives_batch;
                }
                w_val_loss = (w*((output.x - labels).pow(2))).mean();
                val_loss = torch::nn::functional::mse_loss(labels, output.x);
                val_mse += val_loss.template item<float>();
                w_val_mse += w_val_loss.template item<float>();
                if (!isNet){
                    val_l_t += val_loss_total.template item<float>();
                    w_val_l_t += val_w_loss_total.template item<float>();
                }
                batch_idx++;
                val_batch_idx++;
            }
            val_mse /= (float) val_batch_idx ;
            w_val_mse /= (float) val_batch_idx;
            if (!isNet) {
                val_l_t /= (float) val_batch_idx;
                w_val_l_t /= (float) val_batch_idx;
            }
            printf("\t  %d %d %d %d \t ",true_positives, false_positives, true_negatives, false_negatives );
            if (!isNet){
                val_c_loss /= (float) val_batch_idx;
                printf("Non zeros in validation set -- %d & val cross-entropy loss -- %f", val_non_zeros,val_c_loss);
            }
            printf("Validation set error : %f", val_mse);
            //If the validation error is decreasing with training, update the model parameters
            if(0 == global_config.s_method.compare("is")){
                if(0 == global_config.train_stop.compare("mse")){
                    loss_to_compare = w_val_mse;
                }
                else if(0 == global_config.train_stop.compare("fn")) {
                    loss_to_compare = false_negatives;
                }
                else
                    loss_to_compare = w_val_l_t;
                if (loss_to_compare < prev_val_mse) {   // w_val_mse < prev_val_mse
                    //torch::save(model, "../best_model.pt");
                    best_mse = loss_to_compare ;           //best_mse = w_val_mse;
                    model = w_model;
                    epoch_t = epoch;
                    printf("Best model updated at epoch : %d",epoch);
                    mse_to_compare = loss_to_compare;
                    prev_val_mse=mse_to_compare;
                    count=0;
                }
                else{
                    count++;
                    if(count>global_config.stop_iter)
                        break;
                }
            }
            else{
                if (val_mse  < prev_val_mse) {
                    //torch::save(model, "../best_model.pt");
                    best_mse = val_mse ;
                    model = w_model;
                    epoch_t = epoch;
                    printf("Best model updated at epoch : %d",epoch);
                    mse_to_compare = val_mse ;
                    prev_val_mse=mse_to_compare;
                    count=0;
                }
                else{
                    count++;
                    if(count>global_config.stop_iter)
                        break;
                }
            }
            printf("Epoch number : %d, train mse : %f ", epoch, mse );
        }
        //out of training
        printf("Out of training --epoch num %d", epoch_t);
        auto elapsed = std::chrono::high_resolution_clock::now() - start;
        long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
        double duration = microseconds*1.0 / 1000000;
        duration /= n_epochs;
        this->WS()->time_Train += duration;
        this->WS()->count_Train++;
        std::string o_file = global_config.out_file + ".txt";
        std::ofstream to_write;
        to_write.open(o_file,std::ios_base::app);
        if (to_write.is_open()){
            printf("writing to the file %s", o_file.c_str());
            if (!isNet)
                to_write << global_config.width_problem << '\t'<< _nArgs << '\t' << epoch_t << '\t' << duration/3600 << '\t'<< global_config.loss_weight_mse << '\t' << global_config.loss_weight << '\t' << mse << '\t' << w_mse << '\t' << c_loss << '\t' << (float) val_mse << '\t' << (float)  w_val_mse << '\t' << val_c_loss<< '\t' ;
            else
                to_write << global_config.width_problem << '\t'<< _nArgs << '\t' << epoch_t << '\t' << duration/3600 << '\t' << mse << '\t' << w_mse << '\t' << val_mse << '\t' << w_val_mse << '\t';
            to_write.close();
        }
        global_config.avg_val_mse = (global_config.avg_val_mse*(this->WS()->count_Train - 1) + val_mse)/ this->WS()->count_Train;
        global_config.avg_val_w_mse = (global_config.avg_val_w_mse*(this->WS()->count_Train - 1) + w_val_mse)/ this->WS()->count_Train;
	}

    void test(DATA_SAMPLES *DS_test){
	    if(model==NULL)
	        printf("MODEL IS NOT TRAINED!!");
        bool isNet=false;
        if (0 == global_config.network.compare("net")) {
            isNet=true;
        }
        size_t batch_idx = 0;
        float mse = 0.,w_mse=0.,c_loss=0.0, non_zeros=0, test_non_zeros=0;
        auto dataset = DS_test->map(torch::data::transforms::Stack<>());
        int64_t batch_size = global_config.batch_size;
        torch::DeviceType device_type;
        if (torch::cuda::is_available()) {
            std::cout << "CUDA available! Training on GPU." << std::endl;
            device_type =  torch::kCUDA;//kCUDA
        } else {
            std::cout << "Training on CPU." << std::endl;
            device_type = torch::kCPU;
        }
        torch::Device device(device_type);
        auto data_loader = torch::data::make_data_loader<torch::data::samplers::RandomSampler>(dataset,batch_size);
        int32_t c=0;
        int total_batch =0, neg_batch =0, neg=0;
        int true_positives=0, false_positives=0, true_negatives=0, false_negatives=0;
        int true_positives_batch=0, false_positives_batch=0, true_negatives_batch=0, false_negatives_batch=0;
        for (auto &batch : *data_loader) {
            c++;
            torch::Tensor loss,w_loss, loss_n;
            auto imgs = batch.data;
            auto labels = batch.target.squeeze();
            imgs = imgs.to(device);
            labels = labels.to(device);
            auto output = model->forward(imgs,true,isNet);
            auto w = labels;
            float w1 = w.sum().template item<float>();
            if(w1!=0){
                w = w/w.sum();
            }
            if (!isNet) {
                auto label_binary = torch::where(labels > 0, torch::ones_like(labels), torch::zeros_like(labels));
                loss_n = torch::binary_cross_entropy(output.masked, label_binary);
                c_loss += loss_n.template item<float>();
                non_zeros = (label_binary==1).sum().template item<int>();
                test_non_zeros += non_zeros;
                neg_batch = (label_binary==0).sum().template item<int>();
                neg = neg + neg_batch;
                confusion_matrix(output.masked, label_binary, true_positives_batch, false_positives_batch, true_negatives_batch, false_negatives_batch);
                true_positives += true_positives_batch;
                false_positives += false_positives_batch;
                true_negatives += true_negatives_batch;
                false_negatives += false_negatives_batch;
            }
            w_loss = (w*((output.x - labels).pow(2))).mean();
            loss = torch::nn::functional::mse_loss(labels, output.x);
            mse += loss.template item<float>();
            w_mse += w_loss.template item<float>();
            batch_idx++;
        }

        mse /= (float) batch_idx;
        w_mse /= (float) batch_idx;
        if(!isNet)
            c_loss /= (float) batch_idx;
        printf("Test mse : %f ", mse );
        printf("\t Test: confusion matrix %d %d %d %d \t ",true_positives, false_positives, true_negatives, false_negatives);
        global_config.avg_test_mse = (global_config.avg_test_mse*(this->WS()->count_Train - 1) + mse)/ this->WS()->count_Train;
        global_config.avg_test_w_mse = (global_config.avg_test_w_mse*(this->WS()->count_Train - 1) + w_mse)/ this->WS()->count_Train;
        std::string o_file = global_config.out_file + ".txt";
        std::ofstream to_write;
        to_write.open(o_file,std::ios_base::app);
        if (to_write.is_open()){
            if(isNet)
                to_write << mse << '\t' << w_mse << '\t' << std::endl ;
            else
                to_write << mse << '\t'  << w_mse << '\t' << c_loss << '\t' << true_positives << '\t' << false_positives << '\t' << true_negatives << '\t' << false_negatives << '\t' << std::endl ;
            to_write.close();
        }
    }

    DATA_SAMPLES * samples_to_data(int32_t** samples_signiture, float* samples_values, int32_t input_size, int32_t sample_size, int &non_zeros)
    {
	    float zero = 0.1*pow(10,-50);
	    for(int32_t i=0; i<sample_size; i++){
            samples_values[i] = (float)samples_values[i]/max_value;
            if (samples_values[i]>=zero){
                non_zeros+=1;
            }
        }
        DATA_SAMPLES *DS;
        DS = new DATA_SAMPLES(samples_signiture, samples_values, input_size, sample_size);
        return DS;
    }

    void load_trained_model(){
    //we need to test the data here
    }

	FunctionNN(void)
		:
		Function()
	{
		// TODO own stuff here...
        printf("The void constuctor");
        _TableData = NULL;


	}

	FunctionNN(Workspace *WS, ARP *Problem, int32_t IDX)
	{
	    printf("I am in the Constructor here");
       Function(WS, Problem, IDX);
       _TableData = NULL;
	}

	virtual ~FunctionNN(void)
	{
		Destroy() ;
	}


} ;

/*
    torch::Tensor FunctionNN::weighted_mse_loss(auto input, auto target, auto weight) {
        torch::Tensor temp = weight*((input - target)^2);
        return temp.sum();
    }
*/
    inline FunctionNN * FunctionNNConstructor(void)
{
//    myData = empty_tensor.data_ptr<float>();
	return new FunctionNN;
}

} // namespace ARE

#endif // FunctionNN_HXX_INCLUDED
