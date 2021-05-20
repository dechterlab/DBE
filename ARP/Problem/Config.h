
#ifndef BESAMPLING_CONFIG_H
#define BESAMPLING_CONFIG_H
//Here I create this struct to set the parameters:
typedef struct Config{

        //NN hyper-parameters
        int32_t sample_size = 500000;
        float lr = 0.001;
        int32_t n_epochs = 100;
        int32_t dev_size = 10000;
        std::string network="masked_net";
        std::string s_method="us";
        float loss_weight = 1;
        float loss_weight_mse =0.5;
        int num_hidden = 2;
        float dropout = 0.5;
        int32_t h_dim = 256;
        int32_t batch_size = 256;   //1024
        std::string train_stop="mse";
        int stop_iter=1;

        //i-bound of the problem
        int width_problem = 20;

        //output files
        std::string out_file, out_file2;

        //evaluation metrics for NN
        float avg_val_mse=0.0,avg_test_mse=0.0, avg_val_w_mse=0.0,avg_test_w_mse=0.0 ;


}Config_NN;

extern Config_NN global_config;

#endif //BESAMPLING_CONFIG_H
