//
// Created by yasaman razeghi
//


//this is the model architectuer for grid-> the simplest version
#ifndef BESAMPLING_NET_H
#define BESAMPLING_NET_H
#include <torch/nn/module.h>
#include <torch/torch.h>
#include <iostream>
#include <Config.h>


struct Masked_RET {
    torch::Tensor x;
    torch::Tensor masked;
};

struct Masked_Net : torch::nn::Module {
    int64_t h_dim = global_config.h_dim;
    torch::nn::Linear linear_in, linear1, linear2, linear3, linear4;
    torch::nn::Sequential linear_out, linear_mask, hidden;
   // torch::nn::ModuleList ;

    Masked_Net(int64_t Input_size)
            : linear_in(torch::nn::Linear(Input_size, h_dim)),
              linear_out(torch::nn::Sequential(torch::nn::Linear(global_config.h_dim,global_config.h_dim), torch::nn::ReLU(),torch::nn::Linear(global_config.h_dim,1))),
              linear_mask(torch::nn::Sequential(torch::nn::Linear(global_config.h_dim,global_config.h_dim), torch::nn::ReLU(),torch::nn::Linear(global_config.h_dim,1))),
              linear1(torch::nn::Linear(Input_size, h_dim)),
              linear2(torch::nn::Linear(h_dim, h_dim)),
              linear3(torch::nn::Linear(h_dim, 1)),
              linear4(torch::nn::Linear(h_dim, 1))
    {
        register_module("linear1", linear1);
        register_module("linear2", linear2);
        register_module("linear3", linear3);
        register_module("linear4", linear4);
        register_module("linear_in", linear_in);
        register_module("linear_out", linear_out);
        register_module("linear_mask",linear_mask);
    }

    Masked_RET forward(torch::Tensor input, bool isTest, bool net) {
        Masked_RET ret_Masked_Net;
        if(net){
            torch::Tensor x = torch::relu(linear1(input));
            x =  torch::relu(linear2(x));
            x = linear3(x);
            ret_Masked_Net.x = x.view(-1);
            return ret_Masked_Net;
        }
        else{
            torch::Tensor x = torch::relu(linear_in(input));
            torch::Tensor mask = torch::sigmoid(linear4(x));
            x = torch::relu(linear2(x));
            x = torch::nn::functional::softplus(linear3(x));
            if (isTest == true){
                mask = torch::floor(mask+0.5);
                ret_Masked_Net.x = x*mask.view(-1);
                ret_Masked_Net.masked = mask.view(-1);
                return ret_Masked_Net;
            }
            ret_Masked_Net.x = x.view(-1);
            ret_Masked_Net.masked = mask.view(-1);
            return ret_Masked_Net;
        }
    }
};
//TORCH_MODULE(Net);

#endif //BESAMPLING_NET_H
