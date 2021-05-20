#!/bin/bash

#Parameters of our algorithm in general
s_size="500000"				# Bucket sample_size for training the NN
#induced width-> n_nuckets trained? 

# Hyper-parameters of the NN
learning_rate="0.0001 0.001"	# Optimization learning rate
n_epochs="10"					# Number of epochs to train the net
batch_size="64 128 256"			# Number of training sequences in each batch
n_hidden="100"				# Number of units in the hidden (recurrent) layer  


input_instance_folder="~/myresearch/GDBE/BucketElimination-mac-version/BE-sampling-project/3testproblems/"
instance_name="grid10x10.f10.wrap.uai" #pedigree20.uai"
input_ordering_folder="~/myresearch/GDBE/BucketElimination-mac-version/BE-sampling-project/3testproblems/"
ordering_file="grid10x10.f10.wrap_var_elim_order-StateSpaceSize.vo" #pedigree20_var_elim_order-StateSpaceSize.vo"

results_prefix="~/myresearch/GDBE/BucketElimination-mac-version/results/"


for b in $batch_size; do 
	for lr in $learning_rate; do
		./BESampling -fUAI $input_instance -fVO $input_ordering -iB 9999 -v2sample 369 -nsamples 10 --h_dim $n_hidden --batch_size $b --lr $lr --n_epochs $n_epochs --sample_size $s_size --out_file ${results_prefix}_${n_epochs}_${b}_${n_hidden}_${lr}
	done
done