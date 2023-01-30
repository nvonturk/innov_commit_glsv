%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_loss_function_wrapper.m
% Author: Nicholas von Turkovich
% Date: 12/16/2021
% Note(s): Wrapper function for the loss function to output the SMM
% criterion alone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loss] = compute_loss_function_wrapper(inputs, targets, hyper_params, weights, flag_entry, flag_strwrite, strwrite)
    
    [loss, ~] = compute_loss_function(inputs, targets, hyper_params, weights, flag_entry, flag_strwrite, strwrite);
        
end