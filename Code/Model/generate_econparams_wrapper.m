%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: generate_econparams_wrapper.m
% Author: Nicholas von Turkovich
% Date: 12/21/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = generate_econparams_wrapper(inputs, hyper_params, flag_entry, disp_flag)
    
    if flag_entry == 0
        econparams = generate_econparams([1/12, hyper_params(3), .03, inputs(3), hyper_params(1), -1, inputs(2), 1, 0.2, 0.2, inputs(1), inputs(1), inputs(4), 0, hyper_params(2), 0, 1, flag_entry, 1, 5e-8, 1e-4], disp_flag);
    else
        econparams = generate_econparams([1/12, hyper_params(3), .03, inputs(3), hyper_params(1), -1, inputs(2), 1, 0.2, 0.2, inputs(1), inputs(5), inputs(4), inputs(6), hyper_params(2), hyper_params(2), 1, flag_entry, 1, 5e-8, 1e-4], disp_flag);
    end
    
end