%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: satisfy_hjb.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vdiff] = satisfy_hjb(econparams, vnew, x_vector_overload)
    
    %% HJB Equations
    
    econparams.value_function = vnew;
    econparams = compute_innovation_decision(econparams);

    if nargin > 2
        vdiff = vnew - compute_hjb(econparams, x_vector_overload);
    else
        vdiff = vnew - compute_hjb(econparams);
    end
    
end
