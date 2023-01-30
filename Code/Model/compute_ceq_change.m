%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_ceq_change.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi] = compute_ceq_change(W_from, W_to, rho)
    
    % CGLS A5
    xi = (exp((W_to - W_from)*rho)) - 1;
    
end
