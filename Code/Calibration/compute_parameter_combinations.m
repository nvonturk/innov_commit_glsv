%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_parameter_combinations.m
% Author: Nicholas von Turkovich
% Date: 1/25/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [hyper_parameters] = compute_parameter_combinations(varargin)
    
    hyper_parameters = [varargin{1}];
    
    for i = 2:nargin
        
        moments_temp = [];
        
        for j = 1:size(varargin{i},1)
            
            row = [hyper_parameters, ones(size(hyper_parameters,1),size(varargin{i},2)).*varargin{i}(j,:)];
            moments_temp = [moments_temp; row];
            
        end
        
        hyper_parameters = moments_temp;
        
    end
    
end
