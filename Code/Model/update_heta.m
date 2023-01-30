%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: update_heta.m
% Author: Nicholas von Turkovich
% Date: 11/8/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = update_heta(econparams, heta, statedep_cutoff)
    
    if nargin > 2 && ~isempty(statedep_cutoff)

        econparams.heta = heta;
        econparams.heta_vector = ones(2*econparams.n + 1, 1);
        econparams.heta_vector(econparams.n + 1 - statedep_cutoff : econparams.n + 1 + statedep_cutoff) = heta(1);
        econparams.heta_vector(1 : econparams.n + 1 - statedep_cutoff - 1) = heta(2);
        econparams.heta_vector(econparams.n + 1 + statedep_cutoff + 1:end) = heta(2);
        econparams.heta_vector(econparams.n + 1) = 0;

    else 
        
        econparams.heta = heta;
        econparams.heta_vector = ones(2*econparams.n + 1, 1);
        econparams.heta_vector(econparams.n + 1) = 0;
        econparams.heta_vector = econparams.heta_vector * heta;
        
    end
        
        
end