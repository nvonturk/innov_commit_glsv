%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: shooting_algorithm.m
% Author: Nicholas von Turkovich
% Date: 5/16/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_vector_path, x_vector_e_path, value_function_path, mu_path] = shooting_algorithm(econparams, value_function_end, mu_start, heta_path, rho_path, final_leg_start)
    
    %% Shoot back
  
    % Starting from the final BGP, we step backwards to determine the
    % innovation decisions for the firms
    [x_vector_path, x_vector_e_path, value_function_path] = shoot_back(econparams, ...
        value_function_end, heta_path, rho_path, final_leg_start);

    %% Shoot forward
    % Starting from the initial mu
    [mu_path] = shoot_forward(econparams, mu_start, heta_path, rho_path, x_vector_path, x_vector_e_path, final_leg_start);

end