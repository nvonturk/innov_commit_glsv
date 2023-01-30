%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: value_function_iteration_robust.m
% Author: Nicholas von Turkovich
% Date: 10/25/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = value_function_iteration_robust(econparams, fhk_flag, x_vector_overload)
    
    % VFI
    if nargin > 2
        econparams = value_function_iteration_outer(econparams, x_vector_overload);
    else
        econparams = value_function_iteration_outer(econparams);
    end
    
    options = optimset('Display', 'off', 'MaxIter', 5e4, 'MaxFunEval', 5e4, ...
            'TolX', eps, 'TolFun', eps, 'Algorithm', 'trust-region-dogleg', ...
            'UseParallel', false);
        
    % Adjust value function 
    if nargin > 2
        [robust_value_function, fval1, flag1, output1] = fsolve(@(v) satisfy_hjb(econparams, v, x_vector_overload), econparams.value_function, options);
    else
        [robust_value_function, fval1, flag1, output1] = fsolve(@(v) satisfy_hjb(econparams, v), econparams.value_function, options);
    end

    % Recalculate the innovation decisions given final value function
    econparams.value_function_vfi = econparams.value_function;
    econparams.x_vector_vfi = econparams.x_vector;
    econparams.value_function = robust_value_function;
    econparams.fsolve_flag = flag1;
    econparams = compute_innovation_decision(econparams);

    % Compute the stationary distribution of states
    econparams = compute_stationary_distribution(econparams);
    
    if min(econparams.mu) < -1e-10
        error("Significantly negative mass in distribution");
    else
        econparams.mu = abs(econparams.mu)/(sum(abs(econparams.mu)));
    end
    
    % To allow interpolation to work in compute_key_moments
    econparams.mu(econparams.mu < eps) = eps;
    if abs(sum(econparams.mu) - 1) > 1e-10
        error("Adjustment to mu is signficiant")
    else
        econparams.mu = econparams.mu ./ sum(econparams.mu);
    end
    
    % Compute key moments (growth and profit share)
    econparams = compute_key_moments(econparams, fhk_flag);
    
    if econparams.disp_flag
        fprintf("VFI Completed, rho = %1.3f percent annually, eta = %1.3f percent annually\n", ...
            100*(power(1 + econparams.rho, 1/econparams.frequency)-1),...
            100*econparams.heta/econparams.frequency);
    end

end
