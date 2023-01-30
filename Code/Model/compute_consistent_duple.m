%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: compute_consistent_duple.m
% Author: Nicholas von Turkovich
% Date: 11/17/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [heta2_consistent, heta2_max, mu_policy_change, W_schedule] = compute_consistent_duple(econparams, econparams_bgp, heta1, heta2_consistent, rho_changes, interval_lengths, heta_lb, heta_ub, options, statedep_cutoff)
    
    % Initialize distance between the consistent policy guess and the
    % optimal response from the second social planner
    dist = 10;
    
    % Initialize the bounds on the search space for the consistent policy
    heta_lb_bisection = heta_lb;
    heta_ub_bisection = heta_ub;
    
    % Instantiate some additional structures to track the progress of the
    % bisection
    iter = 1;
    heta2_maxes = [];
    heta2_cons = [];
    heta_lbs = [];
    heta_ubs = [];
    dists = [];
    
    % Continue the search until the consistent policy guess and the second
    % social planner's best response are sufficiently close
    while abs(dist) > options.OptimalityTolerance
       
        % Collect data from the last guess
        if iter > 1
            heta2_maxes = [heta2_maxes, heta2_max];
            heta2_cons = [heta2_cons, heta2_consistent];
            heta_lbs = [heta_lbs, heta_lb_bisection];
            heta_ubs = [heta_ubs, heta_ub_bisection];
            dists = [dists, dist];
        end
        
        if iter > 1
            
            if dist < 0
                % Guess was higher than response to shrink the search space
                % from above
                heta_ub_bisection = heta2_consistent;
            else
                % Guess was lower than response so shrink the search space
                % from below
                heta_lb_bisection = heta2_consistent;
            end 
            
            if iter > 2
                
                % For additional speed, attempt to more accurately guess
                % the consistent solution using splines and break if the
                % spline fails; can only use spline for guess after first
                % two iterations
                try
                    heta2_consistent = spline(dists, heta2_cons, 0);
                catch
                    test.dists = dists;
                    test.heta2_cons = heta2_cons;
                    test.heta2_maxes = heta2_maxes;
                    test.heta_lbs = heta_lbs;
                    test.heta_ubs = heta_ubs;
                    test.econparams_bgp = econparams_bgp;
                    test.heta1 = heta1;
                    save("../Data/Intermediate/test.mat", "-struct", "test");
                    disp("Hit spline error")
                end
                
                % Safeguards to ensure the guess using splines falls within
                % the valid search space, guesses are not too close to the
                % search bounds, and that consecutive guesses are not too
                % close together
                if heta2_consistent > heta_ub_bisection || heta2_consistent < heta_lb_bisection
                    heta2_consistent = mean([heta_lb_bisection, heta_ub_bisection]);
                elseif (heta_ub_bisection - heta_lb_bisection) < options.OptimalityTolerance*10
                    heta2_consistent = mean([heta_lb_bisection, heta_ub_bisection]);
                elseif abs(heta2_consistent - heta2_cons(end)) < options.OptimalityTolerance
                    heta2_consistent = mean([heta_lb_bisection, heta_ub_bisection]);
                end
                                
            else
                heta2_consistent = mean([heta_lb_bisection, heta_ub_bisection]);
            end
                
        end
                    
        % Determine distance from optimal choice for second planner given
        % guess
        if nargin == 10 && ~isempty(statedep_cutoff)
            [dist, heta2_max, mu_policy_change, W_schedule] = find_consistent_heta(econparams, econparams_bgp, heta1, heta2_consistent, rho_changes, interval_lengths, heta_lb_bisection, heta_ub_bisection, options, statedep_cutoff);
        else
            [dist, heta2_max, mu_policy_change, W_schedule] = find_consistent_heta(econparams, econparams_bgp, heta1, heta2_consistent, rho_changes, interval_lengths, heta_lb_bisection, heta_ub_bisection, options);
        end
        
        iter = iter + 1;
    end
                
    if abs(heta2_max - heta_ub)*100/econparams.frequency < 1
        fprintf("Warning: consistent solution close to upper bound of %2.2f\n", heta_ub*100/econparams.frequency)
        error("Halting code given consistent solution found artificially")
    end

end
