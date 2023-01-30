%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: find_consistent_triple.m
% Author: Nicholas von Turkovich
% Date: 3/14/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [loss] = find_consistent_interp_triple(econparams, econparams_bgp, heta1, heta_intended, rho_changes, interval_lengths, heta_lb, heta_ub, heta_grid, options, triple_reopt_filename)

    tic;

    % First compute the full transition to determine mu_T
    transition = compute_transition(econparams, econparams_bgp, [heta1, heta_intended], rho_changes, interval_lengths);
        
    % Given the mu_T of this transition, replace the mu for the starting
    % point for the consistent duple calculation
    econparams_bgp_temp = econparams;
    econparams_bgp_temp.mu = transition.mu_path(:,interval_lengths(1)+1);
    econparams_bgp_temp.Q0 = transition.Q(interval_lengths(1)+1);
    econparams_bgp_temp.Y0 = transition.C(interval_lengths(1)+1);
    
    % Now given this triple and the resulting mu, we need to find an
    % optimal duple for SP2
    heta3_consistent = nan(size(heta_grid));
    W_consistent = nan(size(heta_grid));
    
    % Solve for the consistent schedule conditional on the starting mu
    % and choice for eta1
    
    % Load prior curves for consistent duple pairs
    log_file = readmatrix("Output/Logs/Main/" + triple_reopt_filename + "_curve.csv", 'NumHeaderLines', 1);
    
    % Check for row with closest overlap
    if size(log_file, 1) > 0
        
        [m, idx] = min(abs(log_file(:,1) - heta_intended(1)) + abs(log_file(:,2) - heta_intended(2)));
        heta3_guess = log_file(idx, 4:(3 + length(heta_grid)));

        if m < eps
            lookup = idx + 1;
            heta3_consistent = heta3_guess;
            W_consistent = log_file(idx,(3+length(heta_grid)+1):end);
            disp("Identical duple found; drawing from lookup table")
        else
            lookup = -1;
        end
    else
        heta3_guess = repmat(0.5*(heta_lb + heta_ub), 1, length(heta_grid));
        lookup = -1;
    end
    
    if lookup == -1
        
        parfor heta2_index = 1:length(heta_grid)
            
            [heta3_consistent(heta2_index), ~, ~, W_consistent(heta2_index)] = compute_consistent_duple(econparams, econparams_bgp_temp, heta_grid(heta2_index), heta3_guess(heta2_index), rho_changes(end-1:end), interval_lengths(end-1:end), heta_lb, heta_ub, options);  

        end

        writematrix([heta_intended, lookup, heta3_consistent, W_consistent], "Output/Logs/Main/" + triple_reopt_filename + "_curve.csv", 'WriteMode', 'append')
        disp("Line added to lookup table")
    end
    
    % Index of highest welfare consistent duple
    [~, idx_con] = max(W_consistent);
    if idx_con == 1
        heta2_interp = linspace(heta_grid(idx_con), heta_grid(idx_con+1), 100);
        warning("Optimal consistent duple on lower boundary")
    elseif idx_con == length(W_consistent)
        heta2_interp = linspace(heta_grid(idx_con-1), heta_grid(idx_con), 100);
        error("Optimal consistent duple on upper boundary")
    else
        heta2_interp = linspace(heta_grid(idx_con-1), heta_grid(idx_con+1), 200);
    end
    
    % Interpolate consistent values
    heta3_consistent_interp = spline(heta_grid, heta3_consistent, heta2_interp);
    W_consistent_interp = nan(1,length(heta2_interp));

    compute_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp_temp, heta_changes, rho_changes(end-1:end), interval_lengths(end-1:end));
    parfor idx = 1:length(W_consistent_interp)
       W_consistent_interp(idx) = -1*compute_transition_anon([heta2_interp(idx), heta3_consistent_interp(idx)]);
    end

    % Identify from interpolated values the optimal consistent policy
    [~, optcons_idx] = max(W_consistent_interp);
    heta2_max = heta2_interp(optcons_idx);
    heta3_max = heta3_consistent_interp(optcons_idx);

    loss = (abs(heta_intended(1) - heta2_max) + abs(heta_intended(2) - heta3_max))*100/econparams.frequency;
    dist2 = heta2_max - heta_intended(1);
    dist3 = heta3_max - heta_intended(2);
    
    % Save triple
    consistent_triple.heta1 = heta1;
    consistent_triple.heta2 = heta2_max;
    consistent_triple.heta3 = heta3_max;
    save("Data/Intermediate/" + triple_reopt_filename + "_consistent_triple.mat", "-struct", "consistent_triple");
    
    disp("========================")
    disp("Triple")
    fprintf("Intended: %2.5f, %2.5f, %2.5f\n", heta1*100/econparams.frequency, heta_intended(1)*100/econparams.frequency, heta_intended(2)*100/econparams.frequency)
    fprintf("Result: %2.5f, %2.5f, %2.5f\n", consistent_triple.heta1*100/econparams.frequency, consistent_triple.heta2*100/econparams.frequency, consistent_triple.heta3*100/econparams.frequency)
    fprintf("Loss: %2.5f\n", loss)
    runtime = toc/60
    writematrix([loss, dist2*100/econparams.frequency, dist3*100/econparams.frequency, runtime, heta1*100/econparams.frequency, heta_intended(1)*100/econparams.frequency, heta_intended(2)*100/econparams.frequency, ...
        consistent_triple.heta1*100/econparams.frequency, consistent_triple.heta2*100/econparams.frequency, consistent_triple.heta3*100/econparams.frequency], "Output/Logs/Main/" + triple_reopt_filename + ".csv", 'WriteMode', 'append')
    disp("========================")
    
    
end
