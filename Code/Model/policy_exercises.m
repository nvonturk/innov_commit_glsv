%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: policy_exercises.m
% Author: Nicholas von Turkovich
% Date: 11/16/2021
% Note(s): Core set of transition, commitment, etc. exercises
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, ...
    var, adjustments)

    %% Determine experiments to run

    % Switches for different experiments
    switch_static_bgp = switches(1);
    switch_static_transition = switches(2);
    switch_duple_commitment = switches(3);
    switch_duple_commitment_statedep = switches(4);
    switch_triple_commitment = switches(5);
    switch_contour = switches(6);
    switch_consistent = switches(7);
    switch_consistent_statedep = switches(8);
    switch_interp_triple = switches(9);
    switch_spline_triple = switches(10);

    %% Generate model parameters for experiments

    if param_file == "calibration_entry"
        flag_entry = 1;
    else
        flag_entry = 0;
    end

    alg = "NLOPT_GN_DIRECT"; % Algorithm used in NLOPT run

    if param_file == "calibration_presclerosis"
        calib = 12; % Lowest loss calibration in the pre-sclerosis economy
    else
        calib = 23; % Lowest loss calibration in the baseline economy
    end

    % Read in the calibration file
    if flag_entry
        calibrations = readtable("Output/Tracked/master_compiled_calibration_entry.csv");
    elseif param_file == "calibration"
        calibrations = readtable("Output/Tracked/master_compiled_calibration.csv");
    elseif param_file == "calibration_presclerosis"
        calibrations = readtable("Output/Tracked/master_compiled_calibration_presclerosis.csv");
    else
        error("Unrecognized calibration specification");
    end

    % Identify the specific parameter set from the calibration
    specification = strcmp(alg, calibrations.algorithm) & ...
        calibrations.calibration == calib;
    current_calibration = calibrations(specification,:);

    % Choice of maximum gap to overload function
    n_alt = 25;

    % Depending on the type of run, choose the econparams and store the
    % calibration output for reference
    if param_file ~= "calibration_entry"
        econparams = generate_econparams_wrapper([current_calibration.B, current_calibration.eta, current_calibration.lambda, current_calibration.phi], [current_calibration.gamma, current_calibration.leap, n_alt], flag_entry, 0);
        econparams.calibration_file = current_calibration;
    else
        econparams = generate_econparams_wrapper([current_calibration.B, current_calibration.eta, current_calibration.lambda, current_calibration.phi, current_calibration.B_e, current_calibration.phi_e], [current_calibration.gamma, current_calibration.leap, n_alt], flag_entry, 0);
        econparams.calibration_file = current_calibration;
    end

    %% Setup auxiliary variables needed to run experiments

    % Bounds for eta in setting up the grid as well as search space
    heta_lb_dec = 0; % Value in annual decimal

    heta_ub = heta_ub_dec*econparams.frequency; % Scaled to solution frequency
    heta_lb = heta_lb_dec*econparams.frequency; % Scaled to solution frequency

    % Grid for set of policies under commitment, time-consistent policy, and
    % the contour plots
    heta_grid = linspace(heta_lb, heta_ub, 100*(heta_ub_dec - heta_lb_dec) + 1);

    % Grid interp dendsity multiplier used to get finer grid through
    % interpolation by spline
    grid_interp_density = 100;

    % Tolerance for how close to ub on search space to throw error (to avoid
    % misinterpretation of solution when just on boundary of solver)
    heta_dist_tol = 0.01*econparams.frequency; % Annual decimal (then scaled)

    % Initial guess for heta in optimization exercises
    heta_guess = 0.01*econparams.frequency; % Annual decimal (then scaled)

    % Fmincon parameters (these also are used in determining convergence in
    % finding consistent schedules, see compute_consistent_duple.m)
    heta_opt_tol = 1e-8;
    heta_step_tol = 1e-11;
    options = optimoptions('fmincon', ...
        'OptimalityTolerance', heta_opt_tol, ...
        'StepTolerance', heta_step_tol,...
        'Display', 'off',...
        'Algorithm', 'sqp');

    % Transition intervals which determine the amount of time for each leg of
    % the transition
    short_leg = 25; % Years
    final_leg = 100; % Years
    econparams.interval_lengths = [short_leg, short_leg, final_leg]/econparams.frequency; % Scaled to solution frequency

    % State dependency cutoff for policy transitions
    statedep_cutoff = 1;

    %% Adjustments to core exercise

    % Var indicates the type of exogenous parameter change we are interested in
    % (note that if var = "" that is our baseline for the majority of analysis)
    if var == ""
        disp("Adjustment is: baseline");
    else
        disp("Adjustment is: " + var);
    end

    % If we are experimenting with adjustments to a parameter, generate an
    % parameter file for each adjustment
    if var ~= ""

        % Save the econparams objects and tag for the adjustment
        econparams_vec = cell(1,size(adjustments,2));
        suffix_vec = cell(1,size(adjustments,2));

        for i = 1:size(adjustments,2)

            econparams_temp = econparams;

            if  var == "phi"
                econparams_temp = update_phi(econparams_temp, adjustments(i));
                suffix_temp = var + round(econparams_temp.phi,2) + "_";
            elseif var == "phi_gamma"
                econparams_temp = update_econparams(econparams_temp, ["phi", "gamma"], adjustments(:,i));
                suffix_temp = "phi" + round(econparams_temp.phi,2) + "_gamma" + round(econparams_temp.gamma, 1) + "_";
            elseif var == "T_series"
                medium_horizon = sum(econparams.interval_lengths(1:2));
                econparams_temp.interval_lengths(1:2) = [adjustments(i)/econparams.frequency, medium_horizon - adjustments(i)/econparams.frequency];
                suffix_temp = var + round(econparams_temp.interval_lengths(1)*econparams.frequency, 0) + "_";
            elseif var == "T_series_alt"
                econparams_temp.interval_lengths(1:2) = [adjustments(i), adjustments(i)];
                suffix_temp = var + round(econparams_temp.interval_lengths(1)*econparams.frequency, 0) + "_";
            elseif var == "rho"
                econparams_temp = update_econparams(econparams, "rho", adjustments(i));
                suffix_temp = var + round(econparams_temp.rho/econparams.frequency, 3) + "_";
            elseif var == "gamma"
                econparams_temp = update_econparams(econparams, "gamma", adjustments(i));
                suffix_temp = var + round(econparams_temp.gamma, 1) + "_";
            elseif var == "B"
                econparams_temp = update_econparams(econparams, "B", adjustments(i));
                suffix_temp = var + round(econparams_temp.B/econparams.frequency, 2) + "_";
            end

            if flag_use_arbitrary_mu
                suffix_temp  = suffix_temp + "arbmu_";
            end

            econparams_vec{i} = econparams_temp;
            suffix_vec{i} = suffix_temp;

        end

    else

        econparams_vec = {econparams};

        if flag_use_arbitrary_mu
            suffix_vec = {"arbmu_"};
        else
            suffix_vec = {""};
        end

    end

    %% Exercises for each economy

    for i = 1:length(econparams_vec)

        %% Save econparams for each adjustment
        fprintf("Iteration %1.0f\n", i);

        econparams = econparams_vec{i};
        suffix = suffix_vec{i};

        if heta_ub_dec ~= 0.2 % Default search space 
            suffix = strcat(suffix, "hub", num2str(100*heta_ub_dec), "_");
        end

        if n_alt ~= 25 % Default choice for maximum technology gap
            suffix = strcat(suffix, "n", num2str(n_alt), "_");
        end

        save("Data/Intermediate/" + param_file + "_" + suffix + "econparams" + ".mat", "-struct", "econparams")
        disp("Saved econparams")   

        %% Solve for the BGP and initialize mu

        % Solve for the BGP for the calibration
        econparams_bgp = value_function_iteration_robust(econparams, false);

        % Option to make the initial mu distribution highly competitive
        if flag_use_arbitrary_mu
            mu = zeros(size(econparams_bgp.mu));
            mu(1) = 1;
            econparams_bgp.mu = mu;
            econparams_bgp.Y0 = econparams_bgp.Q0; % When all firms are in the tied state, there are no markups
        end

        save("Data/Intermediate/" + param_file + "_" + suffix + "econparams_bgp" + ".mat", "-struct", "econparams_bgp")
        disp("Saved econparams_bgp")

        %% Optimal BGP eta

        if switch_static_bgp

            tic

            % Determine the BGP with optimal eta
            optimal_bgp_anon = @(heta) compute_optimal_bgp(econparams, heta);
            [heta_static_bgp, ~] = fmincon(@(heta) optimal_bgp_anon(heta), [heta_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);

            if (heta_static_bgp + heta_dist_tol) > heta_ub
                error("Hitting boundary on search for eta")
            end

            save("Data/Intermediate/" + param_file + "_" + suffix + "heta_static_bgp" + ".mat", "heta_static_bgp")

            fprintf("Optimal BGP Done. Elapsed time is %1.1f minutes.\n", toc/60);

        end

        %% Optimal static transition

        if switch_static_transition

            tic

            optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp, heta_changes, econparams.rho, sum(econparams.interval_lengths));
            [heta_static_transition, ~] = fmincon(@(heta_changes) optimal_transition_anon(heta_changes), [heta_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);

            if (heta_static_transition + heta_dist_tol) > heta_ub
                error("Hitting boundary on search for eta")
            end

            save("Data/Intermediate/" + param_file + "_" + suffix + "heta_static_transition" + ".mat", "heta_static_transition")

            fprintf("Optimal Static Transition Done. Elapsed time is %1.1f minutes.\n", toc/60);

        end

        %% Optimal duple

        if switch_duple_commitment

            tic

            optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);
            [heta_duple_commitment, ~] = fmincon(@(heta_changes) optimal_transition_anon(heta_changes), repmat([heta_guess],1,2), [], [], [], [], repmat([heta_lb],1,2), repmat([heta_ub],1,2), [], options);

            if sum((heta_duple_commitment + heta_dist_tol) > heta_ub) > 0
                error("Hitting boundary on search for eta")
            end

            % Store the optimal duple in a struct and set up the grid for
            % optimal commitment policy fixing eta1
            commitment.heta_duple_commitment = heta_duple_commitment;
            commitment.heta_grid = heta_grid;

            heta2_opt = nan(size(heta_grid));

            % Solve for the optimal choice of eta2 given a fixed eta1 (with
            % commitment)
            parfor i = 1:length(heta_grid)
                [heta2_opt_temp, ~] = fmincon(@(heta_changes) optimal_transition_anon([heta_grid(i), heta_changes]), heta_guess, [], [], [], [], heta_lb, heta_ub, [], options);
                heta2_opt(i) = heta2_opt_temp;
            end

            commitment.heta2_opt = heta2_opt;

            % Interpolate values to form a finer grid of optimal duples with
            % commitment
            heta_grid_interp = linspace(min(heta_grid), max(heta_grid), grid_interp_density*length(heta_grid));
            heta2_opt_interp = spline(heta_grid, heta2_opt, heta_grid_interp);

            % While we interpolate to solve for the commitment policies,
            % compute the transition for each interpolated duple to arrive at
            % the welfare (compute welfare from scratch)
            W_opt_interp = nan(size(heta_grid_interp));
            W_opt_prod_interp = nan(size(heta_grid_interp));
            W_opt_muplab_interp = nan(size(heta_grid_interp));
            compute_transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

            parfor i = 1:length(heta_grid_interp)
                transition_temp = compute_transition_anon([heta_grid_interp(i), heta2_opt_interp(i)]);
                W_opt_interp(i) = transition_temp.W;
                W_opt_prod_interp(i) = transition_temp.W_prod;
                W_opt_muplab_interp(i) = transition_temp.W_muplab;
            end

            % Store the interpolated grid and policy values
            commitment.heta_grid_interp = heta_grid_interp;
            commitment.heta2_opt_interp = heta2_opt_interp;
            commitment.W_opt_interp = W_opt_interp;
            commitment.W_opt_prod_interp = W_opt_prod_interp;
            commitment.W_opt_muplab_interp = W_opt_muplab_interp;

            % Compute the optimal policy from the interpolated grid (should be
            % close to the solution solved globally by fmincon)
            [heta_opt_duple_W,idx] = max(W_opt_interp);
            commitment.heta_opt_duple = [heta_grid_interp(idx), heta2_opt_interp(idx)];
            commitment.heta_opt_duple_W = heta_opt_duple_W;
            commitment.heta_opt_duple_W_prod = W_opt_prod_interp(idx);
            commitment.heta_opt_duple_W_muplab = W_opt_muplab_interp(idx);

            % Given the optimal duple with commitment...
            duple_commitment_transition = compute_transition_anon([commitment.heta_duple_commitment]);

            % ...determine the state of the economy at time T...
            econparams_T = econparams;
            econparams_T.mu = duple_commitment_transition.mu_path(:,econparams.interval_lengths(1)+1);
            econparams_T.Q0 = duple_commitment_transition.Q(econparams.interval_lengths(1) + 1);
            econparams_T.Y0 = duple_commitment_transition.C(econparams.interval_lengths(1) + 1);

            % ... and solve the second social planner's problem
            optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_T, heta_changes, econparams.rho, sum(econparams.interval_lengths(end-1:end)));
            [commitment.sp2_best_response, ~] = fmincon(@(heta_changes) optimal_transition_anon(heta_changes), [heta_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);

            % Repeat systematically for the grid of optimal duples with
            % commitment
            sp2_best_response_grid = nan(size(heta_grid));

            optimal_transition_anon = @(econparams_bgp_changes, heta_changes) compute_optimal_transition(econparams, econparams_bgp_changes, heta_changes, econparams.rho, sum(econparams.interval_lengths(end-1:end)));

            parfor i = 1:length(heta_grid)

                econparams_T_temp = econparams;

                tx_temp = compute_transition_anon([heta_grid(i), heta2_opt(i)]);

                econparams_T_temp.mu = tx_temp.mu_path(:,econparams.interval_lengths(1)+1);
                econparams_T_temp.Q0 = tx_temp.Q(econparams.interval_lengths(1) + 1);
                econparams_T_temp.Y0 = tx_temp.C(econparams.interval_lengths(1) + 1);

                [sp2_best_response_grid(i), ~] = fmincon(@(heta_changes) optimal_transition_anon(econparams_T_temp, heta_changes), [heta_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);

            end

            commitment.sp2_best_response_grid = sp2_best_response_grid;

            fprintf("Optimal Duple Done. Elapsed time is %1.1f minutes.\n", toc/60);
            save("Data/Intermediate/" + param_file + "_" + suffix + "commitment" + ".mat", "-struct", "commitment")


        end

        %% Optimal duple with state dependent policy

        if switch_duple_commitment_statedep

            tic

            optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], statedep_cutoff);
            [heta_duple_commitment_statedep, ~] = fmincon(@(heta_changes) optimal_transition_anon([[heta_changes(1); heta_changes(2)],[heta_changes(3); heta_changes(3)]]), repmat([heta_guess],1,3), [], [], [], [], repmat([heta_lb],1,3), repmat([heta_ub],1,3), [], options);

            if sum((heta_duple_commitment_statedep + heta_dist_tol) > heta_ub) > 0
                heta_duple_commitment_statedep*100/econparams.frequency
                error("Hitting boundary on search for eta")
            end
            
            fprintf("Completed optimal state dependent commitment policy. Elapsed time is %1.1f minutes.\n", toc/60);

            % Store the optimal duple in a struct and set up the grid for
            % optimal commitment policy fixing eta1
            commitment_statedep.heta_duple_commitment_statedep = heta_duple_commitment_statedep;
            commitment_statedep.heta_grid = heta_grid;

            heta2_opt_statedep = nan(length(heta_grid), length(heta_grid));

            % Solve for the optimal choice of eta2 given a fixed eta1 for
            % competitive and uncompetitive states
            for i = 1:length(heta_grid)
                tic
                parfor j = 1:length(heta_grid)
                    [heta2_opt_temp_statedep, ~] = fmincon(@(heta_changes) optimal_transition_anon([[heta_grid(i); heta_grid(j)], [heta_changes; heta_changes]]), heta_guess, [], [], [], [], heta_lb, heta_ub, [], options);
                    heta2_opt_statedep(i,j) = heta2_opt_temp_statedep;
                end
                fprintf("Completed loop iteration %1.1f in grid. Elapsed time is %1.1f minutes.\n", i, toc/60);
            end

            disp("Completed double for loop for grid of state dependent" + ...
                " commitment duples")
            tic

            commitment_statedep.heta2_opt_statedep = heta2_opt_statedep;

            % Interpolate values to form a finer grid of optimal duples with
            % commitment
            heta_grid_interp = linspace(min(heta_grid), max(heta_grid), grid_interp_density/10*length(heta_grid));

            [heta_grid_comp, heta_grid_uncomp] = ndgrid(heta_grid, heta_grid);
            [heta_grid_interp_comp, heta_grid_interp_uncomp] = ndgrid(heta_grid_interp, heta_grid_interp);

            % Requires meshgrid so after interpolating flip the output back
            heta2_opt_interp_statedep = interp2(heta_grid_comp', heta_grid_uncomp', heta2_opt_statedep', heta_grid_interp_comp', heta_grid_interp_uncomp', 'spline');
            heta2_opt_interp_statedep = heta2_opt_interp_statedep';

            % While we interpolate to solve for the commitment policies,
            % compute the transition for each interpolated duple to arrive at
            % the welfare (compute welfare from scratch)
            W_opt_interp_statedep = nan(length(heta_grid_interp), length(heta_grid_interp));
            compute_transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], [], statedep_cutoff);

            for i = 1:length(heta_grid_interp)
                parfor j = 1:length(heta_grid_interp)
                    transition_temp = compute_transition_anon([[heta_grid_interp(i); heta_grid_interp(j)], repmat(heta2_opt_interp_statedep(i,j), 2, 1)]);
                    W_opt_interp_statedep(i,j) = transition_temp.W;
                end
            end

            % Store the interpolated grid and policy values
            commitment_statedep.heta_grid_interp = heta_grid_interp;
            commitment_statedep.heta2_opt_interp_statedep = heta2_opt_interp_statedep;
            commitment_statedep.W_opt_interp_statedep = W_opt_interp_statedep;

            % Compute the optimal policy from the interpolated grid (should be
            % close to the solution solved globally by fmincon)
            [heta_opt_duple_W_statedep,idx] = max(W_opt_interp_statedep(:));
            [idx_i, idx_j] = ind2sub(size(W_opt_interp_statedep), idx);
            commitment_statedep.heta_opt_duple_statedep = [heta_grid_interp(idx_i), heta_grid_interp(idx_j), heta2_opt_interp_statedep(idx_i,idx_j)];
            commitment_statedep.heta_opt_duple_W_statedep = heta_opt_duple_W_statedep;

            fprintf("Commitment statedependency interpolation done. Elapsed time is %1.1f minutes.\n", toc/60);
            save("Data/Intermediate/" + param_file + "_" + suffix + "commitment_statedep" + num2str(statedep_cutoff) + ".mat", "-struct", "commitment_statedep")


        end

        %% Optimal triple

        if switch_triple_commitment

            tic

            optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,3), econparams.interval_lengths(end-2:end));
            [heta_triple_commitment, ~] = fmincon(@(heta_changes) optimal_transition_anon(heta_changes), repmat([heta_guess],1,3), [], [], [], [], repmat([heta_lb],1,3), repmat([heta_ub],1,3), [], options);

            if sum((heta_triple_commitment + heta_dist_tol) > heta_ub) > 0
                error("Hitting boundary on search for eta")
            end

            save("Data/Intermediate/" + param_file + "_" + suffix + "heta_triple_commitment" + ".mat", "heta_triple_commitment")

            fprintf("Optimal Triple Done. Elapsed time is %1.1f minutes.\n", toc/60);


        end

        %% Contour plot

        if switch_contour

            tic

            % Instantiate matrices for the welfare of each duple with
            % commitment, the mu at the point of transition, and the optimal
            % static transition decision of the second social planner given the
            % duple chosen by the first social splanner
            W_grid = nan(length(heta_grid));
            heta2_opt = nan(size(heta_grid));
            mu_grid = nan(length(heta_grid), length(heta_grid), length(econparams_bgp.mu));
            sp2_static = nan(length(heta_grid));

            % Anonymous transition function
            compute_transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

            % Double for loop to determine welfare of each duple in the grid
            for heta1_index = 1:length(heta_grid)

                heta1 = heta_grid(heta1_index);

                parfor heta2_index = 1:length(heta_grid)

                    heta2 = heta_grid(heta2_index);

                    % Store welfare of transition
                    transition_temp = compute_transition_anon([heta1, heta2]);
                    W_grid(heta1_index, heta2_index) = transition_temp.W;

                    % Store the mu at the change over in social planners
                    transition_points = cumsum(econparams.interval_lengths(end-2:end))+1;
                    mu_change = transition_temp.mu_path(:,transition_points(1));
                    mu_grid(heta1_index, heta2_index, :) = mu_change;

                    % Separate computation - given the mu resulting from this
                    % path, what would social planner 2 choose to do
                    econparams_bgp_temp = econparams_bgp;
                    econparams_bgp_temp.mu = mu_change;
                    econparams_bgp_temp.Q0 = transition_temp.Q(transition_points(1));
                    econparams_bgp_temp.Y0 = transition_temp.C(transition_points(1));

                    optimal_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp_temp, heta_changes, econparams.rho, sum(econparams.interval_lengths(end-1:end)));
                    [sp2_static(heta1_index, heta2_index), ~] = fmincon(@(heta_changes) optimal_transition_anon(heta_changes), [heta_guess], [], [], [], [], [heta_lb], [heta_ub], [], options);

                end

                [~, idx] = max(W_grid(heta1_index, :));
                heta2_opt(heta1_index) = heta_grid(idx);

            end

            % Store the structures
            contour_mat.W_grid = W_grid;
            contour_mat.heta2_opt = heta2_opt;
            contour_mat.heta_grid = heta_grid;
            contour_mat.mu_grid = mu_grid;
            contour_mat.sp2_static = sp2_static;

            save("Data/Intermediate/" + param_file + "_" + suffix + "contour" + ".mat", "-struct", "contour_mat")

            fprintf("Contour Done. Elapsed time is %1.1f minutes.\n", toc/60);

        end

        %% Dynamic policies without commitment (consistent policies)

        if switch_consistent

            tic

            % Create matrices for the consistent choices for eta 2, the mu on
            % change over, the max welfare choice of the second social planner
            % (should be within the tolerance level of the consistent policy),
            % and the welfare associated with the consistent duple
            heta2_consistent = nan(size(heta_grid));
            mu_policy_change = nan(length(econparams_bgp.mu), length(heta_grid));
            heta2_max = nan(size(heta_grid));
            W_consistent = nan(size(heta_grid));

            % Solve for the consistent schedule conditional on the starting mu
            % and choice for eta1
            parfor heta1_index = 1:length(heta_grid)

                [heta2_consistent(heta1_index), heta2_max(heta1_index), mu_policy_change(:,heta1_index), W_consistent(heta1_index)] = compute_consistent_duple(econparams, econparams_bgp, heta_grid(heta1_index), heta_guess, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], heta_lb, heta_ub, options);  

            end

            % Store consistent values without interpolation
            consistent.heta_grid = heta_grid;
            consistent.heta2_consistent = heta2_consistent;
            consistent.mu_policy_change = mu_policy_change;
            consistent.heta2_max = heta2_max;
            consistent.W_consistent = W_consistent;

            % Interpolate consistent values
            consistent.heta_grid_interp = linspace(min(heta_grid), max(heta_grid), grid_interp_density*length(heta_grid));
            consistent.heta2_consistent_interp = spline(consistent.heta_grid, consistent.heta2_consistent, consistent.heta_grid_interp);
            W_consistent_interp = nan(1,length(consistent.heta_grid_interp));
            W_consistent_prod_interp = nan(1,length(consistent.heta_grid_interp));
            W_consistent_muplab_interp = nan(1,length(consistent.heta_grid_interp));

            compute_transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);
            parfor idx = 1:length(W_consistent_interp)
               transition_temp = compute_transition_anon([consistent.heta_grid_interp(idx), consistent.heta2_consistent_interp(idx)]);
               W_consistent_interp(idx) = transition_temp.W;
               W_consistent_prod_interp(idx) = transition_temp.W_prod;
               W_consistent_muplab_interp(idx) = transition_temp.W_muplab;
            end

            consistent.W_consistent_interp = W_consistent_interp;
            consistent.W_consistent_prod_interp = W_consistent_prod_interp;
            consistent.W_consistent_muplab_interp = W_consistent_muplab_interp;

            % Identify from interpolated values the optimal consistent policy
            [heta_optcons_duple_W, optcons_idx] = max(consistent.W_consistent_interp);
            heta_optcons_duple = [consistent.heta_grid_interp(optcons_idx), consistent.heta2_consistent_interp(optcons_idx)];

            % Save the max interpolated value
            consistent.heta_optcons_duple = heta_optcons_duple;
            consistent.heta_optcons_duple_W = heta_optcons_duple_W;
            consistent.heta_optcons_duple_W_prod = consistent.W_consistent_prod_interp(optcons_idx);
            consistent.heta_optcons_duple_W_muplab = consistent.W_consistent_muplab_interp(optcons_idx);

            % Compute a check on the interpolation by solving for a sample of
            % the interpolated heta1 values

            heta_grid_interp_check = consistent.heta_grid_interp(grid_interp_density/2:grid_interp_density:length(consistent.heta_grid_interp));
            heta2_consistent_interp_check = nan(size(heta_grid_interp_check));

            parfor heta1_index = 1:length(heta_grid_interp_check)

                [heta2_consistent_interp_check(heta1_index), ~, ~, ~] = compute_consistent_duple(econparams, econparams_bgp, heta_grid_interp_check(heta1_index), heta_guess, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], heta_lb, heta_ub, options);  

            end

            consistent.heta_grid_interp_check = heta_grid_interp_check;
            consistent.heta2_consistent_interp_check = heta2_consistent_interp_check;

            save("Data/Intermediate/" + param_file + "_" + suffix + "consistent" + ".mat", "-struct", "consistent")

            fprintf("Consistent Duples Done. Elapsed time is %1.1f minutes.\n", toc/60);

        end

        if switch_consistent_statedep

            tic

            % Create matrices for the consistent choices for eta 2, the mu on
            % change over, the max welfare choice of the second social planner
            % (should be within the tolerance level of the consistent policy),
            % and the welfare associated with the consistent duple
            heta2_consistent_statedep = nan(length(heta_grid), length(heta_grid));
            mu_policy_change_statedep = nan(length(econparams_bgp.mu), length(heta_grid), length(heta_grid));
            heta2_max_statedep = nan(length(heta_grid), length(heta_grid));
            W_consistent_statedep = nan(length(heta_grid), length(heta_grid));

            % Solve for the consistent schedule conditional on the starting mu
            % and choice for eta1
            for heta1_index_comp = 1:length(heta_grid)
                tic
                parfor heta1_index_uncomp = 1:length(heta_grid)

                    [heta2_consistent_statedep(heta1_index_comp, heta1_index_uncomp), heta2_max_statedep(heta1_index_comp, heta1_index_uncomp), mu_policy_change_statedep(:,heta1_index_comp, heta1_index_uncomp), W_consistent_statedep(heta1_index_comp, heta1_index_uncomp)] = compute_consistent_duple(econparams, econparams_bgp, [heta_grid(heta1_index_comp); heta_grid(heta1_index_uncomp)], heta_guess, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], heta_lb, heta_ub, options, statedep_cutoff);  

                end
                fprintf("Completed loop iteration %1.1f in grid. Elapsed time is %1.1f minutes.\n", i, toc/60);
            end
            
            disp("Completed double for loop for grid of state dependent" + ...
                " consistent duples")
            tic

            % Store consistent values without interpolation
            consistent_statedep.heta_grid = heta_grid;
            consistent_statedep.heta2_consistent_statedep = heta2_consistent_statedep;
            consistent_statedep.mu_policy_change_statedep = mu_policy_change_statedep;
            consistent_statedep.heta2_max_statedep = heta2_max_statedep;
            consistent_statedep.W_consistent_statedep = W_consistent_statedep;

            % Interpolate consistent values
            consistent_statedep.heta_grid_interp = linspace(min(heta_grid), max(heta_grid), grid_interp_density/10*length(heta_grid));

            [heta_grid_comp, heta_grid_uncomp] = ndgrid(heta_grid, heta_grid);
            [heta_grid_interp_comp, heta_grid_interp_uncomp] = ndgrid(consistent_statedep.heta_grid_interp, consistent_statedep.heta_grid_interp);
            heta2_consistent_interp_statedep = interp2(heta_grid_comp', heta_grid_uncomp', consistent_statedep.heta2_consistent_statedep', heta_grid_interp_comp', heta_grid_interp_uncomp', 'spline');
            consistent_statedep.heta2_consistent_interp_statedep = heta2_consistent_interp_statedep';

            W_consistent_interp_statedep = nan(length(consistent_statedep.heta_grid_interp),length(consistent_statedep.heta_grid_interp));

            compute_transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))], [], statedep_cutoff);

            for i = 1:length(consistent_statedep.heta_grid_interp)
                parfor j = 1:length(consistent_statedep.heta_grid_interp)
                   transition_temp = compute_transition_anon([[consistent_statedep.heta_grid_interp(i); consistent_statedep.heta_grid_interp(j)], [consistent_statedep.heta2_consistent_interp_statedep(i,j); consistent_statedep.heta2_consistent_interp_statedep(i,j)]]);
                   W_consistent_interp_statedep(i,j) = transition_temp.W;
                end
            end

            consistent_statedep.W_consistent_interp_statedep = W_consistent_interp_statedep;

            % Identify from interpolated values the optimal consistent policy
            [heta_optcons_duple_W_statedep, optcons_idx] = max(consistent_statedep.W_consistent_interp_statedep(:));
            [optcons_idx_i, optcons_idx_j] = ind2sub(size(consistent_statedep.W_consistent_interp_statedep), optcons_idx);
            consistent_statedep.heta_optcons_duple_statedep = [consistent_statedep.heta_grid_interp(optcons_idx_i), consistent_statedep.heta_grid_interp(optcons_idx_j), consistent_statedep.heta2_consistent_interp_statedep(optcons_idx_i, optcons_idx_j)];
            consistent_statedep.heta_optcons_duple_W_statedep = heta_optcons_duple_W_statedep;

            save("Data/Intermediate/" + param_file + "_" + suffix + "consistent_statedep" + num2str(statedep_cutoff) + ".mat", "-struct", "consistent_statedep")

            fprintf("Consistent state-dependent interpolation done. Elapsed time is %1.1f minutes.\n", toc/60);

        end

        %% Dynamic triple without commitment

        if switch_interp_triple

            for i = 1:2:length(heta_grid)

                heta_grid_idx = i;

                % Name of triple reopt file
                triple_reopt_filename = param_file + "_" + suffix + "tripconslog";

                % Set the bounds on the search space to the "factory"
                % settings
                heta2_lb = heta_lb;
                heta2_ub = heta_ub;
                heta3_lb = heta_lb;
                heta3_ub = heta_ub;

                % Create filename for the log
                triple_reopt_filename = triple_reopt_filename + "_eta" + round(heta_grid(heta_grid_idx)*100/econparams.frequency,2);

                % Create the log for this particular value of eta1
                writematrix(["loss", "dist2", "dist3", "runtime", "heta1_intended", "heta2_intended", "heta3_intended", "heta1", "heta2", "heta3"], "Output/Logs/Main/" + triple_reopt_filename + ".csv"); 
                if exist("Output/Logs/Main/" + triple_reopt_filename + "_curve.csv", 'file') ~= 2
                    writematrix(["heta2_int", "heta3_int", "look_up", repmat("heta3_con", 1, length(heta_grid)), repmat("Wcon", 1, length(heta_grid))], "Output/Logs/Main/" + triple_reopt_filename + "_curve.csv"); 
                end

                % Generate anonymous function for finding the consistent triple
                find_consistent_interp_triple_anon = @(heta1, heta2, heta3) find_consistent_interp_triple(econparams, econparams_bgp, heta1, [heta2, heta3], repmat(econparams.rho,1,3), econparams.interval_lengths, heta_lb, heta_ub, heta_grid, options, triple_reopt_filename);

                % Start point for search
                heta_region_start = 0.5*([heta2_lb, heta3_lb] + [heta2_ub, heta3_ub]);

                [heta_min, loss] = fmincon(@(heta_intended) find_consistent_interp_triple_anon(heta_grid(heta_grid_idx), heta_intended(1), heta_intended(2)), heta_region_start, [], [], [], [], [heta2_lb, heta3_lb], [heta2_ub, heta3_ub], [], options);
                
            end

        end

        if switch_spline_triple

            tic;

            % Name of triple reopt file
            triple_reopt_filename = param_file + "_" + suffix + "tripconslog";
            system('Rscript Code/Model/tripcons_compile.R ' + triple_reopt_filename);
            compiled = readtable("Output/Logs/Main/" + triple_reopt_filename + "_compiled.csv");

            compute_transition_anon = @(heta_changes) compute_optimal_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,3), econparams.interval_lengths(end-2:end));

            % Load data from compiled file
            consistent_triple_interp.heta1 = compiled.heta1_intended/100*econparams.frequency;
            consistent_triple_interp.heta2_consistent = compiled.heta2_intended/100*econparams.frequency;
            consistent_triple_interp.heta3_consistent = compiled.heta3_intended/100*econparams.frequency;

            W_consistent_triple = nan(1,length(consistent_triple_interp.heta1));

            parfor idx = 1:length(W_consistent_triple)
               W_consistent_triple(idx) = -1*compute_transition_anon([consistent_triple_interp.heta1(idx), consistent_triple_interp.heta2_consistent(idx), consistent_triple_interp.heta3_consistent(idx)]);
            end

            consistent_triple_interp.W_consistent_triple = W_consistent_triple;

            % Interpolate consistent values
            consistent_triple_interp.heta1_interp = linspace(min(consistent_triple_interp.heta1), max(consistent_triple_interp.heta1), grid_interp_density*length(heta_grid));
            consistent_triple_interp.heta2_consistent_interp = spline(consistent_triple_interp.heta1, consistent_triple_interp.heta2_consistent, consistent_triple_interp.heta1_interp);
            consistent_triple_interp.heta3_consistent_interp = spline(consistent_triple_interp.heta1, consistent_triple_interp.heta3_consistent, consistent_triple_interp.heta1_interp);        

            W_consistent_triple_interp = nan(1,length(consistent_triple_interp.heta1_interp));

            parfor idx = 1:length(W_consistent_triple_interp)
               W_consistent_triple_interp(idx) = -1*compute_transition_anon([consistent_triple_interp.heta1_interp(idx), consistent_triple_interp.heta2_consistent_interp(idx), consistent_triple_interp.heta3_consistent_interp(idx)]);
            end

            consistent_triple_interp.W_consistent_interp = W_consistent_triple_interp;

            % Identify from interpolated values the optimal consistent policy
            [heta_optcons_duple_W, optcons_idx] = max(consistent_triple_interp.W_consistent_interp);
            heta_optcons_triple = [consistent_triple_interp.heta1_interp(optcons_idx), consistent_triple_interp.heta2_consistent_interp(optcons_idx), consistent_triple_interp.heta3_consistent_interp(optcons_idx)];

            % Save the max interpolated value
            consistent_triple_interp.heta_max_optcons = heta_optcons_triple;
            consistent_triple_interp.heta_max_optcons_W = heta_optcons_duple_W;

            save("Data/Intermediate/" + param_file + "_" + suffix + "tripcons" + ".mat", "-struct", "consistent_triple_interp")

            fprintf("Consistent Triples Done. Elapsed time is %1.1f minutes.\n", toc/60);


        end

    end
    
end

