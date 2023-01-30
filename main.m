%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: main.m
% Author: Nicholas von Turkovich
% Date: 4/19/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%% Setup

% Generate directories

folders = ["Data", "Data/Raw", "Data/Intermediate", ...
    "Output", "Output/Tables", "Output/Figures", "Output/Logs", "Output/Tables", ...
    "Output/Figures/Main", "Output/Figures/Robustness",...
    "Output/Logs/Main", "Output/Logs/Calibration_Deterministic", "Output/Logs/Robustness", ...
    "Output/Logs/Calibration_Deterministic/Compiled"];

for i = 1:length(folders)
    if ~isfolder(folders(i))
        mkdir(folders(i))
    end
end

% Paths to necessary code

% This path will need to point to the Matlab module for your version of
% NLOPT
nlopt_path = "/home/software/nlopt/2.7.1-matlab/lib64/matlab/";

try
    addpath(nlopt_path)
catch
    warning("NLOPT library not linked. If intending to run calibration code," + ...
        " will need to add NLOPT library to path.")
end

addpath(genpath("./Code/Model"))
addpath(genpath("./Code/Calibration"))
addpath(genpath("./Code/Output"))

% Setup parallel pool

% Determine number of threads available and launch parallel pool; note that
% this assumes code is being run in a SLURM environment. If not the
% try-catch will assign the parallel cluster number of workers
sz = str2num(getenv('SLURM_NPROCS'));

c = parcluster;
if ~isempty(sz)
    c.NumWorkers = sz;
else
    sz = c.NumWorkers;
end

disp(c.NumWorkers);

try
    pool = parpool(sz, 'IdleTimeout', Inf);
    fprintf("Pool size is %d\n", pool.NumWorkers)
catch
    disp("Pool exists")
end

%% Replication switches

% Plotting

% Generate plots
switch_plots = true;

% Policy exercises

% Baseline calibration (order of hours)
switch_baseline = false;

% Pre-sclerosis calibration (order of hours)
switch_presclerosis = false;

% State dependent policy (order of a day)
switch_statedep = false;

% Double reoptimization (order of days given the nesting optimizations)
switch_triple = false;

% Growth decomposition (order of hours - 1 day)
switch_decomposition = false;

% Time pattern of profits (order of hours with high memory requirements
% given simulation averages over many random draws)
switch_profits = false;

% Adjusting baseline calibration parameters to study time-consistency
% problem (order of hours)
switch_robust = false;

% Adjusting the length of the initial period T (order of hours)
switch_T = false;

% Calibration (time on order of days to converge)

% Run baseline calibration
switch_calibration = false;

% Run pre-sclerosis calibration
switch_calibration_presclerosis = false;

%% Code execution

% Baseline calibration (output used in majority of figures and tables)
if switch_baseline

    disp("switch_baseline");

    switches = [true  ;... % Optimal static BGP patent expiry rate
                true  ;... % Optimal static patent expiry rate incorporating transition dynamics
                true  ;... % Duples with commitment
                false ;... % Duples with commitment and state dependent policy 
                false ;... % Optimal triple with commitment
                true  ;... % Welfare of duples and second social planner's best response
                true  ;... % Time-consistent duples
                false ;... % Time-consistent duples with state dependent policy
                false ;... % Time-consistent triples
                false];    % Interpolation of time-consistent triples

    param_file = "calibration";

    heta_ub_dec = 0.2; % Upper bound on search space for eta (annual decimal)
    flag_use_arbitrary_mu = false; % Flag for overriding initial BGP mu with perfectly competitive vector
    var = ""; % Type of parameter adjustment
    adjustments = []; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
    
end

% Presclerosis calibration
if switch_presclerosis

    disp("switch_presclerosis");

    switches = [true  ;... % Optimal static BGP patent expiry rate
                true  ;... % Optimal static patent expiry rate incorporating transition dynamics
                true  ;... % Duples with commitment
                false ;... % Duples with commitment and state dependent policy 
                false ;... % Optimal triple with commitment
                true  ;... % Welfare of duples and second social planner's best response
                true  ;... % Time-consistent duples
                false ;... % Time-consistent duples with state dependent policy
                false ;... % Time-consistent triples
                false];    % Interpolation of time-consistent triples

    param_file = "calibration_presclerosis";

    heta_ub_dec = 0.2; % Upper bound on search space for eta (annual decimal)
    flag_use_arbitrary_mu = false; % Flag for overriding initial BGP mu with perfectly competitive vector
    var = ""; % Type of parameter adjustment
    adjustments = []; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
    
end

% State dependent policy exercises
if switch_statedep

    disp("switch_statedep");

    switches = [false ;... % Optimal static BGP patent expiry rate
                false ;... % Optimal static patent expiry rate incorporating transition dynamics
                false ;... % Duples with commitment
                true  ;... % Duples with commitment and state dependent policy 
                false ;... % Optimal triple with commitment
                false ;... % Welfare of duples and second social planner's best response
                false ;... % Time-consistent duples
                true  ;... % Time-consistent duples with state dependent policy
                false ;... % Time-consistent triples
                false];    % Interpolation of time-consistent triples

    param_file = "calibration";

    heta_ub_dec = 0.4; % Upper bound on search space for eta (annual decimal)
    flag_use_arbitrary_mu = false; % Flag for overriding initial BGP mu with perfectly competitive vector
    var = ""; % Type of parameter adjustment
    adjustments = []; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
    
end

% Exercises with double reoptimization
if switch_triple

    disp("switch_triple");

    switches = [false  ;... % Optimal static BGP patent expiry rate
                false  ;... % Optimal static patent expiry rate incorporating transition dynamics
                false  ;... % Duples with commitment
                false  ;... % Duples with commitment and state dependent policy 
                true   ;... % Optimal triple with commitment
                false  ;... % Welfare of duples and second social planner's best response
                false  ;... % Time-consistent duples
                false  ;... % Time-consistent duples with state dependent policy
                true   ;... % Time-consistent triples
                true];    % Interpolation of time-consistent triples

    param_file = "calibration";

    heta_ub_dec = 0.2; % Upper bound on search space for eta (annual decimal)
    flag_use_arbitrary_mu = false; % Flag for overriding initial BGP mu with perfectly competitive vector
    var = ""; % Type of parameter adjustment
    adjustments = []; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
    
end

% Adjusting parameters to study the nature of the time consistency problem
if switch_robust

    disp("switch_robust");

    var_vec = ["phi", "B", "rho", "", "gamma", "phi_gamma"];
    adjustments_vec = {0.25, [0.12], [0.015], [], [0.5], [0; 2/3]};
    arbmu_vec = [false, false, false, true, false, true];
    
    for i = 1:length(var_vec)
    
        switches = [false ;... % Optimal static BGP patent expiry rate
                    false ;... % Optimal static patent expiry rate incorporating transition dynamics
                    true  ;... % Duples with commitment
                    false ;... % Duples with commitment and state dependent policy 
                    false ;... % Optimal triple with commitment
                    true  ;... % Welfare of duples and second social planner's best response
                    true  ;... % Time-consistent duples
                    false ;... % Time-consistent duples with state dependent policy
                    false ;... % Time-consistent triples
                    false];    % Interpolation of time-consistent triples

        param_file = "calibration";

        heta_ub_dec = 0.2; % Upper bound on search space for eta (annual decimal)
        flag_use_arbitrary_mu = arbmu_vec(i); % Flag for overriding initial BGP mu with perfectly competitive vector
        var = var_vec(i); % Type of parameter adjustment
        adjustments = adjustments_vec{i}; % Vector containing the parameter adjustments associated with 'var'

        policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
        
    end
    
end

% Varying the initial patent expiry period T
if switch_T
 
    disp("switch_T");

    % First for T = 50 years with the contour exercise
    switches = [false ;... % Optimal static BGP patent expiry rate
                false ;... % Optimal static patent expiry rate incorporating transition dynamics
                true  ;... % Duples with commitment
                false ;... % Duples with commitment and state dependent policy 
                false ;... % Optimal triple with commitment
                true  ;... % Welfare of duples and second social planner's best response
                true  ;... % Time-consistent duples
                false ;... % Time-consistent duples with state dependent policy
                false ;... % Time-consistent triples
                false];    % Interpolation of time-consistent triples

    param_file = "calibration";

    heta_ub_dec = 0.2; % Upper bound on search space for eta (annual decimal)
    flag_use_arbitrary_mu = false; % Flag for overriding initial BGP mu with perfectly competitive vector
    var = "T_series"; % Type of parameter adjustment
    adjustments = 50; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
    
    % Second, iterate for series of T values without contour plots
    switches = [false ;... % Optimal static BGP patent expiry rate
                false ;... % Optimal static patent expiry rate incorporating transition dynamics
                true  ;... % Duples with commitment
                false ;... % Duples with commitment and state dependent policy 
                false ;... % Optimal triple with commitment
                false ;... % Welfare of duples and second social planner's best response
                true  ;... % Time-consistent duples
                false ;... % Time-consistent duples with state dependent policy
                false ;... % Time-consistent triples
                false];    % Interpolation of time-consistent triples

    heta_ub_dec = 0.4; % Upper bound on search space for eta (annual decimal)
    adjustments = [10, 15, 20, 25, 30, 35, 50]; % Vector containing the parameter adjustments associated with 'var'

    policy_exercises(switches, param_file, heta_ub_dec, flag_use_arbitrary_mu, var, adjustments);
            
end

% Exercises relying on baseline exercises

% Growth decomposition in response to a small change in the terminal expiry
% rate
if switch_decomposition

    disp("switch_decomposition");

    param_file = "calibration";
    transition_type = "consistent";
    
    strategic_decomposition(param_file, transition_type);
    
end

% Time pattern of profits associated with an innovation
if switch_profits

    disp("switch_profits");

    param_file = "calibration";
    transition_type = "consistent";
    
    simulate_profits_experiment(param_file, transition_type);
    
end

% Plots
if switch_plots
    
    disp("switch_plots");

    plot_global_settings();
    param_file = "calibration";
    
    plot_core(param_file);
    
end

% Calibration exercises

if switch_calibration

    disp("switch_calibration");

    flag_entry = 0;
    flag_presclerosis = 0;
    
    % Calibrate baseline model
    calibration_deterministic(flag_entry, flag_presclerosis);
    
    % Compile the results of the calibration
    system('Rscript --vanilla Code/Calibration/calibration_compile.R --args calibration');
    
end

if switch_calibration_presclerosis

    disp("switch_presclerosis");

    flag_entry = 0;
    flag_presclerosis = 1;
    
    % Calibrate pre-sclerosis model
    calibration_deterministic(flag_entry, flag_presclerosis);
    
    % Compile the results of the calibration
    system('Rscript --vanilla Code/Calibration/calibration_compile.R --args calibration_presclerosis');
    
end



