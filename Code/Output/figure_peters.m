%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_peters.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_peters(specification, econparams)
    
    % Extract calibration moments and targets with positive weight
    calibration_colnames = econparams.calibration_file.Properties.VariableNames;
    calibration_colnames = calibration_colnames(endsWith(calibration_colnames, "_w"));

    weights = table2array(econparams.calibration_file(:,calibration_colnames));

    moment_names = strrep(calibration_colnames, "_w", "");

    target_names = strcat(moment_names, "_target", "");
    target_values = table2array(econparams.calibration_file(1,target_names));

    param_names = ["B", "eta", "lambda", "phi"];
    param_names_plt = ["B", "$\eta$", "$\lambda$", "$\phi$"];

    if specification == "calibration_entry"
        param_names = [param_names, "B_e", "phi_e"];
        param_names_plt = [param_names_plt, "B_e", "$\phi_e$"];
    end

    param_values = table2array(econparams.calibration_file(1,param_names));

    moments_interest = ["tfp_growth", "rd_gdp", "mean_markup"];
    moments_interest_plt = ["$g^Q$", "R\&D to GDP", "Mean net markup"];
    [~,~,moments_interest_idx] = intersect(moments_interest, moment_names, 'stable');

    % Number of points to generate
    grid_pts = 100;
    B_grid = linspace(0.1,6,grid_pts);
    eta_grid = linspace(0,0.2,grid_pts);
    lambda_grid = linspace(1.001,1.05,grid_pts);
    phi_grid = linspace(0,1,grid_pts);
    B_e_grid = linspace(0.1,6,grid_pts);
    phi_e_grid = linspace(0,1,grid_pts);

    % Data structures to house loss/moments at each altered paameter
    % combination
    loss_grid = nan(length(param_names), grid_pts);
    moments_grid = nan(length(param_names), length(moment_names), grid_pts);

    % Anonymous function for computing loss/moments
    compute_loss_function_anon = @(param_v) compute_loss_function(param_v, target_values, [econparams.gamma, econparams.leap, econparams.n], weights, econparams.entry_flag, 0, "");

    % Calibration loss
    [loss_calibration, ~] = compute_loss_function_anon(param_values);

    % Loop through each parameter and fill in the loss/moments structures
    for i = 1:length(param_names)
        parfor j = 1:grid_pts

            param_temp = param_values;
            if i == 1
                param_temp(i) = B_grid(j);
            elseif i == 2
                param_temp(i) = eta_grid(j);
            elseif i == 3
                param_temp(i) = lambda_grid(j);
            elseif i == 4
                param_temp(i) = phi_grid(j);
            elseif i == 5
                param_temp(i) = B_e_grid(j);
            elseif i == 6
                param_temp(i) = phi_e_grid(j);
            else
                error("Set up entry")
            end

            [loss_grid(i,j), moments_grid(i,:,j)] = compute_loss_function_anon(param_temp);


        end

    end

    % Extract the moments of interest for the plot
    moments_grid_interest = moments_grid(:,moments_interest_idx, :);
    target_interest = target_values(moments_interest_idx);

    % Consolidate parameter grids
    param_grids = [B_grid; eta_grid; lambda_grid; phi_grid];

    f = figure('name', 'peters');
    hold on
    iter = 1;
    for j = 1:length(moments_interest_idx)

        for i = 1:length(param_names)

            subplot(length(moments_interest_idx) + 1, length(param_names), iter)

            plot(param_grids(i,:), squeeze(moments_grid_interest(i,j,:)), 'k-')

            xline(param_values(i), 'Color', [0.5, 0.5, 0.5])
            yline(target_interest(j),'r--')

            xlabel(param_names_plt(i))
            ylabel(moments_interest_plt(j))

            iter = iter + 1;

        end

    end

    for i = 1:length(param_names)
        subplot(length(moments_interest_idx) + 1, length(param_names), iter)
        plot(param_grids(i,:), loss_grid(i,:), 'k-')

        xline(param_values(i), 'Color', [0.5, 0.5, 0.5])
        yline(loss_calibration, 'r--')

        xlabel(param_names_plt(i))
        ylabel("Loss")

        iter = iter + 1;
    end

    hold off


end



