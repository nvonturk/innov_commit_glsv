%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_isowelfare_contours_robust.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_isowelfare_contours_robust(filestart_data, econparams, commitment, consistent)

    % Adjustments
    adjustments = ["phi0.25", "B0.12", "rho0.015", "arbmu", "gamma0.5", "phi0_gamma0.7_arbmu"];

    % Load necessary data for each plot
    adj_econparams = {};
    adj_econparams_bgp = {};
    adj_commitment = {};
    adj_contour_mat = {};
    adj_consistent = {};

    for i = 1:length(adjustments)
        adj_econparams{i} = load(filestart_data + "_" + adjustments(i) + "_econparams.mat");
        adj_econparams_bgp{i} = load(filestart_data + "_" + adjustments(i) + "_econparams_bgp.mat");
        adj_commitment{i} = load(filestart_data + "_" + adjustments(i) + "_commitment.mat");
        adj_contour_mat{i} = load(filestart_data + "_" + adjustments(i) + "_contour.mat");
        adj_consistent{i} = load(filestart_data + "_" + adjustments(i) + "_consistent.mat");
        
    end

    % Note that we report the effective B (i.e. assuming no R&D subsidy)
    adjustments_clean = [strcat("$\phi = ", num2str(round(adj_econparams{1}.phi,2)), "$"),  ...
        strcat("$B = ", num2str(round(adj_econparams{2}.B/((1-adj_econparams{2}.tau_rd)^(adj_econparams{2}.gamma))/adj_econparams{2}.frequency,2)), "$"), ...
        strcat("$\rho = ", num2str(round((power(adj_econparams{3}.rho + 1, 1/adj_econparams{3}.frequency)-1)*100,2)), "$"), ...
        strcat("$\mu_0(0) = ", num2str(adj_econparams_bgp{4}.mu(1)), "$"),  ...
        strcat("$\gamma = ", num2str(adj_econparams{5}.gamma), "$"), ...
        strcat("$\phi = ", num2str(round(adj_econparams{6}.phi,2)), "$, ", "$\gamma = ", num2str(round(adj_econparams{6}.gamma,2)), "$, ", "$\mu_0(0) = ", num2str(adj_econparams_bgp{6}.mu(1)), "$")];

    f = figure('name', 'isowelfare_contours_robust');

    for i = 1:length(adjustments)

        subplot(2,3,i)

        % Lay down the contours first
        [heta1_grid, heta2_grid] = ndgrid(adj_contour_mat{i}.heta_grid, adj_contour_mat{i}.heta_grid);
        contour(heta1_grid*100/econparams.frequency, heta2_grid*100/econparams.frequency, adj_contour_mat{i}.W_grid, 'DisplayName', 'Isowelfare contours');

        % Plot the baseline
        hold on
        plot(consistent.heta_grid_interp*100/econparams.frequency, consistent.heta2_consistent_interp*100/econparams.frequency, 'k-', 'DisplayName', 'Time-consistent duples')
        plot(consistent.heta_optcons_duple(1)*100/econparams.frequency, consistent.heta_optcons_duple(2)*100/econparams.frequency, 'k.', 'MarkerSize', 40, 'DisplayName', 'Optimal time-consistent policy')
        plot(commitment.heta_duple_commitment(1)*100/econparams.frequency, commitment.heta_duple_commitment(2)*100/econparams.frequency, 'g.', 'MarkerSize', 40, 'DisplayName', 'Optimal policy with commitment');

        % Plot the consistent values
        plot(adj_consistent{i}.heta_grid_interp*100/econparams.frequency, adj_consistent{i}.heta2_consistent_interp*100/econparams.frequency,  'k--', 'DisplayName', 'Time-consistent duples')
        plot(adj_consistent{i}.heta_optcons_duple(1)*100/econparams.frequency, adj_consistent{i}.heta_optcons_duple(2)*100/econparams.frequency, 'kx', 'MarkerSize', 15, 'DisplayName', 'Optimal time-consistent policy')
        plot(adj_commitment{i}.heta_duple_commitment(1)*100/econparams.frequency, adj_commitment{i}.heta_duple_commitment(2)*100/econparams.frequency, 'gx', 'MarkerSize', 15, 'DisplayName', 'Optimal policy with commitment');    
        hold off

        if i+1 > 4
            xlabel(["$\eta_1$"])
        end

        if i == 1 || i+1 == 5
            ylabel(["$\eta_2$"])        
        end

        title(adjustments_clean(i))
        ylim([0,12])
        yticks(0:3:12)
        xlim([0,20])
        xticks([0:5:20])


    end
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 10, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")

end



