%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_isowelfare_contours_adjusting_t.m
% Author: Nicholas von Turkovich
% Date: 07/1/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = figure_isowelfare_contours_adjusting_t(filestart_data, econparams, commitment, consistent)

    % Adjustments
    adjustment = "T_series50";

    % Load necessary data
    adj_econparams = load(filestart_data + "_" + adjustment + "_econparams.mat");
    adj_econparams_bgp = load(filestart_data + "_" + adjustment + "_econparams_bgp.mat");
    adj_commitment = load(filestart_data + "_" + adjustment + "_commitment.mat");
    adj_contour_mat = load(filestart_data + "_" + adjustment + "_contour.mat");
    adj_consistent = load(filestart_data + "_" + adjustment + "_consistent.mat");
    
    adjustment_clean = strcat("$T = ", num2str(adj_econparams.interval_lengths(1)*adj_econparams.frequency), "$ years");

    f = figure('name', "isowelfare_contours_adjusting_t");

    % Lay down the contours first
    [heta1_grid, heta2_grid] = ndgrid(adj_contour_mat.heta_grid, adj_contour_mat.heta_grid);
    contour(heta1_grid*100/econparams.frequency, heta2_grid*100/econparams.frequency, adj_contour_mat.W_grid, 'DisplayName', 'Isowelfare contours');

    % Plot the baseline
    hold on
    plot(consistent.heta_grid_interp*100/econparams.frequency, consistent.heta2_consistent_interp*100/econparams.frequency, 'k-', 'DisplayName', 'Time-consistent duples')
    plot(consistent.heta_optcons_duple(1)*100/econparams.frequency, consistent.heta_optcons_duple(2)*100/econparams.frequency, 'k.', 'MarkerSize', 40, 'DisplayName', 'Optimal time-consistent policy')
    plot(commitment.heta_duple_commitment(1)*100/econparams.frequency, commitment.heta_duple_commitment(2)*100/econparams.frequency, 'g.', 'MarkerSize', 40, 'DisplayName', 'Optimal policy with commitment');

    % Plot the consistent values
    plot(adj_consistent.heta_grid_interp*100/econparams.frequency, adj_consistent.heta2_consistent_interp*100/econparams.frequency,  'k--', 'DisplayName', 'Time-consistent duples')
    plot(adj_consistent.heta_optcons_duple(1)*100/econparams.frequency, adj_consistent.heta_optcons_duple(2)*100/econparams.frequency, 'kx', 'MarkerSize', 15, 'DisplayName', 'Optimal time-consistent policy')
    plot(adj_commitment.heta_duple_commitment(1)*100/econparams.frequency, adj_commitment.heta_duple_commitment(2)*100/econparams.frequency, 'gx', 'MarkerSize', 15, 'DisplayName', 'Optimal policy with commitment');    
    hold off

    xlabel(["$\eta_1$"])

    ylabel(["$\eta_2$"])        

    title(adjustment_clean)
    ylim([0,12])
    yticks(0:3:12)
    xlim([0,20])
    xticks([0:5:20])

    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 10, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")    

end



