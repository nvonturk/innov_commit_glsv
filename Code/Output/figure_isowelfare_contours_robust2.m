%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_isowelfare_contours_robust2.m
% Author: Nicholas von Turkovich
% Date: 06/17/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = figure_isowelfare_contours_robust2(filestart_data, econparams, commitment, consistent)

    % Adjustment versions
    adjustment_versions = [["phi0", "phi0_gamma0.7", "phi0_gamma0.7_arbmu"]];
        
    figs = {};
    
    for j = 1:size(adjustment_versions,1)

        % Adjustments
        adjustments = adjustment_versions(j,:);

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

        adjustments_clean = [strcat("$\phi = ", num2str(round(adj_econparams{1}.phi,2)), "$"),  ...
            strcat("$\phi = ", num2str(round(adj_econparams{1}.phi,2)), "$, ", "$\gamma = ", num2str(round(adj_econparams{2}.gamma,2)), "$"), ...
            strcat("$\phi = ", num2str(round(adj_econparams{3}.phi,2)), "$, ", "$\gamma = ", num2str(round(adj_econparams{3}.gamma,2)), "$, ", "$\mu_0(0) = ", num2str(adj_econparams_bgp{3}.mu(1)), "$")];

        f = figure('name', "isowelfare_contours_robust2_v" + num2str(j));

        for i = 1:length(adjustments)

            subplot(1,3,i)

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

            xlabel(["$\eta_1$"])

            if i == 1 
                ylabel(["$\eta_2$"])        
            end

            title(adjustments_clean(i))
            ylim([0,12])
            yticks(0:3:12)
            xlim([0,20])
            xticks([0:5:20])


        end
        
        figs{j} = f;
        
    end
    
    for j = 1:length(figs)
        f = figs{j};
        f.Units = 'centimeters';
        orient(f, 'landscape')
        f.PaperPosition = [-1, -1, 10, 4];
        saveas(f, "../../Output/Figures/Main/Figure_" + f.Name + ".png")
    end
    
    

end



