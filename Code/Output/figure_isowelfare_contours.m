%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_isowelfare_contorus.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_isowelfare_contours(econparams, commitment, contour_mat, consistent, specification)

    f = figure("name", strcat("isowelfare_contours_", specification));

    % Extract contour data first before plotting
    [heta1_grid, heta2_grid] = ndgrid(contour_mat.heta_grid, contour_mat.heta_grid);
    
    % Extract contours for multiple isowelfare curves
    M = contour(heta1_grid, heta2_grid, contour_mat.W_grid);
    initial_contours = custom_contour(M, [contour_mat.heta_grid(1), contour_mat.heta_grid(end)], [commitment.heta_duple_commitment(2), contour_mat.heta_grid(end)]);
    
    % Extract contour specifically tangent to the optimal consistent duple
    M = contour(heta1_grid, heta2_grid, contour_mat.W_grid, [consistent.heta_optcons_duple_W, consistent.heta_optcons_duple_W]);
    tangent_contour = custom_contour(M, [contour_mat.heta_grid(1), contour_mat.heta_grid(end)], [commitment.heta_duple_commitment(2), contour_mat.heta_grid(end)]);

    % Plot the curve of consistent duple policies
    plot(consistent.heta_grid_interp*100/econparams.frequency, consistent.heta2_consistent_interp*100/econparams.frequency, 'k-')
    
    hold on
    
    % Plot the optimal consistent / optimal commitment policies
    plot(consistent.heta_optcons_duple(1)/econparams.frequency*100, consistent.heta_optcons_duple(2)/econparams.frequency*100, 'k.', 'MarkerSize', 30)
    plot(commitment.heta_duple_commitment(1)/econparams.frequency*100, commitment.heta_duple_commitment(2)/econparams.frequency*100, 'g.', 'MarkerSize', 30)

    for i = 1:length(initial_contours)

        c = initial_contours{i};
        plot(c(:,1)*100/econparams.frequency, c(:,2)*100/econparams.frequency, 'b-', 'LineWidth', 1);

    end

    leg = legend('Time-consistent duples', 'Optimal time-consistent policy', 'Optimal policy with commitment', 'Isowelfare contours');
    set(leg, 'AutoUpdate', 'off', 'color', 'white', 'edgecolor', 'white')
    c = tangent_contour{1};
    plot(c(:,1)*100/econparams.frequency, c(:,2)*100/econparams.frequency, 'b-', 'LineWidth', 1)
    hold off
    legend boxoff

    xlim([0, 20])
    xticks(0:2:20)
    xlabel("Initial patent expiry rate ($\eta_1$)")
    ylabel("Long-run patent expiry rate ($\eta_2$)")
    ylim([0, 10])
    grid off
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



