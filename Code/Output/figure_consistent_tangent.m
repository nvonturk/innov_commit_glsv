%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_consistent_tangent.m
% Author: Nicholas von Turkovich
% Date: 06/1/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_consistent_tangent(econparams, commitment, consistent)

    f = figure("name", "consistent_tangent");
    
    % Show only local neighborhood around consistent duple; cutoff_value
    % gives the lowest eta1 value to be plotted
    cutoff_value = 10;

    % Extract the decomposition of welfare for the portion of the grid
    % above the cutoff value
    W_consistent_prod_interp_span = consistent.W_consistent_prod_interp(consistent.heta_grid_interp*100/econparams.frequency >= cutoff_value);
    W_consistent_muplab_interp_span = consistent.W_consistent_muplab_interp(consistent.heta_grid_interp*100/econparams.frequency >= cutoff_value);
    
    % Find the index of the optimal consistent duple and compute the slope
    % (decrease in welfare attributed to growth for an
    % increase in welfare attributed to markups/distortions)
    idx = find(consistent.W_consistent_interp == consistent.heta_optcons_duple_W);
    slope_cons = (consistent.W_consistent_prod_interp(idx + 1) - consistent.W_consistent_prod_interp(idx - 1))/(consistent.W_consistent_muplab_interp(idx + 1) - consistent.W_consistent_muplab_interp(idx - 1));

    % Similar calculation but for the policy with commitment
    idx = find(commitment.W_opt_interp == commitment.heta_opt_duple_W);
    slope_opt = (commitment.W_opt_prod_interp(idx + 1) - commitment.W_opt_prod_interp(idx - 1))/(commitment.W_opt_muplab_interp(idx + 1) - commitment.W_opt_muplab_interp(idx - 1));

    plot(W_consistent_muplab_interp_span, W_consistent_prod_interp_span, 'k-')
    hold on
    plot([consistent.heta_optcons_duple_W_muplab - 1, consistent.heta_optcons_duple_W_muplab + 1], [consistent.heta_optcons_duple_W_prod - 1*slope_cons, consistent.heta_optcons_duple_W_prod + 1*slope_cons], 'r-')
    plot(consistent.heta_optcons_duple_W_muplab, consistent.heta_optcons_duple_W_prod, 'k.', 'MarkerSize', 30)
    plot([commitment.heta_opt_duple_W_muplab + 4, commitment.heta_opt_duple_W_muplab + 8], [commitment.heta_opt_duple_W_prod + 4*slope_opt, commitment.heta_opt_duple_W_prod + 8*slope_opt],'r--')
    hold off
    grid off
    set(gca, 'xtick', [], 'ytick', [])
    xlabel('Markup distortions and disutility of labor')
    ylabel('Productivity index component')

    idx_start = length(W_consistent_muplab_interp_span)*0.2;
    idx_end = length(W_consistent_muplab_interp_span)*0.7;
    xoffset = -0.2;
    yoffset = -0.5;

    xcoords = W_consistent_muplab_interp_span(idx_start:idx_end) + xoffset;
    ycoords = W_consistent_prod_interp_span(idx_start:idx_end) + yoffset;

    hold on
    ha = annotation('arrow', 'color', 'red');
    ha.Parent = gca;
    ha.X = [consistent.heta_optcons_duple_W_muplab + 0.5, consistent.heta_optcons_duple_W_muplab + 2];
    ha.Y = [consistent.heta_optcons_duple_W_prod + 0.5, consistent.heta_optcons_duple_W_prod + 2];
    text(consistent.heta_optcons_duple_W_muplab, consistent.heta_optcons_duple_W_prod + 2, 'Increasing welfare', 'color', 'black')
    plot(xcoords, ycoords, 'k-')
    ha = annotation('arrow');
    arrow_slope = (ycoords(end) - ycoords(end-2))/(xcoords(end) - xcoords(end-2));
    ha.Parent = gca;
    ha.X = [xcoords(end), xcoords(end) + 0.1];
    ha.Y = [ycoords(end), ycoords(end) + 0.1*arrow_slope];
    text(xcoords(end)-1.5, ycoords(end)+0.3, ["Increasing initial patent", "expiry rate"], 'HorizontalAlignment', 'center')


    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



