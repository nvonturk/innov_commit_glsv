%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_duple_contour.m
% Author: Nicholas von Turkovich
% Date: 05/26/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_duple_contour(econparams, heta_static_bgp, heta_static_transition, commitment, contour_mat, consistent, ...
    nochange_ceq, bgp_ceq, static_ceq, duple_commitment_ceq, optcons_duple_ceq, consistent_duple_ceq, specification)
    
    f = figure('name', strcat("duple_contour_", specification));
    
    % Create baseline contour plot
    [heta1_grid, heta2_grid] = ndgrid(contour_mat.heta_grid/econparams.frequency*100, contour_mat.heta_grid/econparams.frequency*100);
    contour(heta1_grid, heta2_grid, contour_mat.W_grid)

    % Overlay the policies of interest
    hold on
    plot(econparams.heta/econparams.frequency*100, econparams.heta/econparams.frequency*100, 'kx')
    plot(heta_static_bgp/econparams.frequency*100, heta_static_bgp/econparams.frequency*100, 'rx')
    plot(heta_static_transition/econparams.frequency*100, heta_static_transition/econparams.frequency*100, 'b*')
    plot(commitment.heta_grid*100/econparams.frequency, commitment.heta2_opt/econparams.frequency*100, 'go')
    plot(commitment.heta_grid_interp*100/econparams.frequency, commitment.heta2_opt_interp*100/econparams.frequency, 'g-')
    plot(commitment.heta_duple_commitment(1)/econparams.frequency*100, commitment.heta_duple_commitment(2)/econparams.frequency*100, 'g.', 'MarkerSize', 30)
    plot(consistent.heta_grid/econparams.frequency*100, consistent.heta2_consistent/econparams.frequency*100, 'ko')
    plot(consistent.heta_grid_interp/econparams.frequency*100, consistent.heta2_consistent_interp/econparams.frequency*100, 'k-')
    plot(consistent.heta_optcons_duple(1)/econparams.frequency*100, consistent.heta_optcons_duple(2)/econparams.frequency*100, 'k.', 'MarkerSize', 30)
    hold off

    xlabel("$\eta_1$ (Annual percent)")
    ylabel("$\eta_2$ (Annual percent)")
    legend("Welfare Contours", "Calibrated BGP", "Optimal BGP", "Static", "Commitment", "Commitment (interp)", "Optimal commitment", "Consistent",...
        "Consistent (interp)", "Optimal consistent (interp)")

    offset = 0.25;

    % Place CEQ welfare differences to each policy
    text(econparams.heta/econparams.frequency*100 + offset, econparams.heta/econparams.frequency*100 + offset, num2str(round(nochange_ceq*100,2))+"\%", 'Interpreter', 'latex')
    text(heta_static_bgp/econparams.frequency*100 + offset, heta_static_bgp/econparams.frequency*100 + offset, num2str(round(bgp_ceq*100,2))+"\%", 'Interpreter', 'latex')
    text(heta_static_transition/econparams.frequency*100 + offset, heta_static_transition/econparams.frequency*100 + offset, num2str(round(static_ceq*100,2))+"\%", 'Interpreter', 'latex')
    text(commitment.heta_duple_commitment(1)/econparams.frequency*100 + offset, commitment.heta_duple_commitment(2)/econparams.frequency*100 + offset, num2str(round(duple_commitment_ceq*100,2))+"\%", 'Interpreter', 'latex')
    text(consistent.heta_optcons_duple(1)/econparams.frequency*100 + offset, consistent.heta_optcons_duple(2)/econparams.frequency*100 + offset, num2str(round(optcons_duple_ceq*100,2))+"\%", 'Interpreter', 'latex')

    for i = 1:length(consistent.heta_grid)
        text(consistent.heta_grid(i)/econparams.frequency*100 + offset, consistent.heta2_consistent(i)/econparams.frequency*100 + offset, num2str(round(consistent_duple_ceq(i)*100,2))+"\%", 'Interpreter', 'latex')
    end
    
end



