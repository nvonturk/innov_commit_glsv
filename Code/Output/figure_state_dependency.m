%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_state_dependency.m
% Author: Nicholas von Turkovich
% Date: 06/29/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = figure_state_dependency(econparams, commitment_statedependent, consistent_statedependent)
    
    [heta_comp, heta_uncomp] = ndgrid(consistent_statedependent.heta_grid_interp, consistent_statedependent.heta_grid_interp);
    
    f = figure('name', 'statedep');
    surfc(heta_comp*100/econparams.frequency, heta_uncomp*100/econparams.frequency, ...
        consistent_statedependent.heta2_consistent_interp_statedep*100/econparams.frequency, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    
    hold on
    plot3(consistent_statedependent.heta_optcons_duple_statedep(1)*100/econparams.frequency,...
        consistent_statedependent.heta_optcons_duple_statedep(2)*100/econparams.frequency,...
        consistent_statedependent.heta_optcons_duple_statedep(3)*100/econparams.frequency, 'k.', 'MarkerSize', 30)
    hold off
    
    xlabel("$\eta_1^c$")
    ylabel("$\eta_1^u$")
    zlabel("$\eta_2$")
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    
end



