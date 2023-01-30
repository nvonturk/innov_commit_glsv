%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_spline_approx.m
% Author: Nicholas von Turkovich
% Date: 06/15/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_spline_approx(econparams, consistent)

    
    f = figure("name", "spline_approx");
    
    plot(consistent.heta_grid*100/econparams.frequency, consistent.heta2_consistent*100/econparams.frequency, 'ko')
    hold on
    plot(consistent.heta_grid_interp*100/econparams.frequency, consistent.heta2_consistent_interp*100/econparams.frequency, 'k-')
    plot(consistent.heta_grid_interp_check*100/econparams.frequency, consistent.heta2_consistent_interp_check*100/econparams.frequency, 'ro') 
    hold off    
  
%     legend("Basis of interpolation", "Interpolated values", "Time-consistent duples solved directly")
%     legend boxoff
    
    xlabel("Initial patent expiry rate ($\eta_1$)")
    ylabel("Long-term patent expiry rate ($\eta_2$)")
    xlim([0, 20])
    ylim([0, 10])
    
    grid off

    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



