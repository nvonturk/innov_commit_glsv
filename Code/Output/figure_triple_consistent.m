%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_triple_consistent.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_triple_consistent(econparams, consistent_triple_interp, heta_triple_commitment)

    f = figure("name", "triple_consistent");

    subplot(2,1,1)
    plot(consistent_triple_interp.heta1_interp/econparams.frequency*100, consistent_triple_interp.heta2_consistent_interp/econparams.frequency*100, 'k-')
    hold on
    plot(consistent_triple_interp.heta_max_optcons(1)/econparams.frequency*100, consistent_triple_interp.heta_max_optcons(2)/econparams.frequency*100, 'k.', 'MarkerSize', 35);
    plot(heta_triple_commitment(1)/econparams.frequency*100, heta_triple_commitment(2)/econparams.frequency*100, 'g.', 'MarkerSize', 35)
    hold off
    xlabel("Initial patent expiry rate, $\eta_1$")
    ylabel(["Medium-run"; "expiry rate, $\eta_2$"])
    legend("Time-consistent duples", "Optimal time-consistent policy", "Optimal policy with commitment")
    ylim([0,20])
    grid off
    legend boxoff

    subplot(2,1,2)
    plot(consistent_triple_interp.heta1_interp/econparams.frequency*100, consistent_triple_interp.heta3_consistent_interp/econparams.frequency*100, 'k-')
    hold on
    plot(consistent_triple_interp.heta_max_optcons(1)/econparams.frequency*100, consistent_triple_interp.heta_max_optcons(3)/econparams.frequency*100, 'k.', 'MarkerSize', 35);
    plot(heta_triple_commitment(1)/econparams.frequency*100, heta_triple_commitment(3)/econparams.frequency*100, 'g.', 'MarkerSize', 35)
    hold off
    xlabel("Initial patent expiry rate, $\eta_1$")
    ylabel(["Long-run"; "expiry rate, $\eta_3$"])
    ylim([0,5])
    grid off
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 7, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



