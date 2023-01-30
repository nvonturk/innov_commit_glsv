%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_menu_xmu.m
% Author: Nicholas von Turkovich
% Date: 06/6/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1] = figure_menu_xmu(econparams, duple_commitment_transition, ...
    optcons_duple_transition)
    
    f = figure('name', 'menu_xmu');
    
    subplot(1,2,1)
    plot(-econparams.n:1:econparams.n, duple_commitment_transition.x_vector_path(:,econparams.interval_lengths(1)*2 + 1)/econparams.frequency, 'g-')
    hold on
    plot(-econparams.n:1:econparams.n, optcons_duple_transition.x_vector_path(:,econparams.interval_lengths(1)*2 + 1)/econparams.frequency, 'k-')
    hold off
    xlabel("Technology position, $\sigma$")
    ylabel("$x_{\sigma}(2T)$")
    ylim([0, 0.22])
    xlim([-15,15])
    title("Cross section of firm innovation rates")
    grid off
    legend("$\eta^{C^*}$", "$\eta^{TC^*}$")
    legend boxoff
    
    subplot(1,2,2)
    plot(0:econparams.n, duple_commitment_transition.mu_path(:,econparams.interval_lengths(1)*2 + 1), 'g-')
    hold on
    plot(0:econparams.n, optcons_duple_transition.mu_path(:,econparams.interval_lengths(1)*2 + 1), 'k-')
    hold off
    xlabel("Technology gap, s")
    ylabel("$\mu_{s}(2T)$")
    title("Distribution of technology gaps")
    grid off
    ylim([0, 0.22])
    xlim([0,20])
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
            
end



