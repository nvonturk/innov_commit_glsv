%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_bgp_moments_eta.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_bgp_moments_eta(econparams, econparams_bgp, consistent)
    
    % Span of eta values to use
    heta_span = linspace(0, 0.05*econparams.frequency, 50);
    heta_span = [heta_span(1:end-1), linspace(0.05*econparams.frequency, 0.2*econparams.frequency, 10)];

    % Create empty data structures for growth and mean net markups
    g_span = nan(size(heta_span));
    mean_net_markup_span = nan(size(heta_span));

    parfor i = 1:length(heta_span)

        econparams_temp = update_heta(econparams, heta_span(i));
        econparams_temp_bgp = value_function_iteration_robust(econparams_temp, false);

        g_span(i) = (power(econparams_temp_bgp.g+1, 1/econparams.frequency)-1)*100;
        mean_net_markup_span(i) = (econparams_temp_bgp.mean_gross_markup - 1)*100;

    end

    f = figure("name", "bgp_moments_eta");

    sgtitle('Growth and mean markups')

    subplot(2,1,1)
    plot(heta_span*100/econparams.frequency, g_span, 'k-')
    hold on
    xline(econparams_bgp.heta*100/econparams.frequency, 'r--')
    xline(consistent.heta_optcons_duple(2)*100/econparams.frequency, 'r:')
    hold off
    legend('$g(\eta)$', '$\eta^{calibrated}$', '$\eta_2^{consistent}$')
    xlabel("$\eta$ (Annual percent)")
    ylabel("$g^{BGP}$ (Annual percent)")
    xlim([0,20])

    subplot(2,1,2)
    plot(heta_span*100/econparams.frequency, mean_net_markup_span, 'k-')
    hold on
    xline(econparams_bgp.heta*100/econparams.frequency, 'r--')
    xline(consistent.heta_optcons_duple(2)*100/econparams.frequency, 'r:')
    hold off
    legend('$meanmarkups(\eta)$', '$\eta_{calibrated}$', '$\eta_2^{consistent}$')
    xlabel("$\eta$ (Annual percent)")
    ylabel("Mean markup (Net percent)")
    xlim([0,20])


end



