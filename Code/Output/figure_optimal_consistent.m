%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_optimal_consistent.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_optimal_consistent(econparams, econparams_bgp, optcons_duple_transition)

    f = figure("name", "optimal_consistent");

    % Add a burn in for this plot to show the state of the economy prior to
    % policy change
    nobs = 100;
    t_burnin = linspace(-10/econparams.frequency, 0, nobs);
    g_burnin = repmat(econparams_bgp.g, 1, nobs);
    mean_gross_markup_burnin = repmat(econparams_bgp.mean_gross_markup, 1, nobs);
    pct75_gross_markup_burnin = repmat(econparams_bgp.pct75_gross_markup, 1, nobs);
    pct90_gross_markup_burnin = repmat(econparams_bgp.pct90_gross_markup, 1, nobs);

    % Create extended series with the burn in portion and transition
    t = [t_burnin(1:end-1), (0:length(optcons_duple_transition.heta_path))]*econparams.frequency;
    g = (power(1 + [g_burnin(1:end-1), optcons_duple_transition.g], 1/econparams.frequency)-1)*100;
    prodg = (power(1 + [g_burnin(1:end-1), optcons_duple_transition.prodg], 1/econparams.frequency)-1)*100;
    mean_gross_markup = ([mean_gross_markup_burnin(1:end-1), optcons_duple_transition.mean_gross_markup]-1)*100;
    pct75_gross_markup = ([pct75_gross_markup_burnin(1:end-1), optcons_duple_transition.pct75_gross_markup]-1)*100;
    pct90_gross_markup = ([pct90_gross_markup_burnin(1:end-1), optcons_duple_transition.pct90_gross_markup]-1)*100;

    subplot(1,2,1)
    plot(t, g, 'k-')
    hold on
    plot(t, prodg, 'r--')
    hold off
    xlabel("Years")
    ylabel("")
    leg = legend("$g^Q$", "$g^{TFP}$");
    set(leg, 'AutoUpdate', 'off')
    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xline_prop = xline(econparams.interval_lengths(1)*econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-10, econparams.interval_lengths(2)*econparams.frequency*3])
    grid off
    legend boxoff
    title("Growth")

    subplot(1,2,2)
    plot(t, mean_gross_markup, 'k-')
    hold on
%     plot(t, pct75_gross_markup, 'r--')
    plot(t, pct90_gross_markup, 'r--')
    hold off
    xlabel("Years")
    ylabel("Percent")
    leg = legend("Mean markup", "$90^{th}$ percentile markup");
    set(leg, 'AutoUpdate', 'off')
    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xline_prop = xline(econparams.interval_lengths(1)*econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-10, econparams.interval_lengths(2)*econparams.frequency*3])
    grid off
    legend boxoff
    title("Net markup")
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



