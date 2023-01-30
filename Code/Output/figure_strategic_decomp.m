%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_strategic_decomp.m
% Author: Nicholas von Turkovich
% Date: 06/27/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = figure_strategic_decomp(strat_decomp)

    figs = {};

    f = figure("name", "strategic_decomp");

    time_span = (0:length(strat_decomp.dg_de_t)-1)*strat_decomp.econparams.frequency;
    
    subplot(2,1,1)
        
    plot(time_span, strat_decomp.dg_de_t, 'k-')
    hold on
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.composition_effect_t, 'r--')
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.intensive_margin_t, 'b-.');
    hold off
    
    leg = legend('$\frac{\partial g^Q}{\partial \eta_2}$', 'Composition effect', 'Intensive margin', 'Location', 'SouthWest');
    set(leg, 'AutoUpdate', 'off')
    set(leg, 'FontSize', 8)
    
    hold on
    plot(time_span(end), strat_decomp.dg_de, 'ko')
    plot(time_span(end), strat_decomp.composition_effect*log(strat_decomp.econparams.lambda), 'ro')
    plot(time_span(end), strat_decomp.intensive_margin*log(strat_decomp.econparams.lambda), 'bo')
    hold off
    
    interval_lengths_cum = cumsum(strat_decomp.interval_lengths);
    for i = 1:length(interval_lengths_cum)-1
        xline_prop = xline(interval_lengths_cum(i)*strat_decomp.econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
        xline_prop.LineWidth = 1;
    end
        
    xlabel("Years")
    xlim([0,75])
    ylim([-0.175, 0.05])
    grid off
    legend boxoff
    
    subplot(2,1,2)
    
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.intensive_margin_t, 'b-')
    hold on
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.kw_effect_t, 'color', [255, 165, 0]./256)
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.strategic_effect_escape_t, 'color', [50, 205, 50]./256)
    plot(time_span, log(strat_decomp.econparams.lambda)*strat_decomp.strategic_effect_trickle_t, 'c-')
    hold off
    
    leg = legend('Intensive margin', 'Own-innovation effect', 'Escape competition', 'Trickle down', 'Location', 'SouthWest');
    set(leg, 'AutoUpdate', 'off')
    set(leg, 'FontSize', 8)
    
    hold on
    plot(time_span(end), strat_decomp.kw_effect*log(strat_decomp.econparams.lambda), 'o', 'color', [255, 165, 0]./256)
    plot(time_span(end), strat_decomp.escape_effect*log(strat_decomp.econparams.lambda), 'o', 'color', [50, 205, 50]./256)
    plot(time_span(end), strat_decomp.trickle_effect*log(strat_decomp.econparams.lambda), 'co')
    plot(time_span(end), strat_decomp.intensive_margin*log(strat_decomp.econparams.lambda), 'bo')
    hold off
    
    for i = 1:length(interval_lengths_cum)-1
        xline_prop = xline(interval_lengths_cum(i)*strat_decomp.econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
        xline_prop.LineWidth = 1;
    end
    
    xlabel("Years")
    xlim([0,75])
    ylim([-0.175, 0.05])
    grid off
    legend boxoff
    
    figs{1} = f;
    
    f = figure('name', 'strategic_decomp_cross');
    
    sigma = 4;
    t_idx = strat_decomp.econparams.interval_lengths(1) + 1 - 1/strat_decomp.econparams.frequency;
    dxs_dxsc = zeros(1, 2*strat_decomp.econparams.n + 1);
    for j = 1:2*strat_decomp.econparams.n
       dxs_dxsc(j) = strat_decomp.dxs_dxs_t_unweighted{j}(strat_decomp.econparams.n + 1 + sigma, t_idx);
    end
    
    plot(-strat_decomp.econparams.n:1:strat_decomp.econparams.n, dxs_dxsc)
    xlabel("Technology position, $\sigma'$")
    title("$\sum_{t' > T - 1}\frac{\partial x_" + num2str(sigma) + "(T-1)}{\partial x_{\sigma'_c}(t')}$")
    grid off
    xlim([-10, 10])
    xline_prop = xline(-sigma, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    
    figs{2} = f;
        
    for i = 1:length(figs)
        f = figs{i};
        f.Units = 'centimeters';
        orient(f, 'landscape')
        f.PaperPosition = [-1, -1, 8, 4];
        saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    end

end



