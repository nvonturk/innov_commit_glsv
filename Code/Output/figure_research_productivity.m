%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_research_productivity.m
% Author: Nicholas von Turkovich
% Date: 01/02/2023
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_research_productivity(econparams, econparams_bgp, duple_commitment_transition, optcons_duple_transition)

    f = figure("name", "research_productivity");

    % Compute starting BGP research productivity
    scale_research_productivity = econparams_bgp.g/econparams_bgp.L_research;

    % Add a burn in for this plot to show the state of the economy prior to
    % policy change
    nobs = 100;
    t_burnin = linspace(-10/econparams.frequency, 0, nobs);
    research_productivity_burnin = repmat(1, 1, nobs);
    prodg_burnin = repmat((power(1+econparams_bgp.g, 1/econparams.frequency)-1)*100, 1, nobs);

    % Create extended series with the burn in portion and transition
    t = [t_burnin(1:end-1), (0:length(optcons_duple_transition.heta_path))]*econparams.frequency;
    research_productivity_commit = [research_productivity_burnin(1:end-1), duple_commitment_transition.research_productivity./scale_research_productivity];
    research_productivity_optcons = [research_productivity_burnin(1:end-1), optcons_duple_transition.research_productivity./scale_research_productivity];
    prodg_commit = [prodg_burnin(1:end-1),(power(1+duple_commitment_transition.prodg,1/econparams.frequency)-1)*100];
    prodg_optcons = [prodg_burnin(1:end-1),(power(1+optcons_duple_transition.prodg,1/econparams.frequency)-1)*100];

    subplot(1,2,1)
    plot(t, prodg_commit - prodg_optcons, 'k-')
    xlabel("Years")
    ylabel("Percentage points")
    leg = legend("$\eta^{C^*} - \eta^{TC^*}$");
    set(leg, 'AutoUpdate', 'off')
    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xline_prop = xline(econparams.interval_lengths(1)*econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-10, econparams.interval_lengths(2)*econparams.frequency*3])
    grid off
    legend boxoff
    title("TFP growth")
    
    subplot(1,2,2)
    plot(t, research_productivity_commit - research_productivity_optcons, 'k-')
    xlabel("Years")
%     ylabel("Aggregate research productivity")
    leg = legend("$\eta^{C^*} - \eta^{TC^*}$");
    set(leg, 'AutoUpdate', 'off')
    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xline_prop = xline(econparams.interval_lengths(1)*econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-10, econparams.interval_lengths(2)*econparams.frequency*3])
    grid off
    legend boxoff
    title("Aggregate efficiency of R\&D")

    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



