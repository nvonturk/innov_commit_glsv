%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_fb_quality_index.m
% Author: Nicholas von Turkovich
% Date: 01/02/2023
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_fb_quality_index(econparams, econparams_bgp, duple_commitment_transition, optcons_duple_transition, duple_commitment_terminal, optcons_duple_terminal)

    f = figure("name", "fb_quality_index");    

    % Compute the first best for the starting BGP (each firm on frontier 
    % gets same amount of R&D labor which maximizes innovation given convex
    % R&D costs)
    gfb_start = fb_growth_rate(econparams_bgp.L_research, econparams_bgp.mu(1), ...
        econparams.B, econparams.gamma, econparams.lambda);
    
    % Again for the terminal BGPs under the two alternative policies
    gfb_commit_terminal = fb_growth_rate(duple_commitment_terminal.L_research, ...
        duple_commitment_terminal.mu(1), ...
        econparams.B, econparams.gamma, econparams.lambda);
    gfb_optcons_terminal = fb_growth_rate(optcons_duple_terminal.L_research, ...
        optcons_duple_terminal.mu(1), ...
        econparams.B, econparams.gamma, econparams.lambda);

    % Again for the paths under the two alternative policies
    gfb_commit_transition = fb_growth_rate(duple_commitment_transition.L_research, ...
        duple_commitment_transition.mu_path(1,:), ...
        econparams.B, econparams.gamma, econparams.lambda);
    gfb_optcons_transition = fb_growth_rate(optcons_duple_transition.L_research, ...
        optcons_duple_transition.mu_path(1,:), ...
        econparams.B, econparams.gamma, econparams.lambda);

    % Add a burn in for this plot to show the state of the economy prior to
    % policy change
    nobs = 100;
    t_burnin = linspace(-10/econparams.frequency, 0, nobs);
    burnin = repmat(0, 1, nobs);

    % Create extended series with the burn in portion and transition
    t = [t_burnin(1:end-1), (0:length(optcons_duple_transition.heta_path))]*econparams.frequency;
    gfb_diff_commit = [burnin(1:end-1), (power(1+gfb_commit_transition,1/econparams.frequency)-1)*100 - ...
        (power(1+duple_commitment_transition.g,1/econparams.frequency)-1)*100];

    gfb_diff_optcons = [burnin(1:end-1), (power(1+gfb_optcons_transition,1/econparams.frequency)-1)*100 - ...
        (power(1+optcons_duple_transition.g,1/econparams.frequency)-1)*100];

    plot(t, gfb_diff_optcons, 'k-')
    hold on
    plot(t, gfb_diff_commit, 'g-')
    hold off
    xlabel("Years")
    ylabel("")
    leg = legend("$\eta^{TC^*}$", "$\eta^{C^*}$");
    set(leg, 'AutoUpdate', 'off')
    xline_prop = xline(0, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xline_prop = xline(econparams.interval_lengths(1)*econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
    xline_prop.LineWidth = 1;
    xlim([-10, econparams.interval_lengths(2)*econparams.frequency*3])
    grid off
    legend boxoff
    title("Difference from first-best quality index growth")
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end

function [gfb] = fb_growth_rate(L_research, mu_tied, B, gamma, lambda)

    rd_labor_demand = L_research./(1+mu_tied);
    innovation_rate = B.*rd_labor_demand.^gamma;
    gfb = log(lambda).*innovation_rate.*(1+mu_tied);

end



