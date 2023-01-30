%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_mu_transition.m
% Author: Nicholas von Turkovich
% Date: 06/21/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_mu_transition(econparams, econparams_bgp, consistent)
    
    transition_anon = @(heta_changes) compute_transition(econparams, econparams_bgp, heta_changes, repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);
    
    [~, heta1_span] = intersect(round(consistent.heta_grid*100/econparams.frequency,0), [1, 3, 6, 10]);
    txs = {};
    iter = 1;
    
    patterns = ["-", "--", "-.", "--x"];
    f = figure('name', 'mu_transition');
    
    for i = heta1_span'
        
        txs{iter} = transition_anon([consistent.heta_grid(i), consistent.heta2_consistent(i)]);
        
        hold on
        plot(0:econparams.n, txs{iter}.mu_path(:,econparams.interval_lengths(1) + 1), patterns(iter), 'MarkerSize', 6, ...
            'DisplayName', "$\eta = (" + num2str(round(consistent.heta_grid(i)*100/econparams.frequency,1)) + "\%, " ...
            + num2str(round(consistent.heta2_consistent(i)*100/econparams.frequency,1)) + "\%)$");
        hold off
        
        iter = iter + 1;
        
    end
    
    legend show
    legend boxoff
    xlim([0,20])
    xticks(0:2:20)
    grid off
    xlabel("Technology gap, $s$")
    ylabel("$\mu_s(T)$")        
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
            
end



