%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_transition_objects.m
% Author: Nicholas von Turkovich
% Date: 06/1/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_transition_objects(econparams, econparams_bgp, commitment)


    duple_commitment_transition_anon = @(heta2) compute_transition(econparams, econparams_bgp, [commitment.heta_duple_commitment(1), heta2], repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

    f = figure("name", "transition_objects");
    
    t_val = (econparams.interval_lengths(1)*econparams.frequency - 1)/econparams.frequency + 1;
    heta2_vals = [1.5, 3]*econparams.frequency/100;
    
    subplot(1,2,1)
    
    for idx = 1:length(heta2_vals)
        
        hold on
        
        tx_temp = duple_commitment_transition_anon(heta2_vals(idx));
        plot(-econparams.n:1:econparams.n, tx_temp.x_vector_path(:,t_val)/econparams.frequency);
        
        hold off
        
    end
    
    xlabel("Firm technology position, $\sigma$")
    title(["Firm innovation strategy,"; "$x_{\sigma}(T-1)$ given $(\eta_1^{\{C,*\}}, \eta_2)$"])
    grid off
    leg = legend(strcat("$\eta_2 = ", string(heta2_vals*100/econparams.frequency), "\%$"), 'Location', 'SouthEast');
    legend boxoff
    set(leg, 'AutoUpdate', 'off')
    xlim([-15, 15])
    
    xline(0, '--', 'color', [0.7, 0.7, 0.7])
    
    subplot(1,2,2)
    
    for idx = 1:length(heta2_vals)
        
        hold on
        
        tx_temp = duple_commitment_transition_anon(heta2_vals(idx));
        plot(0:1:econparams.n, tx_temp.mu_path(:,t_val));
        
        hold off
        
    end
    
    xlabel("Industry technology gap, $s$")
    title(["Distribution of technology gaps,"; "$\mu_s(T-1)$ given $(\eta_1^{\{C,*\}}, \eta_2)$"])
    grid off
    xlim([0,15])
    

end



