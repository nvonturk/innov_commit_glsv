%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_bgp_objects.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_bgp_objects(specification, econparams_bgp)
    
    f = figure("name", "bgp_objects");

    sgtitle('BGP $\mu$ and $x_{\sigma}$')

    subplot(1,2,1)
    plot(0:1:econparams_bgp.n, econparams_bgp.mu, 'k-')
    xlabel("$\sigma$")
    ylabel("$\mu_{\sigma}$ (Decimal share)")

    subplot(1,2,2)
    plot(-econparams_bgp.n:1:econparams_bgp.n, econparams_bgp.x_vector/econparams_bgp.frequency, 'k-')

    if specification == "calibration_entry"
        hold on
        plot(-econparams_bgp.n:1:0, econparams_bgp.x_vector_e/econparams_bgp.frequency, 'b--')
        hold off
        legend("Incumbent", "Entrant")
    end

    xlabel("$\sigma$")
    ylabel("$x_{\sigma}$ (Annual decimal)")

end



