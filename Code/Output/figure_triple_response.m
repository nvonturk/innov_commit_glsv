%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_triple_response.m
% Author: Nicholas von Turkovich
% Date: 06/15/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_triple_response(specification)

    responses = readtable("Output/Logs/Main/" + specification + "_tripconslog_compiled_long.csv");

    heta1_vals = [0, 6, 12, 18];

    f = figure("name", strcat("triple_response_", specification));
    set(f, 'DefaultTextFontSize', 10);
    set(f, 'DefaultAxesFontSize', 10);
    
    for i = 1:length(heta1_vals)
        
        responses_filt = responses(responses.heta1_intended == heta1_vals(i),:);
        F = scatteredInterpolant(responses_filt.heta2_intended, responses_filt.heta3_intended, responses_filt.loss, 'linear', 'none');

        Xq = linspace(min(responses_filt.heta2_intended), max(responses_filt.heta2_intended), 100);
        Yq = linspace(min(responses_filt.heta3_intended), max(responses_filt.heta3_intended), 100);
        
        [heta2_interp_grid, heta3_interp_grid] = ndgrid(Xq, Yq);
        loss_interp = F(heta2_interp_grid, heta3_interp_grid);

        [~, min_idx] = min(responses_filt.loss);
        
        subplot(2, 2, i)
        surfc(heta2_interp_grid, heta3_interp_grid, loss_interp, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        hold off
        xlim([0,20])
        ylim([0,20])
        xlabel(["$\eta_2$"])
        ylabel(["$\eta_3$"])
        zlabel("$|\eta_2 - \tilde{\eta}_2| + |\eta_3 - \tilde{\eta}_3|$")
        title([strcat("Initial expiry rate, ", "$\eta_1 = ", sprintf("%1.0f", heta1_vals(i)), "\%$")])
        hold on
        plot3(responses_filt.heta2_intended(min_idx), responses_filt.heta3_intended(min_idx), responses_filt.loss(min_idx), 'k.', 'MarkerSize', 30);
        hold off 
        view([-45, 30])
        
        
    end
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 6, 6];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



