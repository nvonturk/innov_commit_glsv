%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_sp2_best_response.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_sp2_best_response(econparams, contour_mat, consistent)

    % Choose a subset of eta 1 values
    [~, heta_span] = intersect(round(contour_mat.heta_grid*100/econparams.frequency,0), [1, 3, 6, 10]);

    heta_span = heta_span';
    
    iter = 1;

    f = figure("name", "sp2_best_response");

    patterns = ["-", "--", "-.", "--x"];
    hold on
    for i = heta_span

        
        leg = "$\eta_1 = " + num2str(round(contour_mat.heta_grid(i)*100/econparams.frequency,2)) + "\%$";

        % Plot on the x axis the grid of eta 2 values and on the y axis the
        % optimal eta 2 value from SP2's perspective
        plot(contour_mat.heta_grid*100/econparams.frequency, contour_mat.sp2_static(i,:)*100/econparams.frequency, ...
            patterns(iter), 'DisplayName', leg, 'MarkerSize', 6);

        iter = iter + 1;

    end
    
    for i = heta_span
        
        % Pick out the eta 1 value
        heta1 = contour_mat.heta_grid(i);
        
        % For that same eta 1 find the consistent eta2 and plot that
        % (should intersect with one of the curves plotted)
        heta2_cons_idx = find(consistent.heta_grid == heta1);
        plot(consistent.heta2_consistent(heta2_cons_idx)*100/econparams.frequency, consistent.heta2_consistent(heta2_cons_idx)*100/econparams.frequency, 'ko',...
            'HandleVisibility', 'off')
        
    end

    % Plot 45 degree line
    plot(0:(max(contour_mat.heta_grid)/econparams.frequency*100), 0:(max(contour_mat.heta_grid)/econparams.frequency*100), '-', ...
        'LineWidth', 1, 'DisplayName', "45 degree line", 'color', [0.7, 0.7, 0.7]);

    hold off
    legend show
    grid off
    legend boxoff
    leg = get(gca, "Legend");
    leg.Location = 'northwest';
    xlabel("Long-run expiry rate, $\eta_2$")
    ylabel("$\tilde{\eta}_2(X_1(\eta_1,\eta_2),\eta_1)$")
    xlim([0, 16])
    ylim([0, 16])
    xticks(0:4:16)
    yticks(0:4:16)
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")

end



