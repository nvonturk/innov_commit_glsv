%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_bgp_objects_eta.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_bgp_objects_eta(specification, econparams, econparams_bgp)
    
    % Span of eta values to use
    heta_span = [2, 4, 6]/100*econparams.frequency;

    % Instantiate matrices to hold objects of interests
    mu_span = nan(length(econparams_bgp.mu), length(heta_span));
    x_span = nan(length(econparams_bgp.x_vector), length(heta_span));
    x_e_span = nan(length(econparams_bgp.x_vector_e), length(heta_span));

    % Colors and line types
    col = {'[1, 0, 0]', '[0, 0, 1]', '[0.4660 0.6740 0.1880]'};
    lstyle = {'--', '-', ':'};

    % Loop through values for eta and store objects
    for i = 1:length(heta_span)

        % Solve for BGP
        econparams_temp = update_heta(econparams, heta_span(i));
        econparams_temp_bgp = value_function_iteration_robust(econparams_temp, false);

        % Store mu and innovation strategies
        mu_span(:,i) = econparams_temp_bgp.mu;
        x_span(:,i) = econparams_temp_bgp.x_vector;

        if specification == "calibration_entry"
            x_e_span(:,i) = econparams_temp_bgp.x_vector_e;
        end

    end

    f = figure("name", "bgp_objects_eta");
    subplot(1,2,1)
    hold on
    for i = 1:length(heta_span)
        plot(0:1:econparams_bgp.n, mu_span(:,i), lstyle{i}, 'color' , col{i})
    end
    hold off

    grid off
    box on
    xlabel("Industry technology gap, $s$")
    title({"Distribution of technology gaps, $\mu_s$", ""})
    legend(cellstr(strcat('$\eta =', num2str(round(heta_span'*100/econparams.frequency,0)), '\%$')))
    legend boxoff

    subplot(1,2,2)
    hold on
    for i = 1:length(heta_span)
        if specification == "calibration"
            plot(-econparams_bgp.n:1:econparams_bgp.n, x_span(:,i)/econparams.frequency, lstyle{i}, 'color' , col{i})
        elseif specification == "calibration_entry"
            plot(-econparams_bgp.n:1:econparams_bgp.n, x_span(:,i)/econparams.frequency, "-", 'color' , col{i})
            plot(-econparams_bgp.n:1:0, x_e_span(:,i)/econparams.frequency, ":", 'color' , col{i})        
        end
    end
    xline(0, ':', 'color', [0.9, 0.9, 0.9])
    legend(cellstr(strcat('$\eta =', num2str(round(heta_span'*100/econparams.frequency,0)), '\%$')))
    legend boxoff
    hold off

    xlim([-10, 20]);

    grid off
    box on
    xlabel("Firm technology position, $\sigma$")
    title({"Cross section of firm innovation rates, $x_\sigma$", ""})


end



