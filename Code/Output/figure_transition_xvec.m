%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_transition_xvec.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_transition_xvec(econparams, econparams_bgp, consistent)
    
    heta_span = [1,3,5]/100*econparams.frequency;
    t_span = [1, 12.5, 24]/econparams.frequency;

    x_vec_mat = nan(length(t_span), length(heta_span), length(econparams_bgp.x_vector));
    optcons_duple_transition_anon = @(heta2) compute_transition(econparams, econparams_bgp, [consistent.heta_optcons_duple(1), heta2], repmat(econparams.rho,1,2), [econparams.interval_lengths(end-2), sum(econparams.interval_lengths(end-1:end))]);

    for i = 1:length(heta_span)

        transition_temp = optcons_duple_transition_anon(heta_span(i));

        for t = 1:length(t_span)

            x_vec_temp = transition_temp.x_vector_path(:,t_span(t) + 1);

            x_vec_mat(t,i,:) = x_vec_temp;

        end

    end

    f = figure('name', 'transition_xvec'); 
    for t = 1:length(t_span)

        subplot(1,length(t_span),t)

        hold on
        for i = 1:length(heta_span)
            plot(-econparams_bgp.n:1:econparams_bgp.n, squeeze(x_vec_mat(t,i,:))/econparams.frequency, ...
                'DisplayName', "$\eta_2 = " + num2str(round(heta_span(i)*100/econparams.frequency,2)) + "\%$")
        end
        hold off
        xlabel("$\sigma$")
        ylabel("Annual decimal")
        title(strcat("t = ", num2str(t_span(t)*econparams.frequency), " years"))
        legend show
    end

end



