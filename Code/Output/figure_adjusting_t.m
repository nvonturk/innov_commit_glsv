%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_adjusting_t.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_adjusting_t(filestart_data, econparams)

    % Adjustments
    adjustments = [10, 15, 20, 25, 30, 35, 50];

    % Load necessary data for each plot
    adj_econparams = {};
    adj_econparams_bgp = {};
    adj_commitment = {};
    adj_consistent = {};
    adj_commitment_txs = {};
    adj_consistent_txs = {};
    
    adj_optcons = nan(2,size(adjustments,2));
    adj_optcommit = nan(2,size(adjustments,2));
    adj_optcommit_interp = nan(2,size(adjustments,2));
        
    adj_optcons_meanmarkup = nan(1,size(adjustments,2));
    adj_optcommit_meanmarkup = nan(1,size(adjustments,2));
    
    adj_optcons_prodg = nan(1,size(adjustments,2));
    adj_optcommit_prodg = nan(1,size(adjustments,2));

    for i = 1:length(adjustments)
        
        % Load the data structures associated with the particular T
        adj_econparams{i} = load(filestart_data + "_T_series" + adjustments(i) + "_hub40_econparams.mat");
        adj_econparams_bgp{i} = load(filestart_data + "_T_series" + adjustments(i) + "_hub40_econparams_bgp.mat");
        adj_commitment{i} = load(filestart_data + "_T_series" + adjustments(i) + "_hub40_commitment.mat");
        adj_consistent{i} = load(filestart_data + "_T_series" + adjustments(i) + "_hub40_consistent.mat");

        % Store the optimal, consistent duple
        adj_optcons(:,i) = adj_consistent{i}.heta_optcons_duple';
        
        % Store the optimal, commitment duple
        adj_optcommit(:,i) = adj_commitment{i}.heta_duple_commitment';
        
        % Store the optimal, commitment duple (from interpolation)
        adj_optcommit_interp(:,i) = adj_commitment{i}.heta_opt_duple';
        
        % Compute the transitions associated with the optimal commitment
        % and optimal consistent policies
        adj_commitment_txs{i} = compute_transition(adj_econparams{i}, adj_econparams_bgp{i}, adj_commitment{i}.heta_duple_commitment, repmat(adj_econparams{i}.rho,1,2), [adj_econparams{i}.interval_lengths(end-2), sum(adj_econparams{i}.interval_lengths(end-1:end))]);
        adj_consistent_txs{i} = compute_transition(adj_econparams{i}, adj_econparams_bgp{i}, adj_consistent{i}.heta_optcons_duple, repmat(adj_econparams{i}.rho,1,2), [adj_econparams{i}.interval_lengths(end-2), sum(adj_econparams{i}.interval_lengths(end-1:end))]);

        % Store the mean gross markup at time 2T for each transition
        adj_optcommit_meanmarkup(i) = adj_commitment_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*2 + 1);
        adj_optcons_meanmarkup(i) = adj_consistent_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*2 + 1);
        
        % Store the TFP growth at time 2T for each transition
        adj_optcommit_prodg(i) = 100*(power(adj_commitment_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*2 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
        adj_optcons_prodg(i) = 100*(power(adj_consistent_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*2 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
    
    end

    f1 = figure('name', 'adjusting_t');
    
    subplot(1,3,1)
    
    hold on
    plot(adjustments, adj_optcons(2,:)*100/econparams.frequency, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, adj_optcommit(2,:)*100/econparams.frequency, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Long-run expiry rate, $\eta_2$")
    xlim([15,50])
    ylim([0,4.5])
    
    grid off
    legend("Time-consistent", "Commitment")
    legend boxoff
    
    subplot(1,3,2)
    
    hold on
    plot(adjustments, adj_optcons_prodg, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, adj_optcommit_prodg, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Productivity growth at $t=2T$")
    xlim([15,50])
    
    grid off
    

    subplot(1,3,3)
    
    hold on
    plot(adjustments, (adj_optcons_meanmarkup - 1)*100, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, (adj_optcommit_meanmarkup - 1)*100, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Average markup at $t=2T$")
    xlim([15,50])
    
    grid off
    
    figlist = [f1];
    
    for i = 1:length(figlist)
        
        f = figlist(i);
        
        f.Units = 'centimeters';
        orient(f, 'landscape')
        f.PaperPosition = [-1, -1, 12, 4];
        saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
        
    end

end



