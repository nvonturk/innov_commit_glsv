%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_adjusting_t_triple.m
% Author: Nicholas von Turkovich
% Date: 06/20/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_adjusting_t_triple(filestart_data, econparams)

    % Adjustments
    adjustments = [10, 15, 20, 25, 30, 35, 50];

    % Load necessary data for each plot
    adj_econparams = {};
    adj_econparams_bgp = {};
    adj_commitment = {};
    adj_consistent = {};
    adj_commitment_txs = {};
    adj_consistent_txs = {};
    
    adj_optcons = nan(3,size(adjustments,2));
    adj_optcommit = nan(3,size(adjustments,2));
        
    adj_optcons_meanmarkup_2T = nan(1,size(adjustments,2));
    adj_optcommit_meanmarkup_2T = nan(1,size(adjustments,2));
    
    adj_optcons_prodg_2T = nan(1,size(adjustments,2));
    adj_optcommit_prodg_2T = nan(1,size(adjustments,2));
    
    adj_optcons_meanmarkup_3T = nan(1,size(adjustments,2));
    adj_optcommit_meanmarkup_3T = nan(1,size(adjustments,2));
    
    adj_optcons_prodg_3T = nan(1,size(adjustments,2));
    adj_optcommit_prodg_3T = nan(1,size(adjustments,2));

    for i = 1:length(adjustments)
        
        % Load the data structures associated with the particular T
        adj_econparams{i} = load(filestart_data + "_T_series_alt" + adjustments(i) + "_hub60_econparams.mat");
        adj_econparams_bgp{i} = load(filestart_data + "_T_series_alt" + adjustments(i) + "_hub60_econparams_bgp.mat");
        adj_commitment{i} = load(filestart_data + "_T_series_alt" + adjustments(i) + "_hub60_heta_triple_commitment.mat");
        adj_consistent{i} = load(filestart_data + "_T_series_alt" + adjustments(i) + "_hub60_tripcons.mat");

        % Store the optimal, consistent duple
        adj_optcons(:,i) = adj_consistent{i}.heta_max_optcons';
        
        % Store the optimal, commitment duple
        adj_optcommit(:,i) = adj_commitment{i}.heta_triple_commitment';
        
        % Compute the transitions associated with the optimal commitment
        % and optimal consistent policies
        adj_commitment_txs{i} = compute_transition(adj_econparams{i}, adj_econparams_bgp{i}, adj_commitment{i}.heta_triple_commitment, repmat(adj_econparams{i}.rho,1,3), adj_econparams{i}.interval_lengths);
        adj_consistent_txs{i} = compute_transition(adj_econparams{i}, adj_econparams_bgp{i}, adj_consistent{i}.heta_max_optcons, repmat(adj_econparams{i}.rho,1,3), adj_econparams{i}.interval_lengths);

        % Store the mean gross markup at time 2T for each transition
        adj_optcommit_meanmarkup_2T(i) = adj_commitment_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*2 + 1);
        adj_optcons_meanmarkup_2T(i) = adj_consistent_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*2 + 1);
        
        % Store the TFP growth at time 2T for each transition
        adj_optcommit_prodg_2T(i) = 100*(power(adj_commitment_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*2 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
        adj_optcons_prodg_2T(i) = 100*(power(adj_consistent_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*2 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
        
        % Store the mean gross markup at time 3T for each transition
        adj_optcommit_meanmarkup_3T(i) = adj_commitment_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*3 + 1);
        adj_optcons_meanmarkup_3T(i) = adj_consistent_txs{i}.mean_gross_markup(adj_econparams{i}.interval_lengths(1)*3 + 1);
        
        % Store the TFP growth at time 3T for each transition
        adj_optcommit_prodg_3T(i) = 100*(power(adj_commitment_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*3 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
        adj_optcons_prodg_3T(i) = 100*(power(adj_consistent_txs{i}.prodg(adj_econparams{i}.interval_lengths(1)*3 + 1) + 1, 1/adj_econparams{i}.frequency) - 1);
    
    end

    f1 = figure('name', 'adjusting_t_triple_2T');
    
    subplot(1,3,1)
    
    hold on
    plot(adjustments, adj_optcons(2,:)*100/econparams.frequency, 'k--')
    plot(adjustments, adj_optcommit(2,:)*100/econparams.frequency, 'g--')
    plot(adjustments, adj_optcons(3,:)*100/econparams.frequency, 'k-')
    plot(adjustments, adj_optcommit(3,:)*100/econparams.frequency, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Medium- and long-run expiry rates")
    xlim([15,50])
    ylim([0,10])
    
    grid off
    legend("Time-consistent $\eta_2$", "Commitment $\eta_2$", "Time-consistent $\eta_3$", "Commitment $\eta_3$")
    legend boxoff
    
    subplot(1,3,2)
    
    hold on
    plot(adjustments, adj_optcons_prodg_2T, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, adj_optcommit_prodg_2T, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Productivity growth at $t=2T$")
    xlim([15,50])
    
    grid off
    

    subplot(1,3,3)
    
    hold on
    plot(adjustments, (adj_optcons_meanmarkup_2T - 1)*100, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, (adj_optcommit_meanmarkup_2T - 1)*100, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Average markup at $t=2T$")
    xlim([15,50])
    
    grid off
    
    f2 = figure('name', 'adjusting_t_triple_3T');
    
    subplot(1,3,1)
    
    hold on
    plot(adjustments, adj_optcons(2,:)*100/econparams.frequency, 'k--')
    plot(adjustments, adj_optcommit(2,:)*100/econparams.frequency, 'g--')
    plot(adjustments, adj_optcons(3,:)*100/econparams.frequency, 'k-')
    plot(adjustments, adj_optcommit(3,:)*100/econparams.frequency, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Medium- and long-run expiry rates")
    xlim([15,50])
    ylim([0,10])
    
    grid off
    legend("Time-consistent $\eta_2$", "Commitment $\eta_2$", "Time-consistent $\eta_3$", "Commitment $\eta_3$")
    legend boxoff
    
    subplot(1,3,2)
    
    hold on
    plot(adjustments, adj_optcons_prodg_3T, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, adj_optcommit_prodg_3T, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Productivity growth at $t=3T$")
    xlim([15,50])
    
    grid off
    

    subplot(1,3,3)
    
    
    hold on
    plot(adjustments, (adj_optcons_meanmarkup_3T - 1)*100, 'k-', 'MarkerFaceColor', 'black')
    plot(adjustments, (adj_optcommit_meanmarkup_3T - 1)*100, 'g-')
    hold off
    
    xlabel("T (years)")
    ylabel("Percent")
    title("Average markup at $t=3T$")
    xlim([15,50])
    
    grid off
    
    figlist = [f1, f2];
    
    for i = 1:length(figlist)
        
        f = figlist(i);
        
        f.Units = 'centimeters';
        orient(f, 'landscape')
        f.PaperPosition = [-1, -1, 12, 4];
        saveas(f, "../../Output/Figures/Main/Figure_" + f.Name + ".png")
        
    end

end



