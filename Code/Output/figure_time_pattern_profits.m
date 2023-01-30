%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_optimal_consistent.m
% Author: Nicholas von Turkovich
% Date: 06/13/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_time_pattern_profits(sims)

    f = figure("name", "time_pattern_profits");

    % Starting position of firms that experience an innovation
    sigmas = [-2,0,2];

    % For each simulation (i.e., t = T - 24, T - 1)
    for sim_idx = 1:length(sims)

        sim = sims{sim_idx};
        sim.sigmas = sigmas;
        sim.firm_profits_diff = nan(length(sigmas), size(sim.firm_profits,2));
        sim.firm_profits_diff_disc = nan(length(sigmas), size(sim.firm_profits,2));
        sim.firm_profits_diff_norm = nan(length(sigmas), size(sim.firm_profits,2));
        sim.val_diff = nan(length(sigmas), 1);

        for sigma_idx = 1:length(sigmas)

            % Starting position sigma
            sigma = sigmas(sigma_idx);

            % Get prob. distribution over new states
            pdist = sim.econparams.innovation_transitions(:,sim.econparams.n + 1 + sigma);

            % Find possible next sigmas
            next_sigmas = find(pdist > 0) - sim.econparams.n - 1;
            next_sigmas_weights = pdist(pdist > 0);
            
            % Determine the index of the next possible states in the
            % results for the simulation
            [~, ~, next_sigmas_idx] = intersect(next_sigmas, sim.sigma_start);
            
            % Determine the index of the current state in the simulation
            start_idx = find(sim.sigma_start == sigma);

            % Compute the expected difference in profits
            sim.firm_profits_diff(sigma_idx,:) = next_sigmas_weights' * sim.firm_profits(next_sigmas_idx,:) - sim.firm_profits(start_idx, :);
            
            % Discounted profit difference
            sim.firm_profits_diff_disc(sigma_idx,:) = sim.firm_profits_diff(sigma_idx,:) .* exp(-(1:length(sim.firm_profits_diff(sigma_idx,:)))*sim.econparams.rho);

            % Compute value function differences at that point in time
            sim.val_diff(sigma_idx) = next_sigmas_weights' * sim.transition.value_function_path(sim.econparams.n + 1 + next_sigmas,sim.t) - ...
                sim.transition.value_function_path(sim.econparams.n + 1 + sigma,sim.t);

            % Compute norm series (nansum gives close to 1 for each
            % position)
            sim.firm_profits_diff_norm(sigma_idx,:) = sim.firm_profits_diff_disc(sigma_idx,:)/sim.val_diff(sigma_idx);

        end

        sims{sim_idx} = sim;

    end

    patterns = ["-", "--", "-."];

    for sim_idx = 1:length(sims)

        burnin = zeros(1,sims{sim_idx}.t-1);

        subplot(1,length(sims),sim_idx)

        hold on

        for sigma_idx = 1:length(sigmas)

            plot((0:length(sims{sim_idx}.transition.heta_path)-1)*sims{sim_idx}.econparams.frequency, [burnin, sims{sim_idx}.firm_profits_diff_disc(sigma_idx,:)], patterns(sigma_idx), ...
                'DisplayName', strcat("$\sigma$ = ", num2str(sims{sim_idx}.sigmas(sigma_idx))))

        end

        hold off

        if sim_idx == 1
            legend show
            leg = legend;
            set(leg, 'AutoUpdate', 'off')
            legend boxoff
        end

        xlabel("Years, t+z")
        ylabel("$D_{\sigma,t,z}$", 'Interpreter', 'latex')
        xline_prop = xline(sims{sim_idx}.econparams.interval_lengths(1)*sims{sim_idx}.econparams.frequency, '--', 'color', [0.2, 0.2, 0.2]);
        xline_prop.LineWidth = 1;
        title(strcat("t = T - ", num2str(sims{sim_idx}.econparams.interval_lengths(1)*sims{sim_idx}.econparams.frequency - round((sims{sim_idx}.t - 1)*sims{sim_idx}.econparams.frequency,1))))
        grid off
        ylim([-0.015, 0.045])
        xlim([0,100])
        xticks([0,sims{sim_idx}.econparams.interval_lengths(1)*sims{sim_idx}.econparams.frequency,50,100])
        xticklabels([0,"T",50,100])

    end
    
    f.Units = 'centimeters';
    orient(f, 'landscape')
    f.PaperPosition = [-1, -1, 8, 4];
    saveas(f, "Output/Figures/Main/Figure_" + f.Name + ".png")
    

end



