%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: tables.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = tables(specification, econparams, optcons_duple_transition, optcons_duple_terminal, duple_commitment_transition, duple_commitment_terminal,...
    optcons_triple_transition, optcons_triple_terminal, triple_commitment_transition, triple_commitment_terminal, econparams_presclerosis, commitment_statedependent, consistent_statedependent)
    
    %% Fill in calibration table

    create_calibration_table("calibration", "calibration", econparams);
    create_calibration_table("calibration", "calibration_presclerosis", econparams_presclerosis);
    
    %% Fill in dynamics table 2
    dynamics_table_v2 = fileread("Output/Tracked/table_dynamics_template_v2.tex");
    T = econparams.interval_lengths(1);
    intervals = [T, 2*T, 3*T, -1];
    intervals_strings = ["T", "2T", "3T", "BGP"];

    for t = 1:length(intervals)

        if t < length(intervals)

            if t == 1
                idx_start = 1;
                idx_end = intervals(t);
            else
                idx_start = intervals(t-1)+1;
                idx_end = intervals(t);
            end

            % Extract the patent expiry rate for each transition type
            heta_val = round(100*mean(optcons_duple_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "cons", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*mean(duple_commitment_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "dup", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*mean(optcons_triple_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "tripcons", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*mean(triple_commitment_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "trip", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val = round(mean(100*(power(optcons_duple_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "cons", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(mean(100*(power(duple_commitment_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "dup", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(mean(100*(power(optcons_triple_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "tripcons", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(mean(100*(power(triple_commitment_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "trip", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            % Extract the mean markup for each transition type
            meanmarkup_val = round(mean(100*(optcons_duple_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "cons", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(mean(100*(duple_commitment_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "dup", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(mean(100*(optcons_triple_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "tripcons", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(mean(100*(triple_commitment_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "trip", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

        elseif t == length(intervals)

            % Extract the patent expiry rate for each transition type
            heta_val = round(100*optcons_duple_terminal.heta/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "cons", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*duple_commitment_terminal.heta/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "dup", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*optcons_triple_terminal.heta/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "tripcons", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            heta_val = round(100*triple_commitment_terminal.heta/econparams.frequency, 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("eta", "trip", intervals_strings(t)), num2str(heta_val, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val = round(100*(power(optcons_duple_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "cons", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(100*(power(duple_commitment_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "dup", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(100*(power(optcons_triple_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "tripcons", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            prodg_val = round(100*(power(triple_commitment_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("prodg", "trip", intervals_strings(t)), num2str(prodg_val, '%2.2f'));

            % Extract the mean markup for each transition type
            meanmarkup_val = round(100*(optcons_duple_terminal.mean_gross_markup-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "cons", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(100*(duple_commitment_terminal.mean_gross_markup-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "dup", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(100*(optcons_triple_terminal.mean_gross_markup-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "tripcons", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));

            meanmarkup_val = round(100*(triple_commitment_terminal.mean_gross_markup-1), 2);
            dynamics_table_v2 = strrep(dynamics_table_v2, strcat("meanmarkup", "trip", intervals_strings(t)), num2str(meanmarkup_val, '%2.2f'));
                
        else
            
            error("Column count off")
            
        end

    end

    fid = fopen("Output/Tables/table_" + specification + "_dynamics_out_v2.tex", 'w');
    fprintf(fid, '%s', dynamics_table_v2);
    fclose(fid);

    %% Fill in dynamics table 3 (single reopt and double reopt versions)
    dynamics_table_v3 = fileread("Output/Tracked/table_dynamics_template_v3.tex");
    T = econparams.interval_lengths(1);
    intervals = [T, 2*T, -1];
    intervals_strings = ["T", "2T", "BGP"];

    for t = 1:length(intervals)

        if t < length(intervals)

            if t == 1
                idx_start = 1;
                idx_end = intervals(t);
            else
                idx_start = intervals(t-1)+1;
                idx_end = intervals(t);
            end

            % Extract the patent expiry rate for each transition type
            heta_val_cons = round(100*mean(optcons_duple_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);        
            heta_val_commit = round(100*mean(duple_commitment_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("eta", "", intervals_strings(t)), num2str(heta_val_cons - heta_val_commit, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val_cons = round(mean(100*(power(optcons_duple_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            prodg_val_commit = round(mean(100*(power(duple_commitment_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("prodg", "", intervals_strings(t)), num2str(prodg_val_cons - prodg_val_commit, '%2.2f'));

            % Extract the mean markups for each transition type
            meanmarkup_val_cons = round(mean(100*(optcons_duple_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);        
            meanmarkup_val_commit = round(mean(100*(duple_commitment_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);        
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("meanmarkup", "", intervals_strings(t)), num2str(meanmarkup_val_cons - meanmarkup_val_commit, '%2.2f'));

        elseif t == length(intervals)

            % Extract the patent expiry rate for each transition type
            heta_val_cons = round(100*optcons_duple_terminal.heta/econparams.frequency, 2);        
            heta_val_commit = round(100*duple_commitment_terminal.heta/econparams.frequency, 2);
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("eta", "", intervals_strings(t)), num2str(heta_val_cons - heta_val_commit, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val_cons = round(100*(power(optcons_duple_terminal.g+1,1/econparams.frequency)-1), 2);        
            prodg_val_commit = round(100*(power(duple_commitment_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("prodg", "", intervals_strings(t)), num2str(prodg_val_cons - prodg_val_commit, '%2.2f'));

            % Extract the mean markups for each transition type
            meanmarkup_val_cons = round(100*(optcons_duple_terminal.mean_gross_markup-1), 2);
            meanmarkup_val_commit = round(100*(duple_commitment_terminal.mean_gross_markup-1), 2);
            dynamics_table_v3 = strrep(dynamics_table_v3, strcat("meanmarkup", "", intervals_strings(t)), num2str(meanmarkup_val_cons - meanmarkup_val_commit, '%2.2f'));

        else
            error("Column count off")
        end

    end

    fid = fopen("Output/Tables/table_" + specification + "_dynamics_out_v3.tex", 'w');
    fprintf(fid, '%s', dynamics_table_v3);
    fclose(fid);

    dynamics_table_v4 = fileread("Output/Tracked/table_dynamics_template_v3.tex");
    T = econparams.interval_lengths(1);
    intervals = [T, 2*T, -1];
    intervals_strings = ["T", "2T", "BGP"];

    for t = 1:length(intervals)

        if t < length(intervals)

            if t == 1
                idx_start = 1;
                idx_end = intervals(t);
            else
                idx_start = intervals(t-1)+1;
                idx_end = intervals(t);
            end

            % Extract the patent expiry rate for each transition type
            heta_val_cons = round(100*mean(optcons_triple_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);        
            heta_val_commit = round(100*mean(triple_commitment_transition.heta_path(idx_start:idx_end))/econparams.frequency, 2);
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("eta", "", intervals_strings(t)), num2str(heta_val_cons - heta_val_commit, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val_cons = round(mean(100*(power(optcons_triple_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            prodg_val_commit = round(mean(100*(power(triple_commitment_transition.prodg(idx_start:idx_end)+1,1/econparams.frequency)-1)), 2);
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("prodg", "", intervals_strings(t)), num2str(prodg_val_cons - prodg_val_commit, '%2.2f'));

            % Extract the mean markups for each transition type
            meanmarkup_val_cons = round(mean(100*(optcons_triple_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);        
            meanmarkup_val_commit = round(mean(100*(triple_commitment_transition.mean_gross_markup(idx_start:idx_end)-1)), 2);        
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("meanmarkup", "", intervals_strings(t)), num2str(meanmarkup_val_cons - meanmarkup_val_commit, '%2.2f'));

        elseif t == length(intervals)

            % Extract the patent expiry rate for each transition type
            heta_val_cons = round(100*optcons_triple_terminal.heta/econparams.frequency, 2);        
            heta_val_commit = round(100*triple_commitment_terminal.heta/econparams.frequency, 2);
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("eta", "", intervals_strings(t)), num2str(heta_val_cons - heta_val_commit, '%2.2f'));

            % Extract the growth rate for each transition type
            prodg_val_cons = round(100*(power(optcons_triple_terminal.g+1,1/econparams.frequency)-1), 2);        
            prodg_val_commit = round(100*(power(triple_commitment_terminal.g+1,1/econparams.frequency)-1), 2);
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("prodg", "", intervals_strings(t)), num2str(prodg_val_cons - prodg_val_commit, '%2.2f'));

            % Extract the mean markups for each transition type
            meanmarkup_val_cons = round(100*(optcons_triple_terminal.mean_gross_markup-1), 2);
            meanmarkup_val_commit = round(100*(triple_commitment_terminal.mean_gross_markup-1), 2);
            dynamics_table_v4 = strrep(dynamics_table_v4, strcat("meanmarkup", "", intervals_strings(t)), num2str(meanmarkup_val_cons - meanmarkup_val_commit, '%2.2f'));

        else
            error("Column count off")
        end

    end

    fid = fopen("Output/Tables/table_" + specification + "_dynamics_out_v4.tex", 'w');
    fprintf(fid, '%s', dynamics_table_v4);
    fclose(fid);
    
    %% Fill in state dependency table
    statedep_table = fileread("Output/Tracked/table_statedep_template.tex");
    parameter_names = ["eta1comp", "eta1unc", "eta2"];
    
    for i = 1:length(parameter_names)
        
        statedep_table = strrep(statedep_table, strcat(parameter_names(i), "commit"), num2str(commitment_statedependent.heta_opt_duple_statedep(i)*100/econparams.frequency,"%2.2f"));
        statedep_table = strrep(statedep_table, strcat(parameter_names(i), "cons"), num2str(consistent_statedependent.heta_optcons_duple_statedep(i)*100/econparams.frequency,"%2.2f"));

    end
    
    fid = fopen("Output/Tables/table_" + specification + "_statedep.tex", 'w');
    fprintf(fid, '%s', statedep_table);
    fclose(fid);

end



