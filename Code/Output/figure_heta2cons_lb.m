%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_heta2cons_lb.m
% Author: Nicholas von Turkovich
% Date: 06/8/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_heta2cons_lb(filestart_data, econparams)
    
    f = figure("name", "heta2cons_lb");
    
    high_heta1_con = load(filestart_data + "_hub400_consistent.mat");    
    
    subplot(2,2,1)
    plot(high_heta1_con.heta_grid_interp*100/econparams.frequency, high_heta1_con.heta2_consistent_interp*100/econparams.frequency, 'k-')
    xlabel("$\eta_1$")
    ylabel("$\eta_2$")
    subplot(2,2,2)
    plot(high_heta1_con.heta_grid_interp(2:end)*100/econparams.frequency, diff(high_heta1_con.heta2_consistent_interp)*100/econparams.frequency, 'k-')
    xlabel("$\eta_1$")
    ylabel("$\Delta\eta_2$")
    subplot(2,2,3)
    plot(high_heta1_con.heta_grid_interp*100/econparams.frequency, high_heta1_con.heta2_consistent_interp*100/econparams.frequency, 'k-')
    xlabel("$\eta_1$")
    ylabel("$\eta_2$")
    xlim([high_heta1_con.heta_grid_interp(end/2), high_heta1_con.heta_grid_interp(end)]*100/econparams.frequency)
    subplot(2,2,4)
    plot(high_heta1_con.heta_grid_interp(2:end)*100/econparams.frequency, diff(high_heta1_con.heta2_consistent_interp)*100/econparams.frequency, 'k-')
    xlabel("$\eta_1$")
    ylabel("$\Delta\eta_2$")
    xlim([high_heta1_con.heta_grid_interp(end/2), high_heta1_con.heta_grid_interp(end)]*100/econparams.frequency)


end



