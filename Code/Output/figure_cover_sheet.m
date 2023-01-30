%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: figure_cover_sheet.m
% Author: Nicholas von Turkovich
% Date: 05/31/2021
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f] = figure_cover_sheet(specification, econparams)
    
    Parameter = {'n'; 'rho'; 'lambda'; 'phi'; 'eta'; 'zeta'; 'tau'; 'tau_rd'; 'gamma'; 'B'; 'interval_lengths'};
    Unit = {'None'; 'Annual %'; 'None'; '% share'; 'Annual %'; '% share'; '%'; '%'; 'None'; 'Annual'; 'Years'};
    Value = num2str([econparams.n; (power(econparams.rho + 1, 1/econparams.frequency)-1)*100;...
        econparams.lambda; econparams.phi*100; econparams.heta*100/econparams.frequency;...
        econparams.zeta*100; econparams.tau*100; econparams.tau_rd*100; econparams.gamma; econparams.B/econparams.frequency], '%2.4f');
    Value = [Value; compose("%3.2f/%3.2f/%3.2f", round(econparams.interval_lengths*econparams.frequency, 1))];

    if specification == "calibration_entry"
        Parameter = [Parameter; 'phi_e'; 'B_e'];
        Unit = [Unit; '% share'; 'Annual'];
        Value = [Value; num2str([econparams.phi_e*100; econparams.B_e/econparams.frequency], '%2.4f')];
    end

    Value = cellstr(Value);

    T = table(Parameter, Value, Unit);
    f = figure("name", "cover sheet");
    uit = uitable(f, 'Data', T{:,:}, 'Units', 'Normalized', 'Position',[0, 0, 1, 1],...
        'ColumnName', T.Properties.VariableNames);
    uit.ColumnWidth = {150, 200, 100};

end



