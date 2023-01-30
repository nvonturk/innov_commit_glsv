%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: update_econparams.m
% Author: Nicholas von Turkovich
% Date: 4/28/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [econparams] = update_econparams(econparams_orig, types, values)

    % If types and values are mismatched or empty, error out
    if isempty(types) || isempty(values) || length(types) ~= length(values)
        error("Issue with updating econparams")
    end
    
    econparams = econparams_orig;

    % For loop over alterations to the original econparams structure
    for i = 1:length(types)
        
        params_new = econparams.params;
        disp_flag_new = econparams.disp_flag;
        
        if types(i) == "heta"
            
            % Heta is in index 7 and is expecting an annual decimal
            params_new(7) = values(i);
            
        elseif types(i) == "rho"
            
            % Rho is in index 3 and is expecting an annual decimal
            params_new(3) = values(i);
            
        elseif types(i) == "frequency"
            
            % Frequency is in index 1 and is expecting 1/number of periods
            % in a year
            params_new(1) = values(i);
            
        elseif types(i) == "phi"
            
            % phi is in index 13 and is expecting a decimal share
            params_new(13) = values(i);   
            
        elseif types(i) == "gamma"
            
            % gamma is in index 5
            params_new(5) = values(i);
            
        elseif types(i) == "B"
            
            % B is in index 11 and is expecting an annual decimal
            params_new(11) = values(i);
            
        else
            error("Unknown parameter change")
        end
        
        econparams = generate_econparams(params_new, disp_flag_new);
        
    end
    
    if isfield(econparams_orig, "interval_lengths")
        econparams.interval_lengths = econparams_orig.interval_lengths;
    end
    
    if isfield(econparams_orig, "calibration_file")
        econparams.calibration_file = econparams_orig.calibration_file;
    end

end