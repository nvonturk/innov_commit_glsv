%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: custom_contour.m
% Author: Nicholas von Turkovich
% Date: 4/26/2022
% Note(s): 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [contour_array] = custom_contour(M, heta1_bounds, heta2_bounds)
    
    idx_start = 1;
    idx_end = -1;
    iter = 1;
    
    contour_array = {};
    
    while idx_end < size(M,2)
        
        info_col = M(:,idx_start);
        idx_end = idx_start + info_col(2);
        
        xcoords = M(1,idx_start+1:idx_end);
        ycoords = M(2,idx_start+1:idx_end);
        
        toremove = ((xcoords < heta1_bounds(1)) + (xcoords > heta1_bounds(2)) + ...
            (ycoords < heta2_bounds(1)) + (ycoords > heta2_bounds(2))) > 0;
        xcoords = xcoords(~toremove);
        ycoords = ycoords(~toremove);
        
        if length(xcoords) > 1
            contour_array{iter} = unique([xcoords',ycoords'], 'rows', 'stable');
            iter = iter + 1;
        end
        
        idx_start = idx_end + 1;
    end
    

end