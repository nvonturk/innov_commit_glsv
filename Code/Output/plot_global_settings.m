%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: plot_global_settings.m
% Author: Nicholas von Turkovich
% Date: 11/17/2021
% Note(s): Sets a series of default plotting parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_global_settings()
    
    close all;
    
    set(0, 'DefaultLineLineWidth', 3);
    set(0, 'DefaultConstantLineLineWidth', 3);
    set(0, 'DefaultLineMarkerSize', 12)
    set(0, 'DefaultTextInterpreter', 'latex');
    set(0, 'DefaultLegendInterpreter', 'latex');
    set(0, 'DefaultAxesYGrid', 'on');
    set(0, 'DefaultTextFontSize', 12);
    set(0, 'DefaultAxesFontSize', 12);
    
end
