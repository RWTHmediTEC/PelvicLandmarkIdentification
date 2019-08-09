clearvars; close all; opengl hardware

load('data/pelvis.mat')

[APP_TFM, Landmarks] = automaticPelvicCS(pelvis, 'visu', 1, 'debug', 0);

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';