clearvars; close all; opengl hardware
% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

load('data\pelvis.mat')

LM = pelvicLandmarks(pelvis, ASIS, 'vis', true);