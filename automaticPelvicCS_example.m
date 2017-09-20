clearvars; close all; opengl hardware

load('data/pelvis.mat')

[fwTFMinput2APCS, CL_input] = automaticPelvicCS(pelvis, 'vis', true);

view(145,10)

[List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
List.f = List.f'; List.p = List.p';