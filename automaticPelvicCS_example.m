clearvars; close all; opengl hardware

load('data/pelvis.mat')

[fwTFMinput2APCS, CL_input] = automaticPelvicCS(pelvis, 'vis', true);

view(145,10)

