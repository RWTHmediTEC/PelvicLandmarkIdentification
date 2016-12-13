clearvars
close all

load('data\pelvis.mat')

ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, 'vis', true);
