clearvars; close all; opengl hardware

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))

%% Data
DB_Path = '..\..\KlinikumDoDatabase\';
load([DB_Path 'SegmentationDatabase.mat'], 'DB')

% Remove some subjects
rmSubjects={'20500577'};
DB(ismember({DB.ID},rmSubjects))=[];

for p=9%1:length(DB)
    %% Load data
    load([DB_Path DB(p).ID '.mat']);
    
    display([num2str(p) ' Processing ID ' DB(p).ID])
    
    clearvars -except DB_Path DB Info Bone plotLandmarks
    
    % Construct the pelvic bone
    if isequal(Info.NoB, 1:5)
        pelvis = concatenateMeshes([Bone(3:5).MeshAPCS]);
    elseif isequal(Info.NoB, [1 2 6])
        pelvis = Bone(6).MeshAPCS;
    else
        error('No pelvic bone')
    end
    
    % Already detected landmarks
    ASIS = Info.ClinicalLandmarks.ASIS;
    ClinicalLandmarks = pelvicLandmarks(pelvis, ASIS, true);
end