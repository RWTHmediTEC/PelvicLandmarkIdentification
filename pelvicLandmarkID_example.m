clearvars; close all; opengl hardware

% Clone example data
if ~exist('VSD', 'dir')
    try
    !git clone https://github.com/RWTHmediTEC/VSDFullBodyBoneModels VSD
    rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/RWTHmediTEC/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

addpath(genpath('src'))

% Select subjects of the VSD
Subjects = [1 9 13 19 23 24 27 35 36 42 46 49 50 55 56 57 61 62 64 66];
Subjects = arrayfun(@(x) ['z' num2str(x, '%03i')], Subjects', 'uni',0);

for s=1%:size(Subjects, 1)
    name = Subjects{s,1};
    
    load(['VSD\Bones\' name '.mat'],'B');
    % Construct the pelvic bone
    [pelvis.vertices, pelvis.faces] = concatenateMeshes(B(1:3).mesh);
    % Get landmarks and pelvic coordinate system
    [TFM2pelvicCS, Landmarks] = pelvicLandmarkID(pelvis, 'visu',1, 'CS','APP' ,'debug',0);
    set(gcf, 'Name',['Subject: ' name], 'NumberTitle', 'Off')
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';