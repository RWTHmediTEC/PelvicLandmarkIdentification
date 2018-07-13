function PSIS = posteriorSuperiorIliacSpine3(pelvis, ASIS, varargin)
% Posterior superior iliac spine (PSIS) % iliac crest (IC) detection

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});
visu = parser.Results.visualization;

proxPelvis(1) = pelvis(1); % left hip
proxPelvis(2) = pelvis(3); % right hip

% Sagittal plane
sagittalPlane = [0 0 0 0 1 0 0 0 1];
% Height of the APP
APPheight = intersectLinePlane(createLine3d(ASIS(1,:), ASIS(2,:)), sagittalPlane);
% Transverse cutting plane
CUT_FACTOR = 0.6;
transversePlane = [0 0 CUT_FACTOR*APPheight(3) 1 0 0 0 1 0];
for s=1:2
    proxPelvis(s) = cutMeshByPlane(proxPelvis(s), transversePlane,'part','above');
end

yMinIdx=zeros(1,2);
PSIS=zeros(2,3);
for s=1:2
    [~, yMinIdx(s)] = min(proxPelvis(s).vertices(:,2));
    PSIS(s,:) = proxPelvis(s).vertices(yMinIdx(s),:);
end


%% Visualization
if visu == true
    % Posterior superior iliac spine (PSIS)
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'r';
    pointProps.MarkerFaceColor = 'r';
    drawPoint3d(PSIS, pointProps)
    text(PSIS(:,1), PSIS(:,2), PSIS(:,3), 'PSIS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

end