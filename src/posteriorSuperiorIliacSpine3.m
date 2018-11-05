function PSIS = posteriorSuperiorIliacSpine3(pelvis, ASIS, varargin)
% Posterior superior iliac spine (PSIS)

% Parsing 
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
addParameter(p,'debugVisu', false, logParValidFunc);
parse(p,varargin{:});
visu = logical(p.Results.visualization);
debugVisu=logical(p.Results.debugVisu);

% Posterior superior iliac spine (PSIS)
pointProps.Linestyle = 'none';
pointProps.Marker = 'o';
pointProps.MarkerEdgeColor = 'r';
pointProps.MarkerFaceColor = 'r';

% Split pelvis in left and right hip bone
proxPelvis(1) = pelvis(1); % left hip bone 
proxPelvis(2) = pelvis(3); % right hip bone

% Height of the APP
APPheight=projPointOnLine3d([0 0 0],createLine3d(ASIS(1,:), ASIS(2,:)));

% Osteophytes or bridging in the PIIS region could be detected as PSIS
% point. For this reason, the size of the Z component of the PSIS vector is
% limited. This does not work if both PIIS regions are affected.

% The vector connecting the PSIS points
PSISvector=[nan nan Inf];
% Maximal absolute value of the z-component of the PSIS vector
MAX_Z_COMPONENT =0.11;
% Factor of the APP height to create the origin of the transverse cutting plane
CUT_FACTOR = 0.55;

% TODO: Actually, the while loop should start with a cut at the wrong PSIS
% point and be terminated before maximal height of the pelvis is reached.
while abs(PSISvector(3)) > MAX_Z_COMPONENT && CUT_FACTOR<=2
    CUT_FACTOR = CUT_FACTOR+0.05;
    % Transverse cutting plane
    transversePlane = [0 0 CUT_FACTOR*APPheight(3) 1 0 0 0 1 0];
    
    % Decide what side is going to be cutted
    if PSISvector(3) > MAX_Z_COMPONENT && PSISvector(3) ~=Inf
        Sides = 1;
    elseif PSISvector(3) < -MAX_Z_COMPONENT
        Sides = 2;
    elseif PSISvector(3) == Inf
        Sides =1:2;
    end
    % Cut the hip bones at the transverse plane
    for s=Sides
        proxPelvis(s) = cutMeshByPlane(proxPelvis(s), transversePlane,'part','above');
    end
    if debugVisu && visu
        patchProps.FaceColor='none';
        patchProps.EdgeColor='k';
        debugHandle = arrayfun(@(x) patch(x, patchProps), proxPelvis);
        delete(debugHandle)
    end
    
    % Get the most posterior point of the cutted hip bones
    yMinIdx=nan(1,2);
    PSIS=zeros(2,3);
    for s=1:2
        [~, yMinIdx(s)] = min(proxPelvis(s).vertices(:,2));
        PSIS(s,:) = proxPelvis(s).vertices(yMinIdx(s),:);
    end
    if debugVisu && visu
        debugHandle = drawPoint3d(PSIS, pointProps);
        delete(debugHandle)
    end
    % Calculate the direction of the vector connecting the PSIS points
    PSISvector=normalizeVector3d(PSIS(2,:)-PSIS(1,:));
end

%% Visualization
if visu
    drawPoint3d(PSIS, pointProps)
    text(PSIS(:,1), PSIS(:,2), PSIS(:,3), 'PSIS','FontWeight','bold',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

end