function PSIS = posteriorSuperiorIliacSpine3(pelvis,varargin)
% Posterior superior iliac spine (PSIS) % iliac crest (IC) detection

parser = inputParser;
addOptional(parser,'visualization',true,@islogical);
parse(parser,varargin{:});
visu = parser.Results.visualization;

% Cut the pelvic bone along the sagittal plane
proxPelvis(1) = pelvis(1);
proxPelvis(2) = pelvis(3);

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