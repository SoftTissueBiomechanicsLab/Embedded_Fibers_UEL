function [BeamMesh, npmax]= CantileverBeamSegmentGen(W, nx, ny, nz, nh, SolidMesh, BeamMesh)

nB = BeamMesh.nB;
nS = SolidMesh.nS;
BCon = BeamMesh.Connectivity;
BNodes = BeamMesh.Nodes;
SCon = SolidMesh.Connectivity;  
SNodes = SolidMesh.Nodes;
SCentroids = SolidMesh.Centroids; % Brick Centroids
nBel = size(BCon,1);
nSel = size(SCon,1);
npmax = 0; % max # of bricks penetrated by a single beam element

% Arc Length Discretization
SolidElems = 0:W/nx:W;
BeamElems = 0:W/nh:W;

    
% Loop beam elements
for i = 1:nBel %Loop over Beam Elements

    SegmentMat = [];

    % 3D Position
    r0 = BNodes(BCon(i,1),:)';
    rf = BNodes(BCon(i,end),:)';

    t = rf - r0;
    le = norm(t); % Beam length
    t = t/norm(t); % Normalize tangent

    Bs0 = BeamElems(i);
    Bsf = BeamElems(i+1);

    SegmentInd = find((SolidElems<0.999*Bsf)&(SolidElems>1.0001*Bs0));

    S_points = [Bs0, SolidElems(SegmentInd), Bsf];

    for ii = 2:length(SegmentInd)+2

        rm = r0 + (S_points(ii-1) + (S_points(ii)-S_points(ii-1))/2 - Bs0)*t;
        %--------------------------
        % Find 6 closest centroids
        if nSel>6
            [IDX, ~] = knnsearch(SCentroids,rm','K',6);
        else
            IDX = 1:nSel;
        end

        ID0 = PointInBrick(rm, SCon(IDX,:), SNodes);
        % Transform from physical space to current element
        xi1 = (2*S_points(ii-1)-(Bs0+Bsf))/(Bsf-Bs0);
        xi2 = (2*S_points(ii)-(Bs0+Bsf))/(Bsf-Bs0);
        
        SegmentMat =[SegmentMat; xi1, xi2, IDX(ID0)];
    end

    BeamMesh.SegmentMat{i} = SegmentMat;
    
    if size(SegmentMat,1)>npmax
        npmax = size(SegmentMat,1);
    end

end

    