function [BeamMesh]= BeamMeshGenerator_v2(BeamAxes, BeamOrder)
% Generates Beam Mesh depending on anything dependable
n_beams = size(BeamAxes,1);

%Export Gauss points for length calculations (quadratic elems)
nGP = 12;
[xis, wGP] = lgwt(nGP,-1,1);
[Phi, dPhi] = ShapeFnc(xis,  BeamOrder);

% Loop over beams
for i = 1:n_beams
    
    BeamDataLoc = {BeamAxes{i,:}};
    hbeam = BeamDataLoc{9};
    lb = CurvePoint(BeamDataLoc); % Beam length
    nh = max([round(lb/hbeam),1]); % Beam elements
    le = lb/nh; % Length element
    r0 = BeamDataLoc{1}';
    rf = BeamDataLoc{2}';

    % Loop over elements
    SParam = zeros(nh,BeamOrder+1);
    ArcLength = 0;
    FE_ArcLength = 0;
    FESParam = zeros(nh,BeamOrder+1);
    [Nodes,~,~] = CurvePoint(BeamDataLoc, ArcLength);
    Nodes = Nodes';
    Connectivity = [];
    
    for ii = 1:nh
    
        if BeamOrder == 1 || BeamOrder==3%Linear Elements OR EulerBernoulli
                        
            SParam(ii,:) = [ArcLength,ArcLength + le];
            ArcLength = ArcLength + le;
            [rc,~,~] = CurvePoint(BeamDataLoc, ArcLength);
            
            % Update discretized arc length
            FEle = norm(Nodes(end,:)-rc');
            FESParam(ii,:) = [FE_ArcLength,FE_ArcLength + FEle];
            FE_ArcLength = FE_ArcLength + FEle;
            
            % Update nodes list and connectivity        
            Nodes = [Nodes; rc'];
            Connectivity = [Connectivity; ii, ii+1];
            
        elseif BeamOrder == 2 %Quadratic
            SParam(ii,:) = [ArcLength,ArcLength+le/2,ArcLength + le];
            MidArcLength =  ArcLength + le/2;
            ArcLength = ArcLength + le;
             
            [rm,~,~] = CurvePoint(BeamDataLoc, MidArcLength);
            [rc,~,~] = CurvePoint(BeamDataLoc, ArcLength);
            
            % Update discretized arc length
            locNodes = [Nodes(end,:); rm'; rc']';      
            fd = (dPhi*locNodes');
            JA2 = (sum(fd.^2,2)).^(1/2);
            FEle = sum(JA2.*wGP);
            
            FESParam(ii,:) = [FE_ArcLength,FE_ArcLength+FEle/2,FE_ArcLength + FEle];
            FE_ArcLength = FE_ArcLength + FEle;
            
            
            % Keep Nodes and Connectivity
            Nodes = [Nodes; rm'; rc'];
            istart = (ii-1)*2+1;
            Connectivity = [Connectivity; istart, istart+1, istart+2];
        end
        
    end
% HERE FOR BEAM NETWORKS. Will need to figure out sth for common nodes for
% the SParam vector, since if I merge the nodes then will also have to
% include indices as well there.
BeamMesh(i).Nodes = Nodes;
BeamMesh(i).SParam = FESParam;
BeamMesh(i).PhysicalSParam = SParam;
BeamMesh(i).Connectivity = Connectivity;
BeamMesh(i).nB = size(Nodes,1);
BeamMesh(i).Order = BeamOrder;

end





    
    
    
    
    
    
    
    
