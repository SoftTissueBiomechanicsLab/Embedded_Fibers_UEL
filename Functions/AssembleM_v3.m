function [Mg] = AssembleM_v3(BeamAxes, BeamMesh, SolidMesh, LMOrder, nGP)
% Mglobal.........: M matrix per coupling element
% CoupleElemCon...: Connectivity matrix of interface elements
% nsmax...........: Max. number of solid nodes associated with coupling
%                   element
%M_check..........: Non sparse format of global M

% Export necessary data
SolCon = SolidMesh.Connectivity;
SolNodes = SolidMesh.Nodes;

BCon = BeamMesh.Connectivity;
BeamOrder = BeamMesh.Order;

nh = size(BeamMesh.Connectivity,1);
nB = size(BeamMesh.Nodes,1);
nLM = nB;
nS = SolidMesh.nS;



[xi_tilde, wGP] = lgwt(nGP,-1,1);

% Initialize
Mg = sparse(3*nLM, 3*nS);

% Length Integral for test
ltest = 0;

% Loop beam elements
for iel = 1:nh

    % Export local data for current beam element
    SegMat = BeamMesh.SegmentMat{iel};
    nSegments = size(SegMat,1);
    
    BeamNodeIDs = BCon(iel,:);
    LMNodeIDs = BeamNodeIDs; % To be changed in future
    
    r_elem = BeamMesh.Nodes(BeamNodeIDs,:)';
   
    % Loop over segments
    for iseg = 1:nSegments

        %Export segment endpoints
        S_seg = SegMat(iseg,1:2);
        
        %InBrick = SegMat(iseg,3:end);
        InBrick = SegMat(iseg,3); % Keep just the first brick, in case beam is aligned with a brick edge

        BrickIDs = SolCon(InBrick,:);
        rk = SolNodes(BrickIDs,:); % Global coordinates of solid
       
        JB2 = (S_seg(2)-S_seg(1))/2;
        
        for igp = 1:nGP
          % Export data needed
          [xi, JB] = Segment2BeamXi_Straight(xi_tilde(igp), [-1,1], S_seg);
          
          %the actual beam domain
          [Phi, dPhi] = ShapeFnc(xi,  LMOrder);
          
          fd = (dPhi*r_elem');
          JA2 = (fd*fd')^(1/2);
                    
          ltest = ltest + JA2*JB2*wGP(igp);
          

          % Find global position r = SUM{ Hi(xi_beam)*r_i}
          if BeamOrder==3
              Hr = ShapeFnc(xi, BeamOrder);
              if BeamOrder~=3
                  Ht = zeros(BeamOrder+1,1); % For B31, B32 Elements
              elseif BeamOrder==3
                  temp = ShapeFnc(xi, BeamOrder);
                  Ht = [temp(3),temp(4)];
              end
              r = r_elem*Hr(1:2)';
              disp('VALID ONLY FOR INITIALLY STRAIGHT BEAMS')
          else
              r = r_elem*Phi';
          end


          % Inverse mapping from r -> (ksi, zeta, eta) in solid
          [xis,zetas,etas] = InverseMapping(r, rk, 1);
          

%           % Calculate Nk(ξ)
%           if xis<-1
%               xis = -1;
%           elseif zetas<-1
%               zetas = -1;
%           elseif etas<-1
%               etas = -1; 
%           elseif xis>1
%               xis = 1;
%           elseif zetas>1 
%               zetas = 1;
%           elseif etas>1
%               etas = 1;
%           end
          Nks = Hex8ShapeFnc(xis,zetas,etas);
         

          % Carry out the Guass integration and Assembly
          for j = 1:LMOrder+1
              row_ID = (BeamNodeIDs(j)-1)*3+1;
              for jj = 1:8
                  col_ID = (BrickIDs(jj)-1)*3+1;
                  Mg(row_ID:row_ID+2,col_ID:col_ID+2) = Mg(row_ID:row_ID+2,col_ID:col_ID+2) + JA2*JB*wGP(igp)*Phi(j)*Nks(jj)*eye(3); 
              end       
     
          end
        end

    end
end


fprintf('Total beam Length (approx. integr.)...%f \n',ltest)