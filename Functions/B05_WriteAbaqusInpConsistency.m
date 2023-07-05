function [] = B05_WriteAbaqusInpConsistency(WorkDir, JobNum, SolidMesh, BeamMesh, PenaltyConst, BeamRadius, param_beam, nSteps, Matl, param_solid, CoupleElemCon, nmax, nT, SolverOptions, FyTraction,nh1,nh2)

check_solid = false;

filename = [WorkDir, '/Job',num2str(JobNum,'%.4d'),'/Model',num2str(JobNum,'%.4d'),'.inp'];
fid = fopen(filename,'w');

%% GEOMETRY SECTION
BeamOrder = BeamMesh.Order;
nS = size(SolidMesh.Nodes,1);
nB = size(BeamMesh.Nodes,1);
nh = size(BeamMesh.Connectivity,1);

Nodes = [SolidMesh.Nodes; BeamMesh.Nodes];
ntotal = nS+nB;
GeometryData = [(1:ntotal)',Nodes];

% find free end bricks
[YCs,S1FaceBricks]=sort(SolidMesh.Centroids(:,3));
S1FaceBricks = S1FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S2FaceBricks]=sort(SolidMesh.Centroids(:,3),'descend');
S2FaceBricks = S2FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S6FaceBricks]=sort(SolidMesh.Centroids(:,1));
S6FaceBricks = S6FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S4FaceBricks]=sort(SolidMesh.Centroids(:,1),'descend');
S4FaceBricks = S4FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S3FaceBricks]=sort(SolidMesh.Centroids(:,2));
S3FaceBricks = S3FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S5FaceBricks]=sort(SolidMesh.Centroids(:,2),'descend');
S5FaceBricks = S5FaceBricks(1:length(YCs)/length(unique(YCs)));



W=max(Nodes(:,1))-min(Nodes(:,1));
D=max(Nodes(:,2))-min(Nodes(:,2));
H=max(Nodes(:,3))-min(Nodes(:,3));

fprintf(fid, '*NODE \n');
fprintf(fid, '\t %i, %2.14e, %2.14e, %2.14e \n', GeometryData');
%% Solid Node Sets
FixSolid = SolidMesh.negx;
n_in_set = length(FixSolid);
fprintf(fid, '*NSET, NSET=NegX \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, FixSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, FixSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, FixSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posx;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosX \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posy;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosY \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.negy;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=NegY \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posz;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosZ \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.negz;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=NegZ \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

% Beam Node Sets
fprintf(fid, '*NSET, NSET=FixBeamSet \n');
fprintf(fid, '%i \n', nS+1);
fprintf(fid, '*NSET, NSET=BeamEndPointSet \n');
fprintf(fid, '%i \n', ntotal);


%% Ref Node Sets
fprintf(fid, '*NSET, NSET=FixRefNode\n');
fprintf(fid, '%i \n', ntotal+1);
fprintf(fid, '*NSET, NSET=FreeRefNode\n');
fprintf(fid, '%i \n', ntotal+2);

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Xp_Face_S4\n');
fprintf(fid, '%i,\n', S4FaceBricks');
fprintf(fid, '*SURFACE, Name=Xp_Face\n');
fprintf(fid, 'Xp_Face_S4, S4\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Xn_Face_S6\n');
fprintf(fid, '%i,\n', S6FaceBricks');
fprintf(fid, '*SURFACE, Name=Xn_Face\n');
fprintf(fid, 'Xn_Face_S6, S6\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Yp_Face_S5\n');
fprintf(fid, '%i,\n', S5FaceBricks');
fprintf(fid, '*SURFACE, Name=Yp_Face\n');
fprintf(fid, 'Yp_Face_S5, S5\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Yn_Face_S3\n');
fprintf(fid, '%i,\n', S3FaceBricks');
fprintf(fid, '*SURFACE, Name=Yn_Face\n');
fprintf(fid, 'Yn_Face_S3, S3\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Zp_Face_S2\n');
fprintf(fid, '%i,\n', S2FaceBricks');
fprintf(fid, '*SURFACE, Name=Zp_Face\n');
fprintf(fid, 'Zp_Face_S2, S2\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Zn_Face_S1\n');
fprintf(fid, '%i,\n', S1FaceBricks');
fprintf(fid, '*SURFACE, Name=Zn_Face\n');
fprintf(fid, 'Zn_Face_S1, S1\n');
n_count = 0;
%% ElEMENT SECTION
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U001, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.A);
fprintf(fid, '1,2,3 \n'); 
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U002, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.B);
fprintf(fid, '1,2,3 \n');
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U003, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.C);
fprintf(fid, '1,2,3 \n');


% Bricks
fprintf(fid, '*ELEMENT, TYPE=C3D8RH, ELSET=BrickElems \n');
SCon = SolidMesh.Connectivity;
nSolEl = size(SCon,1); 
SolidConnectivity = [(1:nSolEl)',SCon];
fprintf(fid, ['%i', repmat(', %i',1,8), ' \n'], SolidConnectivity');

% Beams
BCon = BeamMesh.Connectivity + nS; % add nS for global connectivity
nBeamEl = size(BCon,1); 
nTotEl = nBeamEl + nSolEl;
BeamConnectivity = [(nSolEl+1:nTotEl)',BCon];
% fprintf(fid, '*ELEMENT, TYPE=B3%d , ELSET=BeamElems \n', BeamOrder);
% if BeamOrder == 1 || BeamOrder == 3
%     fprintf(fid, ['%i', repmat(', %i',1,2), ' \n'], BeamConnectivity');
% else
%     fprintf(fid, ['%i', repmat(', %i',1,3), ' \n'], BeamConnectivity');
% end
fprintf(fid, '*ELEMENT, TYPE=B3%dH, ELSET=BeamElem1\n', [BeamOrder]);
for i = 1:nh1
    if BeamOrder == 1 || BeamOrder == 3
        fprintf(fid, ['%i', repmat(', %i',1,2), ' \n'], BeamConnectivity(i,:)');
    else
        fprintf(fid, ['%i', repmat(', %i',1,3), ' \n'], BeamConnectivity(i,:)');
    end
end

fprintf(fid, '*ELEMENT, TYPE=B3%dH, ELSET=BeamElem2\n', [BeamOrder]);
for i = nh1+1:nh1+nh2
    if BeamOrder == 1 || BeamOrder == 3
        fprintf(fid, ['%i', repmat(', %i',1,2), ' \n'], BeamConnectivity(i,:)');
    else
        fprintf(fid, ['%i', repmat(', %i',1,3), ' \n'], BeamConnectivity(i,:)');
    end
end


                
    

% Define normals nq
fprintf(fid, '*NORMAL \n');
for i = 1:nh
    fprintf(fid, ['%d, %d, %.6f, %.6f, %.6f \n'], [BeamConnectivity(i,1),BeamConnectivity(i,2),BeamMesh.n2(i,:)]);
%     % Assign the same normal to the end point as the beginning of the next
%     % element
%     if i<nh
%         if BeamConnectivity(i,end)==BeamConnectivity(i+1,1)
%             fprintf(fid, ['%d, %d, %.6f, %.6f, %.6f \n'], [BeamConnectivity(i,1),BeamConnectivity(i,end),BeamMesh.n2(i+1,:)]);
%         end
%     end
            
end

% UELA
ncplA = nmax.A;
for i = 1:size(CoupleElemCon.A,1)
    fprintf(fid, '*ELEMENT, TYPE=U001 , ELSET=CplEl%d \n',i);    
    UEL_Con = [nTotEl+i, CoupleElemCon.A(i,:)];
    
    if (ncplA+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.A), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplA+1)/8);
        nrest = mod((ncplA+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% UELB
ncplB = nmax.B;
for i = 1:size(CoupleElemCon.B,1)
    fprintf(fid, '*ELEMENT, TYPE=U002 , ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1));    
    UEL_Con = [nTotEl+i+size(CoupleElemCon.A,1), CoupleElemCon.B(i,:)];
    
    if (ncplB+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.B), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplB+1)/8);
        nrest = mod((ncplB+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% UELC
ncplC = nmax.C;
for i = 1:size(CoupleElemCon.C,1)
    fprintf(fid, '*ELEMENT, TYPE=U003 , ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1));    
    UEL_Con = [nTotEl+i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1), CoupleElemCon.C(i,:)];
    
    if (ncplC+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.C), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplC+1)/8);
        nrest = mod((ncplC+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% if crosslinks.merge
%     
%     %fprintf(fid, '*ELEMENT, TYPE=CONN3D2, ELSET=Crosslinks \n');
%     
%     
%     % Loop through merged nodes groups
%     n_count = 0;
%     for i = 1:size(BeamMesh.Nodes2Merge,1)
%         
%         n_loc_merge = BeamMesh.Nodes2Merge{i, 2};
%         
%         for iii = 1:crosslinks.dof
%             fprintf(fid, '*ELEMENT, TYPE=SPRING2, ELSET=Crosslinks%d\n',iii);
%         
%             for ii = 1:n_loc_merge
% 
%                 n_count = n_count + 1;
%                 n_loc_nodes = BeamMesh.Nodes2Merge{i, 1} + nS;
% 
%                 if ii<n_loc_merge
% 
%                     fprintf(fid, '%d, %d, %d\n',[nTotEl+nh+n_count,n_loc_nodes(ii),n_loc_nodes(ii+1)]);
% 
%                 elseif ii>=n_loc_merge && n_loc_merge>2
% 
%                     fprintf(fid, '%d, %d, %d\n',[nTotEl+nh+n_count,n_loc_nodes(ii),n_loc_nodes(1)]);
% 
%                 end
%             end
%         end
%     end
% end

for i = 1:size(CoupleElemCon.A,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i);
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end

for i = 1:size(CoupleElemCon.B,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1));
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end

for i = 1:size(CoupleElemCon.C,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1));
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end


%% SECTIONS
fprintf(fid, '** \n');
fprintf(fid, '** SECTIONS \n');
fprintf(fid, '** \n');
if strcmp(Matl, 'NHK') || strcmp(Matl, 'OGD') || strcmp(Matl, 'MNRIV')
    fprintf(fid, '*SOLID SECTION, ELSET=BrickElems, material=SolidMatl \n');
    fprintf(fid, ',\n');
elseif strcmp(Matl, 'SVK')
    fprintf(fid, '*SOLID SECTION, ELSET=BrickElems, material=UANISO_SVK, orientation=DefaultOri \n');
    fprintf(fid, ',\n');
    fprintf(fid, '*Hourglass stiffness\n');
    fprintf(fid, '1.0,1.0,1.0\n');
else
    error('ERROR... MATERIAL NOT SUPPORTED (YET)')
end

fprintf(fid, '*BEAM SECTION, ELSET=BeamElem1, material=NeoHookeanBeam, section=CIRC \n');

  
    fprintf(fid, '%2.7e \n', BeamRadius);
    % fprintf(fid, '%2.14e, %2.14e, %2.14e \n', Na);
    fprintf(fid, '%.6f, %.6f, %.6f \n',BeamMesh.n1(1,:));


fprintf(fid, '*BEAM SECTION, ELSET=BeamElem2, material=NeoHookeanBeam, section=CIRC \n');

  
    fprintf(fid, '%2.7e \n', BeamRadius);
    % fprintf(fid, '%2.14e, %2.14e, %2.14e \n', Na);
    fprintf(fid, '%.6f, %.6f, %.6f \n',BeamMesh.n1(1,:));





%% Materials
fprintf(fid, '** \n');
fprintf(fid, '** Materials \n');
fprintf(fid, '** \n');

if strcmp(Matl, 'NHK')%Solid Neo Hookean
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, neo hooke \n');
    fprintf(fid, '%2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'OGD')% Ogden Incompressible
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, OGDEN, N=1 \n');
    fprintf(fid, '%2.3e, %2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'MNRIV')% Ogden Incompressible

    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, MOONEY-RIVLIN \n');
    fprintf(fid, '%2.3e, %2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'SVK') %Solid SVK
    
    fprintf(fid, '*orientation, name=DefaultOri, local directions=3 \n');
    fprintf(fid, '1.0,0.0,0.0,0.0,1.0,0.0 \n');
    fprintf(fid, '3,0.0 \n');
    fprintf(fid, '1,0,0 \n');
    fprintf(fid, '0,1,0 \n');
    fprintf(fid, '0,0,1 \n');
    
    fprintf(fid, '*Material, name=UANISO_SVK \n');
    fprintf(fid, '*density \n');
    fprintf(fid, '1e-12 \n');
    
    fprintf(fid, '*ANISOTROPIC HYPERELASTIC, USER, FORMULATION=STRAIN, TYPE=COMPRESSIBLE, PROPERTIES=9 \n');
%     fprintf(fid, '%2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e \n', param);
    fprintf(fid, '%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f \n', param);
    
end
    
% Beam Elastic Material
% fprintf(fid, '*Material, name=BeamMatl \n');
% fprintf(fid, '*Hyperelastic, neo hooke \n');
% fprintf(fid, '%2.3e, %2.3e \n', [2.5*mod_ratio, 0.6/mod_ratio]);
    
% Beam Elastic Material
% Beam Elastic Material
fprintf(fid, '*Material, name=NeoHookeanBeam \n');
fprintf(fid, '*Hyperelastic, neo hooke \n');
fprintf(fid, '%2.3e, %2.3e \n', param_beam);


% TIME POINTS
fprintf(fid,'*TIME POINTS, name=MustPoints, GENERATE\n');
fprintf(fid,'0., 1., %f\n',1/SolverOptions.n_incr);
%% STEP
InitStep = SolverOptions.InitStep;
minStep = SolverOptions.minStep;
maxStep = SolverOptions.maxStep;
StepTraction = linspace(0, FyTraction,nSteps);

for i = 1:nSteps
    fprintf(fid, '*Step, name=Step-%d, nlgeom=YES \n',i);
    if SolverOptions.STBL
        fprintf(fid, '*Static, STABILIZE=%.4e\n',SolverOptions.STBLfac);
    else
        fprintf(fid, '*Static \n');
    end
    
    fprintf(fid, '%.2e, 1., %.2e, %.2e \n', [InitStep, minStep, maxStep]);
    % Boundary Conditions
    fprintf(fid, '*Boundary \n');
    fprintf(fid, 'NegX, ENCASTRE \n');
    fprintf(fid, '*DLOAD \n');
    fprintf(fid, 'BeamElem1, PZ, %2.12e\n', StepTraction(i));
    fprintf(fid, '*DLOAD \n');
    fprintf(fid, 'BeamElem2, PZ, %2.12e\n', -StepTraction(i));
    
    % Output
    fprintf(fid, '*Restart, write, frequency=0 \n');
    fprintf(fid, '*Output, history, variable=PRESELECT, time points=MustPoints \n');
    fprintf(fid, '*Output, field, variable=PRESELECT, time points=MustPoints \n');
    fprintf(fid, '*ELEMENT Output\n');
    fprintf(fid, 'S,SE,SF,SM,SK \n');
    fprintf(fid, '*End Step \n');
    
end