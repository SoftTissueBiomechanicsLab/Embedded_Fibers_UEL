function [GlobalBeamMesh] = GlobalBeamMeshGen_v2(BeamMesh, dist_tol)
% Receives the meshed beam domains (branches) and generates a global mesh.
% Merges common nodes (aka crosslinks) (JUST the nodes, not the elements).
% For constant stress transfer test (consistency test) merging overlapping
% nodes is not applied.

n_beams = length(BeamMesh);

% Loop beam branches
elem_count = 0;
node_count = 0;
GlobalNodes = [];
GlobalCon = [];
GlobalSParam = [];
GlobalTParam = [];
GlobalSegmentMat = {};
IEN = []; %first column: beam branch ID per element
Global_n1 = [];
Global_n2 = [];

% Different n1, n2 for quadratic beams
if BeamMesh(1).Order==2
    [~,dN,dN2] = ShapeFnc(-1, 2); % Export only for the first node
end

for i = 1:n_beams
    
    % Export Beam branch data
    nB = BeamMesh(i).nB;
    Nodes = BeamMesh(i).Nodes  ;
    Connectivity = BeamMesh(i).Connectivity;  
    
    % Merge matrices
    GlobalNodes = [GlobalNodes;Nodes];
    GlobalCon = [GlobalCon; Connectivity+node_count];
    
    % Update nodecounter
    node_count = node_count+nB;

    % Global Element quantities (are NOT merged, no matter what)
    nel = size(Connectivity,1);
    GlobalSParam = [GlobalSParam;BeamMesh(i).SParam];
    IEN = [IEN;i*ones(nel,1)]; % Index of original beam branch
    
    for j = 1:nel
        elem_count = elem_count + 1;

%         % Export t-paramter at element nodes
%         t_loc = [];
%         for jj = 1:size(BeamMesh(i).SParam,2)
%             s_current = BeamMesh(i).SParam(j,jj);
%             [~,t,~] = CurvePoint(BeamDataLoc, s_current);
%             t_loc =[t_loc,t];
%         end
%         GlobalTParam = [GlobalTParam;t_loc];

        % Build Global Segment Matrix
        GlobalSegmentMat{elem_count} = BeamMesh(i).SegmentMat{j};

        % Compute Cross section Orientation
        s1 = BeamMesh(i).SParam(j,1);
        Node1 = Nodes(Connectivity(j,1),:)';
        Node2 = Nodes(Connectivity(j,2),:)';
        Node3 = Nodes(Connectivity(j,end),:)'; % for linear is Node2

        % Abaqus tangent vector
        t = Node2-Node1;
        t = t/norm(t);

        % Cross Section Orientation with no curvature
        n1_approx = [0,0,-1]';
%         
%         % Special cases
%         if strcmp(BeamDataLoc{8}, 'helix')
%             [Node1Test,tparam1,~] = CurvePoint(BeamDataLoc,s1);
%             r0 = BeamDataLoc{1}';
%             rf = BeamDataLoc{2}';
%             axis_length = norm(rf-r0);
%             ezloc = (rf-r0)/axis_length;
%             nloops = BeamDataLoc{7};
%             cparam = axis_length/(2*nloops*pi);
%             zparam = cparam*tparam1;
%             rcenter = r0+zparam*ezloc;
%             n1_approx = rcenter-Node1;
%             n1_approx = n1_approx/norm(n1_approx);
%         end
        n2 = cross(t,n1_approx);
        
        if norm(n2) > 1e-12
            n2 = n2/norm(n2);
            n1 = cross(n2,t);
            n1 = n1/norm(n1);
        else

            n1 = [-1,0,0]';
            n2 = cross(t,n1);
            n2 = n2/norm(n2);
        end
        
        % For curved beams of second degree, define normals as follows
        if BeamMesh(i).Order==2 
            Xs = [Node1,Node2,Node3]';
            fd = (dN*Xs);
            fdd = (dN2*Xs);
             if sum(fdd)~=0 && norm(fdd)>1e-12 % There's a chance the segment is straight
                 n1 = fdd/(fd*fd')^(1/2)-fd*(fd*fdd')/(fd*fd')^(3/2);
                 n1 = n1'/norm(n1);
                 n2 = cross(t,n1);
             end
        end
        if isnan(n1(1))
            error('CHECK GEOMETRY')
        end
        % Fill in the vectors
        Global_n1 = [Global_n1;n1'];
        Global_n2 = [Global_n2;n2'];
%         figure(5)
%         hold on
%         Xs = [Node1';Node2'; Node3';];
%         plot3(Xs(:,1),Xs(:,2),Xs(:,3),'ko-','linewidth',1.2);
%         quiver3(Node1(1),Node1(2),Node1(3),n1(1),n1(2),n1(3),0.1,'r', 'linewidth',1.2)
%         quiver3(Node1(1),Node1(2),Node1(3),n2(1),n2(2),n2(3),0.1,'g', 'linewidth',1.2)
%         quiver3(Node1(1),Node1(2),Node1(3),t(1),t(2),t(3),0.1,'b', 'linewidth',1.2)
%         axis equal
%         hold off
      
%         
    end      

end

% Export Data
GlobalBeamMesh.Nodes = GlobalNodes;
GlobalBeamMesh.SParam = GlobalSParam;
GlobalBeamMesh.TParam = GlobalTParam;
GlobalBeamMesh.Connectivity = GlobalCon;
GlobalBeamMesh.SegmentMat=GlobalSegmentMat;
GlobalBeamMesh.n1 = Global_n1;
GlobalBeamMesh.n2 = Global_n2;
GlobalBeamMesh.IEN = IEN;


%% Merge overlapping nodes if requested
Group2Merge = {};
counter = 0;
nBg = size(GlobalNodes,1);
for i = 1:nBg
    [~, D] = knnsearch(GlobalNodes(i,:), GlobalNodes);
    IDX = find(D<dist_tol); %IDX always includes [i]-node

    if length(IDX)>1 % two or more nodes to be merged
        counter = counter+1;
        Group2Merge{counter} = sort(IDX)';
    end
end

% Keep only unique Groups
a = Group2Merge(:);
[a1,b,c] = unique(cellfun(@char,a,'un',0));
lo = histc(c,1:max(c));
loo = lo(:) > 1;
UniqueGroups = [a(b(loo)), num2cell(lo(loo))];
% 
    % Replace that in connectivity matrix and Node matrix
    IDs = (1:nBg)';
    rows2erase = [];
    for i = 1:size(UniqueGroups,1)
        mergedIDs = UniqueGroups{i,1};
        n_merged = length(mergedIDs);
        fprintf('Merged %d beam nodes \n',n_merged)
        for ii = 2:n_merged
            IDs(mergedIDs(ii)) = mergedIDs(1);
            rows2erase = [rows2erase,mergedIDs(ii)];
        end
    end
    % Find Deleted nodes
    IDXs = find(IDs~=(1:nBg)');
    newIDs = IDs;
    newIDs(IDXs) = [];
    
    for i = 1:length(newIDs)
        ID2keep = find(newIDs(i)==IDs);
        IDs(ID2keep)=i;
    end
    
    %Finally, define new connectivity
    NewGlobalCon = IDs(GlobalCon);

    % Update Global nodes
    NewGlobalNodes = GlobalNodes;
    NewGlobalNodes(IDXs,:) = [];

    % Update Output structure
    GlobalBeamMesh.Nodes = NewGlobalNodes;
    if n_beams==1 && nel==1
        GlobalBeamMesh.Connectivity = NewGlobalCon';
    else
        GlobalBeamMesh.Connectivity = NewGlobalCon;
    end
GlobalBeamMesh.Nodes2Merge = UniqueGroups;





