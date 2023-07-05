function [CubeMesh]= CubeMeshGenerator(nx,ny,nz,W,D,H)
%Returns a structure that contains data for HEX mesh of the unit cube 
% n_elems^3 elements

x = 0:1/nx:1;
y = 0:1/ny:1;
z = 0:1/nz:1;
dx = 1/nx;
dy = 1/ny;
dz = 1/nz;
%Counters
c      = 0; %Nodes
epi_c  = 0; %Epi Nodes
endo_c = 0; %Endo Surface Nodes
base_c = 0; %Base surface Nodes
apex_c = 0; %Apex surface Nodes
posx_c = 0; %Positive x Surface Nodes
negx_c = 0; %Negative x Surface Nodes

for i = 1:nx+1
    for j = 1:ny+1
        for k = 1:nz+1
            c = c+1;
            Nodes(c,:) = [x(i) y(j) z(k)];
            
            %Define Nodesets for Surfaces
            if x(i) == 1
                posx_c = posx_c+1;
                posx(posx_c) = c;
            end
            if x(i) == 0
                negx_c = negx_c+1;
                negx(negx_c)= c;
            end
            if y(j) == 1
                apex_c = apex_c+1;
                posy(apex_c) = c;
            end
            if y(j) == 0
                base_c = base_c+1;
                negy(base_c) = c;
            end    
            if z(k) == 1
                epi_c = epi_c+1;
                posz(epi_c) = c;
            end    
            if z(k) == 0
                endo_c = endo_c+1;
                negz(endo_c) = c;
            end
            
        end
    end
end

c = 0;
Connectivity=[];
for i = 1:nx
    for j = 1:ny
        for k = 1:nz

            c = c+1;
            node1 = [x(i) y(j) z(k)];
            node5 = [x(i) y(j) z(k+1)];
            node4 = [x(i) y(j+1) z(k)];
            node8 = [x(i) y(j+1) z(k+1)];
            node2 = [x(i+1) y(j) z(k)];
            node6 = [x(i+1) y(j) z(k+1)];
            node3 = [x(i+1) y(j+1) z(k)];
            node7 = [x(i+1) y(j+1) z(k+1)];
            

            for l=1:8
               
                eval(['node =' 'node' num2str(l) ';']);
                [ID,~] = find(node(1)==Nodes(:,1) & node(2)==Nodes(:,2) & node(3)==Nodes(:,3));
                Connectivity(c,l)=ID;
               
            end                              
        end
    end
end

Nodes(:,1) = W*Nodes(:,1);
Nodes(:,2) = D*Nodes(:,2);
Nodes(:,3) = H*Nodes(:,3);

% Calculate centroids of bricks
for i = 1:size(Connectivity,1)
    temp = Nodes(Connectivity(i,:),:);
    SCentroids(i,:) = mean(temp);
end


CubeMesh.Nodes = Nodes;
CubeMesh.posx  = posx;
CubeMesh.negx  = negx;
CubeMesh.posy  = posy;
CubeMesh.negy  = negy;
CubeMesh.posz   = posz;
CubeMesh.negz  = negz;
CubeMesh.Connectivity = Connectivity;
CubeMesh.Centroids = SCentroids;
CubeMesh.nS = size(Nodes,1);



