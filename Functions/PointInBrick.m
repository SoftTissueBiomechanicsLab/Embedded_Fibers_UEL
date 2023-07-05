function ID = PointInBrick(r, SCon, SNodes)

ID = [];
nSel_loc = size(SCon,1); %number of elements to search 

%Scan all bricks
fv.faces = [1,4,2;4,3,2;3,7,2;7,6,2;6,8,5;7,8,6;5,8,1;8,4,1;3,4,7;4,8,7;6,5,1;1,2,6];

for i = 1:nSel_loc
    
    fv.vertices = SNodes(SCon(i,:),:);
%     trimesh(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),'facecolor', 'none', 'edgecolor', 'k')  
    
    if inpolyhedron(fv,r')
        
        ID = [ID, i];
        
    end    
end
