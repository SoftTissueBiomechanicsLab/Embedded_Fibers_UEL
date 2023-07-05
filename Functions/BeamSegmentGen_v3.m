function [BeamMesh, npmax]= BeamSegmentGen_v3(W, D, H, nx, ny, nz, SolidMesh, BeamMesh, s_resol, plot_bool)

n_beams = length(BeamMesh);

% Export Solid Mesh Data
nS = SolidMesh.nS;
SCon = SolidMesh.Connectivity;  
SNodes = SolidMesh.Nodes;
SCentroids = SolidMesh.Centroids; % Brick Centroids
nSel = size(SCon,1);
npmax = 0; % max # of bricks penetrated by a single beam element
r_out  = []; % beam points out of solid domain

% Global Planes
nXYZ_Global = [nx,ny,nz];
WDH_Global = [W,D,H];

% Function handle to retrieve only IDX-component of the 3D curve
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-10);
function y = CurveComp(xis, Nodes , IDX)
    n = length(xis);
    y = zeros(n,1);
    nOrder = size(Nodes,1)-1;
    Phis = ShapeFnc(xis,nOrder);
    rtemp = Phis*Nodes;
    y(:,1) = rtemp(:,IDX);
end

for i = 1:n_beams
    elem_cnt(i) = size(BeamMesh(i).Connectivity,1);
end
tot_elem_cnt = sum(elem_cnt);

% Initi progress bar
bar_count = 0;
c1 = 0;
fprogressbar = waitbar(0,'Beam Segmentation...');
for i = 1:n_beams
    
            
    % Export Current Beam Domain Data
    nB = BeamMesh(i).nB;
    BCon = BeamMesh(i).Connectivity;
    nBOrder = size(BCon,2)-1;
    BNodes = BeamMesh(i).Nodes;
    BSParam = BeamMesh(i).SParam;
    nBel = size(BCon,1);
    %BeamDataLoc = {BeamAxes{i,:}};
    %hbeam = BeamDataLoc{9};
    nelems = size(BCon,1);
    out_bool = false;
    
    for jel = 1:nelems
        
        % Progress bar
        bar_count = bar_count +1;
        waitbar(bar_count/tot_elem_cnt,fprogressbar)
        
        % Export local element geometry
        locNodes = BNodes(BCon(jel,:),:);
        
        % Define scan interval
        Xi_Star_Points = [];
        xitot = 2; % searching on the ksi domain
        xi_intrv = xitot/s_resol;
        xi_scan = [-1:xi_intrv:1];

        % Loop through Global planes
        for ii = 1:3
            n_planes = nXYZ_Global(ii);
            x = (0:1/n_planes:1)*WDH_Global(ii); % sections
            for j = 1:n_planes+1
                x0 = x(j);
                Intersect_eqn = @(xi) CurveComp(xi,locNodes, ii)- x0;

                % Step 1: Look for sign changes in the curve-plane intersection

                f_scan = Intersect_eqn(xi_scan);

                % Step 2: If sign changes, solve and keep the solution
                for jj = 1:length(f_scan)-1
                    if f_scan(jj)*f_scan(jj+1)<=0

    %                     [s_star,~,exitflag] = fsolve(@(s)Intersect_eqn , 5,options);
                        [xi_star,~,exitflag] = fzero(@(x)CurveComp(x,locNodes, ii)- x0, [xi_scan(jj), xi_scan(jj+1)]);

                        if xi_star<=-1 || xi_star>=1
                            exitflag=0;
                        end

                        % If solution is okay, keep it
                        if exitflag==1
                            Xi_Star_Points = [Xi_Star_Points;xi_star];
                        end
                    end
                end
            end
        end
    
%     if isempty(S_Star_Points)
%         % The entire beam branch in a single brick
%         S_Star_Points = [0, stot];
%     end
        % Sort Segments and keep unique
        Xi_points = uniquetol(Xi_Star_Points,xi_intrv*1e-6);

        SegmentMat = [];
        SegmentInd = find((Xi_points<0.999)&(Xi_points>-1.0001));
       
        Xi_points_loc = [-1, Xi_points(SegmentInd)', 1];
    
        for ii = 2:length(SegmentInd)+2
            Xmid = Xi_points_loc(ii-1) + (Xi_points_loc(ii)-Xi_points_loc(ii-1))/2;
            rm = (ShapeFnc(Xmid,nBOrder)*locNodes)';

            % Find 6 closest centroids
            if nSel>6
                [IDX, ~] = knnsearch(SCentroids,rm','K',6);
            else
                IDX = 1:nSel;
            end

            ID0 = PointInBrick(rm, SCon(IDX,:), SNodes);

            if ~isempty(ID0)
                SegmentMat =[SegmentMat; Xi_points_loc(ii-1), Xi_points_loc(ii), IDX(ID0(1))];
            else
                fprintf('\nWARNING: Beam outside solid domain, Branch %d\n',i)
                r_out = [r_out; rm'];
                out_bool = true;
            end
        end
        
        if out_bool
            c1 = c1 + 1;
        end

        
        % Store Segment Matrix
        BeamMesh(i).SegmentMat{jel} = SegmentMat;
        
        % Also, check max number of bricks penetrated by a single element
        if size(SegmentMat,1)>npmax
            npmax = size(SegmentMat,1);
        end
    
    


    %% PLOT
        if plot_bool
            figure(1)
            hold on
            if i==1
                
                fv.faces = [1,4,2;4,3,2;3,7,2;7,6,2;6,8,5;7,8,6;5,8,1;8,4,1;3,4,7;4,8,7;6,5,1;1,2,6];
                for j = 1:nSel
                    fv.vertices = SNodes(SCon(j,:),:);
                %     trimesh(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),'facecolor', 'none', 'edgecolor', 'k')
                    aa = [1,5,7,3,2,6,8,4];
                    fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
                    patch('Vertices',fv.vertices,'Faces',fac,'FaceColor','#ffa527','Facealpha',0.1)
                    %plot3(SCentroids(j,1),SCentroids(j,2),SCentroids(j,3),'ro','markersize',10);
                end
            end

            xiv = -1:xi_intrv/10:1;
            rcv = ShapeFnc(xiv,nBOrder)*locNodes;
            plot3(rcv(:,1),rcv(:,2),rcv(:,3),'.','color','#B22727')


            rcv = ShapeFnc(Xi_points,nBOrder)*locNodes;
            plot3(rcv(:,1),rcv(:,2),rcv(:,3),'x','MarkerSize',20,'LineWidth',2,'color','#006E7F')


            rcv = ShapeFnc([-1,0,1],nBOrder)*locNodes;
            plot3(rcv(:,1),rcv(:,2),rcv(:,3),'s','MarkerSize',20,'LineWidth',2,'color','#55aaff')

            axis equal
            view(17,35)
            axis off
            box off
        end
        
    end

end

if ~isempty(r_out)
    fprintf('Beam Nodes Outside of solid domain')
    % plot points out of solid
    %plot3(r_out(:,1),r_out(:,2),r_out(:,3),'cx','markersize',10)
end
hold off

close(fprogressbar)
close all
end