function [D_g, k_g_inv]=AssembleDnK_v3(BeamMesh, LMOrder, nGP)
% D_lobal.........: M matrix per coupling element
% k_global_inv....: scaling matrix (inverse
% D_check.........: Non sparse format of global D 
% k_inv_check.....: Non sparse format of global k_global_inv
nh = size(BeamMesh.Connectivity,1);
BeamOrder = BeamMesh.Order;
BCon = BeamMesh.Connectivity;

D_g = sparse(3*size(BeamMesh.Nodes,1), 3*size(BeamMesh.Nodes,1));
k_g = sparse(3*size(BeamMesh.Nodes,1), 3*size(BeamMesh.Nodes,1));
k_g_inv = sparse(3*size(BeamMesh.Nodes,1), 3*size(BeamMesh.Nodes,1));


[xi_tilde, wGP] = lgwt(nGP,-1,1);

% Loop beam elements
for iel = 1:nh

    % Export local data for current beam element
    % NEED TO CHECK THIS OUT
    BeamNodeIDs = BCon(iel,:); % Because of the connectivity in  Couple element 
    LMNodeIDs = BeamNodeIDs; % To be changed in future
    
    r_elem = BeamMesh.Nodes(BeamNodeIDs,:)';

    d_count = 0;
    % Loop through local LM nodes
    for j = 1:length(LMNodeIDs)

        % Initialize K
        kjj = 0;

        %Loop through beam Nodes
        for l = 1:length(BeamNodeIDs)
            % Initialize
            d_count = d_count+1;
            Djlr = 0;
            Djlt = 0;

            %Gauss quadrature
            for igp = 1:nGP
                % Export data needed
                xi = xi_tilde(igp);
                           
                [Phi, dPhi] = ShapeFnc(xi, LMOrder);
                fd = (dPhi*r_elem');
                JA2 = (fd*fd')^(1/2);
                JB2 = 1;
                
%                 fprintf('detJ1 = %f',JA*JB)
%                 fprintf(' , detJ2 = %f \n', JA2*JB2)
                
                Hr = ShapeFnc(xi, BeamOrder);
                if BeamOrder~=3
                    Ht = zeros(BeamOrder+1,1); % For B31, B32 Elements
                    %Ht = Hr;
                elseif BeamOrder==3
                    temp = ShapeFnc(xi, BeamOrder);
                    Ht = [temp(3),temp(4)];
                end

                % Gauss integration
                Djlr = Djlr + JB2*JA2*wGP(igp)*Phi(j)*Hr(l);
                Djlt = Djlt + JB2*JA2*wGP(igp)*Phi(j)*Ht(l);

                if l==1 % kjj doesnot iterate through beam points
                    kjj = kjj + JB2*JA2*wGP(igp)*Phi(j);
                end

            end
             % Sanity Check
            ID_LM = LMNodeIDs(j);
            ID_B = BeamNodeIDs(l);
            row_id = (ID_LM-1)*3+1;
            col_id = (ID_B-1)*3+1;
            D_g(row_id:row_id+2 , col_id:col_id+2) = D_g(row_id:row_id+2 , col_id:col_id+2) + [Djlr*eye(3)];
   
            
        end
        % Check matrix (non sparse format)
        k_g(row_id:row_id+2,row_id:row_id+2) = k_g(row_id:row_id+2,row_id:row_id+2) + kjj*eye(3);
    end
end
k_g_inv = inv(k_g);
