function Energy = create_vtkRPL(cur_path, JobNum, AbaqusData, SolidMesh, GlobalBeamMesh, nplotbeams, e_pen, M_ns, D_ns, kinv_ns, LM_bool, BCS,RPL)

% Create paths needed
% Define paths
IDStr = num2str(JobNum,'%.4d');
full_workdir = [cur_path,'/AbaqusWorkDir/Job',IDStr,'/'];
inp_dir = [full_workdir,'Model',IDStr,RPL,'.inp'];
vtkDir = [full_workdir,'VTKs',RPL,'/'];
if ~exist(vtkDir, 'dir')
    mkdir(vtkDir)
else
    % Remove old results
    eval(['rmdir ' vtkDir ' s']);
    % Re-make directory
    mkdir(vtkDir)
end

solid_file = [vtkDir,'S_',IDStr];
beam_file = [vtkDir,'B_',IDStr];
seg_file = [vtkDir,'Seg_',IDStr];
elem_file = [vtkDir,'Elem_',IDStr];


%% Read Abaqus input file
fid2 = fopen(inp_dir,'r');
line2 = fgetl(fid2);

flagnode=0;
flagbrick=0;
flagtruss=0;
C3D8R = [];
C3D8I = [];
P = [];
B32 = [];
B31 = [];

while line2>=0

    if contains(line2,'*NODE')
        temp = fscanf(fid2,'%d,%f,%f,%f\n'); 
        temp = reshape(temp,[4 length(temp)/4 ]);
        temp = temp';
        P = [P;temp];
        temp = []; 
    end     
    
    if contains(line2,'C3D8R')
        temp = fscanf(fid2,'%d,%f,%f,%f,%f,%f,%f,%f,%f\n');
        temp = reshape(temp,[9 length(temp)/9]);
        temp = temp';
        C3D8R = [C3D8R; temp];
        temp = [];         
    end 
    
    if contains(line2,'C3D8I')
        temp = fscanf(fid2,'%d,%f,%f,%f,%f,%f,%f,%f,%f\n');
        temp = reshape(temp,[9 length(temp)/9]);
        temp = temp';
        C3D8I = [C3D8I; temp];
        temp = [];         
    end 
    if contains(line2,'B31')
        temp = fscanf(fid2,'%d,%f,%f\n');
        temp = reshape(temp,[3 length(temp)/3]);
        temp = temp';
        B31 = [B31;temp];
        temp = [];
    end 
    
    if contains(line2,'B32')
        temp = fscanf(fid2,'%d,%f,%f,%f\n');
        temp = reshape(temp,[4 length(temp)/4]);
        temp = temp';
        B32 = [B32;temp];
        temp = [];
    end 
    
    line2 = fgetl(fid2);
    
    if isempty(line2)
        line2=fgetl(fid2);
    end 

end 

%B32 = B32(~all(B32 == 0, 2),:); % Removes (0,0) entries from T3D2 

temp = C3D8R; 
for i=1:size(temp,1)
    C3D8R(temp(i,1),:) = temp(i,:);
end 

% Take care of beams
if ~isempty(B31)
    BeamOrder = 1;
    BeamCon = B31;
else
    BeamOrder = 2;
    BeamCon = B32;
end
xis = linspace(-1,1,nplotbeams); % for visualising polylines
Phis = ShapeFnc(xis', BeamOrder);

nS = length(unique(C3D8R(:,2:end)));
nB = length(unique(BeamCon(:,2:end)));

[xGP, wGP] = lgwt(BeamOrder,-1,1);
[NGP,dNGP] = ShapeFnc(xGP, BeamOrder);

%% WRITE VTK Steps

n_steps = size(AbaqusData.U,3);
count = -1;
BeamVisNodes = [];

%Loop steps
for i = 1:n_steps
    
        count = count+1;
        Utemp=[];
        Utemp(:,:) = AbaqusData.U(:,:,i);
        Stemp(:,:) = AbaqusData.S(:,:,i);
        if max(Stemp(:,2))>1
            error('Solid Stress valid only for Reduced Integration Elements');
        end
        
        dS = Utemp(1:nS,2:4);
        dBr = Utemp(nS+1:nS+nB,2:4);
        dB = [dBr]; % Add zeros: rotations (only contribute at the Euler-Bernoulli coupling which is not implemented yet)
        %Flatten the global nodal values vectors
        dS = reshape(dS',3*size(dS,1),1);
        dB = reshape(dB',3*size(dB,1),1);
        
        if LM_bool
        LM_force = e_pen*kinv_ns*[-M_ns,D_ns]*[dS;dB];
        LM_force = reshape(LM_force,3,nB)'; % Reshape for (x,y,z) values at each Lagrange multiplier node
        
        Gc = [-M_ns,D_ns]*[dS;dB]; % Relative beam-solid displacement
        Gc = reshape(Gc,3,nB)'; 
        end
        
        % Solid stress values
        sqmat = [];
        for ii = 1:size(Stemp,1)
            sqmat = [sqmat;diag(Stemp(ii,3:5)) + squareform(Stemp(ii,6:8))]; 
        end
                       
        
        %% Solid Elements
        file_ext = ['_',num2str(count,'%.5d'),'.vtk'];
        fid3 = fopen([solid_file,file_ext],'w');
        fprintf(fid3,'# vtk DataFile Version 2.0\n');
        fprintf(fid3,'test\n');
        fprintf(fid3,'ASCII\n');
        fprintf(fid3,'DATASET UNSTRUCTURED_GRID\n');

        fprintf(fid3,['POINTS ',num2str(size(P,1)),' float\n']);
        fprintf(fid3,'%.7e %.7e %.7e\n',(P(:,2:4)+ Utemp(:,2:4))');

        fprintf(fid3,['CELLS ',num2str(size(C3D8R,1)),' ',num2str(9*size(C3D8R,1)),'\n']);
        fprintf(fid3,'8 %d %d %d %d %d %d %d %d\n',C3D8R(:,2:9)'-1);

        fprintf(fid3,['CELL_TYPES ',num2str(size(C3D8R,1)),'\n']);
        fprintf(fid3,'%d\n',12*ones(1,size(C3D8R,1))); 
        
        fprintf(fid3,['CELL_DATA ',num2str(size(C3D8R,1)),'\n']);
        fprintf(fid3,'TENSORS S float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',sqmat');
        
         
        fprintf(fid3,['POINT_DATA ',num2str(size(P,1)),'\n']);
        fprintf(fid3,'VECTORS U float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',Utemp(:,2:4)');


        fclose(fid3);

         
        %% Beam Elements (Polylines)
        
        % Step one create new node matrix
        X = P(:,2:4)+Utemp(:,2:4);
        Xele = [];
        Uele = [];
        LMele = [];
        Gcele = [];
        SegmsEle = [];
        ElemEle = [];
        
        SFele = [];
        SEele = [];
        SKele = [];
        SMele = [];
        Lamele = [];
        
        for iii = 1:size(BeamCon,1) %loop beam elements
            locIDs = BeamCon(iii,2:end);
            Xele = [Xele; Phis*X(locIDs,:)];
            Uele = [Uele; Phis*Utemp(locIDs,2:4)];
            
            if i==1
                % Segment Data
                SegmentMat = GlobalBeamMesh.SegmentMat{1, iii}; 
                nSegms = size(SegmentMat,1);
                xiSegms = unique(SegmentMat(:,1:2));
                xis = xiSegms(xiSegms~=1 & xiSegms~=-1); %Exclude element endpoints
                if ~isempty(xis)
                    [NSegms,~] = ShapeFnc(xis, BeamOrder);
                    SegmsEle = [SegmsEle; NSegms*X(locIDs,:)];
                end
                % Repeat for element nodes
                if BeamOrder==2
                    [NElems,~] = ShapeFnc([-1,0,1], BeamOrder);
                elseif BeamOrder==1
                    [NElems,~] = ShapeFnc([-1,1], BeamOrder);
                end
                    ElemEle = [ElemEle; NElems*X(locIDs,:)];
            end
            
            
            if LM_bool
                LMele = [LMele; Phis*LM_force(locIDs-nS,:)];
                Gcele = [Gcele; Phis*Gc(locIDs-nS,:)]; 
            end
            
            % Compute Averaged Element SF, SE, SK, SM (Note: this not very
            % accurate: Averaging flux on Gauss points is questionable)
            ID0 = find(BeamCon(iii,1)==AbaqusData.SF(:,1,i));
            SFtemp = squeeze(AbaqusData.SF(ID0(1):ID0(1)+BeamOrder-1,3:end,i));
            SEtemp = squeeze(AbaqusData.SE(ID0(1):ID0(1)+BeamOrder-1,3:end,i));
            SMtemp = squeeze(AbaqusData.SM(ID0(1):ID0(1)+BeamOrder-1,3:end,i));
            SKtemp = squeeze(AbaqusData.SK(ID0(1):ID0(1)+BeamOrder-1,3:end,i));
            Lam1(:,1) = exp(SEtemp(:,1));
            Lam1(:,2) = 1;
            Lam1(:,3) = 1;
            
            SFele = [SFele; repmat(mean(SFtemp),[nplotbeams,1])];
            SEele = [SEele; repmat(mean(SEtemp),[nplotbeams,1])];
            SMele = [SMele; repmat(mean(SMtemp),[nplotbeams,1])];
            SKele = [SKele; repmat(mean(SKtemp),[nplotbeams,1])];
            Lamele = [Lamele; repmat(mean(Lam1), [nplotbeams,1])];
            
            % Compute Element Strain Energy
            UmL(i,iii) = 0; % Umebrane, linear
            UmNL(i,iii) = 0; % Umebrane, linear
            Ub1L(i,iii) = 0; % U bending axis1, linear
            Ub2L(i,iii) = 0; % U bending axis2, linear
            Us1L(i,iii) = 0; % U shear axis1, linear
            Us2L(i,iii) = 0; % U shear axis2, linear
            UtrL(i,iii) = 0; % U shear axis1, linear
            
            Jel = sqrt(sum((dNGP*P(locIDs,2:4)).^2,2));
            for iv = 1:BeamOrder
                I1 = Lam1(iv,1)^2+2*Lam1(iv,1)^(-1);
                
                UmL(i,iii) = UmL(i,iii) + 1/2*wGP(iv)*SFtemp(iv,1)*SEtemp(iv,1)*Jel(iv)*Lam1(iv,1);
                UmNL(i,iii) = UmNL(i,iii) + 1/2*BCS.A*BCS.mu*wGP(iv)*(I1-3)*Jel(iv)*Lam1(iv,1);
                
                Ub1L(i,iii) = Ub1L(i,iii)+ 1/2*wGP(iv)*SKtemp(iv,1)*SMtemp(iv,1)*Jel(iv)*Lam1(iv,1);
                Ub2L(i,iii) = Ub2L(i,iii)+ 1/2*wGP(iv)*SKtemp(iv,2)*SMtemp(iv,2)*Jel(iv)*Lam1(iv,1);
                
                Us1L(i,iii) = Us1L(i,iii) + 1/2*wGP(iv)*SFtemp(iv,2)*SEtemp(iv,2)*Jel(iv)*Lam1(iv,1);
                Us2L(i,iii) = Us2L(i,iii) + 1/2*wGP(iv)*SFtemp(iv,3)*SEtemp(iv,3)*Jel(iv)*Lam1(iv,1);
                
                UtrL(i,iii) = UtrL(i,iii) + 1/2*wGP(iv)*SKtemp(iv,3)*SMtemp(iv,3)*Jel(iv)*Lam1(iv,1);
                
            end
            % Nonlinear counterparts: differ only by the initial modulus
            UtrNL(i,iii) = 4*UtrL(i,iii);
            Us1NL(i,iii) = Us1L(i,iii);
            Us2NL(i,iii) = Us2L(i,iii);
            
        end
        
        Energy.UmL(i) = sum(UmL(i,:));
        %Energy.UmNL(i) = sum(UmNL(i,:));
        Energy.Ub1L(i) = sum(Ub1L(i,:));
        Energy.Ub2L(i) = sum(Ub2L(i,:));
        Energy.Us1L(i) = sum(Us1L(i,:));
        Energy.Us2L(i) = sum(Us2L(i,:));
        %Energy.Us1NL(i) = sum(Us1NL(i,:));
        %Energy.Us2NL(i) = sum(Us2NL(i,:));
        Energy.UtrL(i) = sum(UtrL(i,:));
        %Energy.UtrNL(i) = sum(UtrNL(i,:));
        Energy.Gc{i} = Gc;
        
        Xcon = reshape(1:size(BeamCon,1)*nplotbeams,nplotbeams, size(BeamCon,1))';     
                
        fid3 = fopen([beam_file,file_ext],'w'); 
        fprintf(fid3,'# vtk DataFile Version 2.0\n');
        fprintf(fid3,'test\n');
        fprintf(fid3,'ASCII\n');
        fprintf(fid3,'DATASET POLYDATA\n');

        fprintf(fid3,['POINTS ',num2str(size(Xele,1)),' float\n']);
        fprintf(fid3,'%.7e %.7e %.7e\n',Xele');

        fprintf(fid3,['LINES ',num2str(size(BeamCon,1)),' ',num2str((1+nplotbeams)*size(BeamCon,1)),'\n']);
        fprintf(fid3,[num2str(nplotbeams), ' ',repmat('%d ', 1, nplotbeams), '\n'],(Xcon-1)');
        
        fprintf(fid3,['POINT_DATA ',num2str(size(Xele,1)),'\n']);
        fprintf(fid3,'VECTORS U float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',Uele');
        
        if LM_bool
        fprintf(fid3,'VECTORS LM float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',LMele');
        
        fprintf(fid3,'VECTORS GC float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',Gcele');
        end
        
         fprintf(fid3,'VECTORS SF float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',SFele');
        
        fprintf(fid3,'VECTORS SE float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',SEele');
        
        fprintf(fid3,'VECTORS SM float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',SMele');
        
        fprintf(fid3,'VECTORS SK float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',SKele');
        
        fprintf(fid3,'VECTORS Lam1 float\n');
        fprintf(fid3,'%.7e %.7e %.7e\n',Lamele');

        fclose(fid3); 
        
        % Write SegmentPoints file
        if i==1
            % Segment Points (Glyphs)
            fid4 = fopen([seg_file,file_ext],'w'); 
            fprintf(fid4,'# vtk DataFile Version 2.0\n');
            fprintf(fid4,'test\n');
            fprintf(fid4,'ASCII\n');
            fprintf(fid4,'DATASET POLYDATA\n');

            fprintf(fid4,['POINTS ',num2str(size(SegmsEle,1)),' float\n']);
            fprintf(fid4,'%.7e %.7e %.7e\n',SegmsEle');
            fclose(fid4);
            
            % Element points (nodes)
            fid5 = fopen([elem_file,file_ext],'w'); 
            fprintf(fid5,'# vtk DataFile Version 2.0\n');
            fprintf(fid5,'test\n');
            fprintf(fid5,'ASCII\n');
            fprintf(fid5,'DATASET POLYDATA\n');

            fprintf(fid5,['POINTS ',num2str(size(ElemEle,1)),' float\n']);
            fprintf(fid5,'%.7e %.7e %.7e\n',ElemEle');
            fclose(fid5); 
        end
    end
end
