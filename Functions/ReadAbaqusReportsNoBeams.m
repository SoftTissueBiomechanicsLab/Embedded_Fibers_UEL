function AbaqusData = ReadAbaqusReportsNoBeams(cur_path, JobNum)
% Reads Abaqus Reports and generates matlab data structures

% Define paths
IDStr = num2str(JobNum,'%.4d');
full_workdir = [cur_path,'/AbaqusWorkDir/Job',IDStr,'/'];
Field_file = [full_workdir,'Model',IDStr,'_Fields.csv'];
History_file = [full_workdir,'Model',IDStr,'_History.csv'];

%% READ DISPLACEMENT ODB DATA 
fid1 = fopen(Field_file,'r');

line = fgetl(fid1);
flag2=0;
count=0;
U = []; % U(n_nodes,[ID, U1, U2, U3], StepNumber, FrameNumber)
S = [];
SE =[];
SF = [];
SM = [];
SK = [];
RF = [];
RF_disp = [];

frame = 0;
Totframe = 1;
step_count = 0;
frame_count = 0;
while line>=0
    
    if contains(line,'Total number of steps = ')
        n_steps = strsplit(line, '=');
        n_steps = str2double(n_steps{2});
    end
    
    if contains(line, 'Step name')
        step_count = step_count+1;
    end
    
    if contains(line,'Total number of frames = ')
        Totframe = strsplit(line, '=');
        Totframe = str2double(Totframe{2});
    end
    if contains(line, 'Frame number')
        frame = strsplit(line, 'number');
        frame = str2double(frame{2});
        if frame==0
            frame=1; %sim has failed
        end
    end
    
    
    % DISPLACEMENTS
    if contains(line,'Instance,  Element ,    Node  , SP, IP,       U1     ,       U2     ,       U3')
        frame_count = frame_count+1;
        stpfrm = (step_count-1)*Totframe + frame_count;
        temp = fscanf(fid1,'\tPART-1-1,\t,%d,\t,\t,%f,%f,%f\n'); 
        temp = reshape(temp,[4 length(temp)/4]); 
        temp = temp';
        temp = sortrows(temp,1);
        tempU(:,:,stpfrm) = temp; 

        for i=1:size(tempU,1)

            U(tempU(i,1),:,stpfrm) = tempU(i,:,stpfrm);

        end
        RF_disp(stpfrm,:) = tempU(end,2:end,stpfrm); % It's alwawys the last node

    end
           
    % STRESSES
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      S11     ,      S22     ,      S33     ,      S12     ,      S13     ,      S23     ')
        temp = fscanf(fid1,'\tPART-1-1,%d,\t,\t,%d,%f,%f,%f,%f,%f,%f\n'); 
        temp = reshape(temp,[8 length(temp)/8]); 
        temp = temp';
        temp = sortrows(temp,1);
        tempS(:,:,stpfrm) = temp; 

        for i=1:size(tempS,1)

            S(tempS(i,1),:,stpfrm) = tempS(i,:,stpfrm);

        end

    end
    
%     % BEAM STRESS COMPONENTS
%     if contains(line,'Instance,  Element ,    Node  , SP, IP,      S11     ,      S22     ,      S12')
%         temp = fscanf(fid1,'\tPART-1-1,%d,\t,%d,%d,%f,%f,%f\n'); 
%         temp = reshape(temp,[6 length(temp)/6]); 
%         temp = temp';
%         temp = sortrows(temp,1);
%         SB(:,:,stpfrm) = temp; 
% 
%     end

    % SECTION Strains
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      SE1     ,      SE2     ,      SE3')
        temp = fscanf(fid1,'\tPART-1-1,%d,\t,\t,%d,%f,%f,%f\n'); 
        temp = reshape(temp,[5 length(temp)/5]); 
        temp = temp';
        temp = sortrows(temp,1);
        SE(:,:,stpfrm) = temp; 

    end
    
    % SECTION FORCES
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      SF1     ,      SF2     ,      SF3')
        temp = fscanf(fid1,'\tPART-1-1,%d,\t,\t,%d,%f,%f,%f\n'); 
        temp = reshape(temp,[5 length(temp)/5]); 
        temp = temp';
        temp = sortrows(temp,1);
        SF(:,:,stpfrm) = temp; 

    end
    
    
    % SECTION MOMENTS
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      SM1     ,      SM2     ,      SM3')
        temp = fscanf(fid1,'\tPART-1-1,%d,\t,\t,%d,%f,%f,%f\n'); 
        temp = reshape(temp,[5 length(temp)/5]); 
        temp = temp';
        temp = sortrows(temp,1);
        SM(:,:,stpfrm) = temp; 

    end
    
    % SECTION CURVATURE
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      SK1     ,      SK2     ,      SK3')
        temp = fscanf(fid1,'\tPART-1-1,%d,\t,\t,%d,%f,%f,%f\n'); 
        temp = reshape(temp,[5 length(temp)/5]); 
        temp = temp';
        temp = sortrows(temp,1);
        SK(:,:,stpfrm) = temp; 

    end
    
    % RIGID FORCE
    if contains(line,'Instance,  Element ,    Node  , SP, IP,      RF1     ,      RF2     ,      RF3')
        temp = fscanf(fid1,'\tPART-1-1,\t,%d,\t,\t,%f,%f,%f\n'); 
        temp = reshape(temp,[4 length(temp)/4]); 
        temp = temp';
        temp = sortrows(temp,1);
        tempRF(:,:,stpfrm) = temp; 
        RF(stpfrm,:) = tempRF(end,2:end,stpfrm); % It's alwawys the last node          

    end

    
    

  line=fgetl(fid1);
    
  if isempty(line)
    line=fgetl(fid1);
  end 

  if isempty(line)
    line=fgetl(fid1);
  end 
end                                                                                                                                     
fclose(fid1);

%% READ HISTORY ODB DATA 
fid1 = fopen(History_file,'r');

line = fgetl(fid1);
flag2=0;
count=0;
ALLSE=[];
ALLSD = [];
ALLIE = [];
ALLWK = [];
frame = 0;
Totframe = 1;
step_count = 0;
frame_count = 0;
while line>=0
         
    % ALLSE
    if contains(line,'History Output ''ALLSE''')
        line=fgetl(fid1);
        line=fgetl(fid1);
        line=fgetl(fid1);
        temp = fscanf(fid1,'\t  %f, %f\n'); 
        temp = reshape(temp,[2 length(temp)/2]); 
        temp = temp';
        temp = sortrows(temp,1);
        ALLSE = temp; 

    end
    
     % ALLWK
    if contains(line,'History Output ''ALLWK''')
        line=fgetl(fid1);
        line=fgetl(fid1);
        line=fgetl(fid1);
        temp = fscanf(fid1,'\t  %f, %f\n'); 
        temp = reshape(temp,[2 length(temp)/2]); 
        temp = temp';
        temp = sortrows(temp,1);
        ALLWK = temp; 

    end
    
    if contains(line,'History Output ''ALLSD''')
        line=fgetl(fid1);
        line=fgetl(fid1);
        line=fgetl(fid1);
        temp = fscanf(fid1,'\t  %f, %f\n'); 
        temp = reshape(temp,[2 length(temp)/2]); 
        temp = temp';
        temp = sortrows(temp,1);
        ALLSD = temp; 

    end
    
    
     if contains(line,'History Output ''ALLIE''')
        line=fgetl(fid1);
        line=fgetl(fid1);
        line=fgetl(fid1);
        temp = fscanf(fid1,'\t  %f, %f\n'); 
        temp = reshape(temp,[2 length(temp)/2]); 
        temp = temp';
        temp = sortrows(temp,1);
        ALLIE = temp; 

    end
   

  line=fgetl(fid1);
    
  if isempty(line)
    line=fgetl(fid1);
  end 

  if isempty(line)
    line=fgetl(fid1);
  end 
end                                                                                                                                     
fclose(fid1);

%% EXPORT

AbaqusData.U = U;
AbaqusData.RFdisp = RF_disp;
AbaqusData.RF = RF;
AbaqusData.S = S;
AbaqusData.SE = SE;
AbaqusData.SF = SF;
AbaqusData.SM = SM;
AbaqusData.SK = SK;
AbaqusData.ALLSE = ALLSE;
AbaqusData.ALLIE = ALLIE;
AbaqusData.ALLSD = ALLSD;
AbaqusData.ALLWK = ALLWK;