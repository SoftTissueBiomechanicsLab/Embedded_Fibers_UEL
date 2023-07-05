function J = ReadDeterminant(cur_path, JobNum, nelems, n_incr, RPL)

timesteps = 0:1/n_incr:1;
timesteps(1) = []; % Does not export first step.

IDStr = num2str(JobNum,'%.4d');
full_workdir = [cur_path,'/AbaqusWorkDir/Job',IDStr,'/'];
dat_file = [full_workdir,'Model',IDStr,RPL,'.dat'];


%% READ DISPLACEMENT ODB DATA 
fid1 = fopen(dat_file,'r');

line = fgetl(fid1);
J = nan(length(timesteps),nelems); % initialize
endflag = false;
step_count = 0;

while ~endflag
     
    if contains(line,'STEP TIME COMPLETED ')
        time_inc = sscanf(line,' STEP TIME COMPLETED %f\n');
        while ~contains(line,'DG11')
             line=fgetl(fid1); % Skip "empty" lines.
        end
        
        if min(abs(time_inc-timesteps))<1e-14 % Check if current time is in mustpoints
            
            step_count = step_count+1;
            line=fgetl(fid1);
            line=fgetl(fid1);
            temp = fscanf(fid1,'\t%d %d %f %f %f %f %f %f %f %f %f\n');
            temp = reshape(temp,[11 length(temp)/11])'; 
            for ii = 1:size(temp,1)
                vec = temp(ii,3:end);
                sqmat = diag(vec(1:3)) + [0,vec(4),vec(5);vec(7),0,vec(6);vec(8),vec(9),0];
                J(step_count,ii) = det(sqmat);
            end
        end
    end
    
  


  line=fgetl(fid1);
  
  if contains(line,'ANALYSIS COMPLETE')
      endflag = true;
  end

  
    
end                                                                                                                                     
fclose(fid1);

