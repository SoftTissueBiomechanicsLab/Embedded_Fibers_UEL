function [] = WriteFortranMatrix(data, matpath, matname, jobnum)

n_rows = size(data,1);
n_cols = size(data,2);

n_4rows = fix(n_cols/4);
n_mod = mod(n_cols, 4);

str_form = '%20.14E';
% every_line = [str_form ', ' str_form ', ' str_form ', ' str_form ' \n'];
% last_line = [str_form, repmat([', ',str_form],1,n_mod-1), ' \n'];
%no comma
every_line = [str_form ' ' str_form ' ' str_form ' ' str_form ' \n'];
last_line = [str_form, repmat([' ',str_form],1,n_mod-1), ' \n'];

matfilename = [matpath,matname,num2str(jobnum,'%.4d'),'.txt'];
fid = fopen(matfilename,'w');
% fprintf(fid, '%i \n', n_rows);
% fprintf(fid, '%i \n', n_cols);

for i = 1:n_rows
    fprintf(fid, every_line, data(i,1:n_4rows*4));
    if n_mod~=0
        fprintf(fid, last_line , data(i,n_4rows*4+1:end));
    end
end
fclose(fid);

% [row col v] = find(data);
% fid = fopen(matfilename,'w');
% fprintf(fid, '%i %i %20.14E\n', [col row v]);
% fclose(fid);
%dlmwrite(matfilename,data, 'delimiter', ' ','precision', '%.7E')
%writematrix(data,matfilename,'delimiter', ' ','precision', '%.14E')