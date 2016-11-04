function true = write_matrix_txt(string,matrix,format)

% This function is to write matrix into txt file.
% Songpeng Zu
% INPUT:STRING of file name;MATRIX;FORMAT(%S/%U/%F...);
% 120610

%% Main function
fp = fopen(string,'w');
[row,column] = size(matrix);
writeformat1 = [format,'\t'];
writeformat2 = [format,'\n'];
if fp>0
    for i = 1:row
        for j = 1:column
            if j==column
                fprintf(fp,writeformat2,matrix(i,j));
            else
                fprintf(fp,writeformat1,matrix(i,j));
            end
        end
    end
    true = 1;
else
    true = 0;
end
fclose(fp);
