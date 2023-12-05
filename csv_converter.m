   function [] = csv_converter()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

d=dir('*.txt');

for i=1:size(d,1),
    a=dlmread(d(i).name);
    c=[(0:1:(size(a,1)-1))',a(:,3),a(:,3)*0.01,a(:,4),a(:,1)*0.116,a(:,2)*0.116];
    %[c]=column_reorder(a);
        
    name=strsplit(d(i).name,'_t');
    
    fileID = fopen(strcat(name{1},name{2},'spoton','.csv'), 'a');
    fprintf(fileID, ',%s,', 'frame');
    fprintf(fileID, '%s,', 't');
    fprintf(fileID, '%s,', 'trajectory');
    fprintf(fileID, '%s,', 'x');
    fprintf(fileID, '%s,', 'y');
    fprintf(fileID, '\n');
    for j=1:size(a,1),
        fprintf(fileID, '%f,%f,%f,%f,%f,%f \n', c(j,:));
       end
    
    
    fclose(fileID);
   
end


end
