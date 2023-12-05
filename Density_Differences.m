 %%  %Compare results from several conditions output from the script:
%"Two_Parameter_SMT_V1"
%Pablo Aurelio Gomez Garcia, 2020.


clear all
close all
clc

warning off
%% Inputs
%------------------------------------------------------------------

%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)
startlocn = {'C:\Users\pgomez.CRG\Dropbox (CRG ADV)\Pablo Gomez\FoxA1 SMT\20200508 New\Dataset\Cell1'};

%Set the name of the REFERENCE
%------In our case if H2B (or any other chromatin constitutional protein)----------------
type.t1 = 'H2B'; %The first dataset must be H2B (or the reference to which we compare)
%-------------------------------------------------------------------------

%Set the name of the Transcription Factors (Chromating Binding proteins)
%----------- (Ex: Oct4, Klf4, Sox2, H1P, FOXA1 ...) -------------------------
type.t2 = 'TF'; %Here it goes the TFs to study
% type.t3 = ''; 
% type.t4 = '';


%Numerical Inputs
%----------------------------------------------------------
sample_size = []; %In order to compare the different conditions with respect to H2B, we randomly subsample the datasets to have the same numeber of tracks on every condition
%If sample_size = [], we will subsample to the condition with the lowest
%number of tracks.

Num_points_diff = 5; %In order to compare the density plots, for each point of the TF, 
%we take the average value of the closest "Num_points_diff" points of H2B.

%Limits for representation of the Two-Parameter scatter plot
x_lim = [0 80]; %Units are in nm
y_lim = [0 80]; %Units are in nm

%Define the limits for the different mobility regions of the Two-parameter scatter
%plot:

%Very low mobility region 
R_Conf_limits_vL = [20 30]; %Units are in nm
Avg_Jump_limits_vL = [15 29];
%Low mobility region
R_Conf_limits_L = [35 50]; %Units are in nm
Avg_Jump_limits_L = [20 35];
%Intermediate mobility region
R_Conf_limits_I = [15 37.5]; %Units are in nm
Avg_Jump_limits_I = [29 36];
%High mobility region
R_Conf_limits_H = [37.5 55]; %Units are in nm
Avg_Jump_limits_H = [37.5 60];
%Upper mobility region
R_Conf_limits_Upp = [55 300]; %Units are in nm
Avg_Jump_limits_Upp = [60 300];


%% Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
datatypes = fieldnames(type);
ntypes = size(datatypes,1);         


%% Load the Data
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%Ask the user to select the files
files = cell(ntypes,1);
for t = 1:ntypes;
    files{t}= Select1DataGroup(type.(datatypes{t}),'*.csv',startlocn{1});
end

%Load the .xlsx files for H2B and TF
% 
for t=1:ntypes;
        
    Two_parameter_SMT{t} = csvread(strcat(files{t}.data{1,2},files{t}.data{1,1}),1,0);

end


%% Subsample the datasets
if isempty(sample_size);
    for t=1:ntypes; 
    Num_Tracks(t) = size(Two_parameter_SMT{t},1); 
    end
sample_size = min(Num_Tracks(t));
end

for t=1:ntypes;  
    
%Column 1: radius of confinement (µm). Colunm 2: average jump (µm).
Two_parameter_SMT_SubSampled{t} = datasample(Two_parameter_SMT{t}(:,[1 2]), sample_size, 'Replace',false);

end


%% Scatter plot of Radius of Confinement vs Avg frame-to-frame Jump for SubSampled datasets

for t=1:ntypes;  
figure()
out{t} = scatplot(Two_parameter_SMT_SubSampled{t}(:,1),Two_parameter_SMT_SubSampled{t}(:,2),'circles',3.2,100,5,3,10);
title(strcat(files{t}.dataname(1:end-9),' - Scatter density plot (downsampled)'),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
xlabel('Radius of confinement (nm)');
ylabel('Average jump (nm)');
axis('square');
% legend(strcat(files{t}.dataname(1:end-9)));
set(gca,'FontSize',20,'FontWeight','bold')
xlim(x_lim);
ylim(y_lim);
end


%% Calculate density differences using H2B as a reference (Make sure that type.t1 is H2B)

Reference_density_H2B = (out{1}.dd)/max(out{1}.dd); 
Reference_position_H2B = [Two_parameter_SMT_SubSampled{1}(:,1) Two_parameter_SMT_SubSampled{1}(:,2)];

strB = '';
for t=2:ntypes;
    
TF_density{t} = (out{t}.dd)/max(out{t}.dd);
TF_position{t} = [Two_parameter_SMT_SubSampled{t}(:,1) Two_parameter_SMT_SubSampled{t}(:,2)];    

    for i=1:size(TF_position{t},1);
        dist_idx{t}{i} = knnsearch(Reference_position_H2B,TF_position{t}(i,:),'K',Num_points_diff);        
        Difference_percentage{t}(i) = (mean(-Reference_density_H2B(dist_idx{t}{i})) + TF_density{t}(i))/mean(Reference_density_H2B(dist_idx{t}{i}));   
    end

strOut = sprintf('Calculating differences in density plots: % 4.1f',100*t/(ntypes));
fprintf([strB strOut '%%\n']);
strB = repmat('\b',1,length(strOut)+2);

end


%% Compare the density scatter plots to H2B, using H2B as a reference
for t=2:ntypes;
figure()
map = colormap(jet);
ind = fix((Difference_percentage{t}-min(Difference_percentage{t}))/(max(Difference_percentage{t})-min(Difference_percentage{t}))*(size(map,1)-1))+1;
x = Two_parameter_SMT_SubSampled{t}(:,1);
y = Two_parameter_SMT_SubSampled{t}(:,2);
for k = 1:size(map,1);
    if any(ind==k);
   line('Xdata',x(ind==k),'Ydata',y(ind==k),'LineStyle','none','Color',map(k,:),'Marker','.','MarkerSize',8);
    end
end

title(strcat('Density differences between',{' '}, files{1}.dataname(1:end-9),{' '}, '&',{' '},files{t}.dataname(1:end-9)),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
xlabel('Radius of confinement (nm)');
ylabel('Average jump (nm)');
axis('square');
% legend(strcat(files{t}.dataname(1:end-9)));
set(gca,'FontSize',20,'FontWeight','bold')
colorbar
xlim(x_lim);
ylim(y_lim);
% caxis([0 5]);

figure()
z = (Difference_percentage{t})';
scatter3(x,y,z,20,z,'filled');
title(strcat('Density differences between',{' '}, files{1}.dataname(1:end-9),{' '}, '&',{' '},files{t}.dataname(1:end-9)),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
xlabel('Radius of confinement (nm)');
ylabel('Average jump (nm)');
zlabel('Density differences');
axis('square');
% legend(strcat(files{t}.dataname(1:end-9)));
set(gca,'FontSize',20,'FontWeight','bold')
% caxis([-1 1]);
xlim(x_lim);
ylim(y_lim);
colorbar

end


%% PLOT ERRORBAR WITH THE DENSITY DIFFERENCES
for t=2:ntypes;

X{t} = [TF_position{t} Difference_percentage{t}'];
x=X{t}(:,1);
y=X{t}(:,2);
z=X{t}(:,3);

%Divide the tracks into 5 populations based on their positions on the
%two-parameter scatter plot
POP1{t} = find((x>R_Conf_limits_vL(1) & x<R_Conf_limits_vL(2) )& (y>Avg_Jump_limits_vL(1) & y<Avg_Jump_limits_vL(2)));
POP2{t} = find((x>R_Conf_limits_L(1) & x<R_Conf_limits_L(2) )& (y>Avg_Jump_limits_L(1) & y<Avg_Jump_limits_L(2)));
POP3{t} = find((x>R_Conf_limits_I(1) & x<R_Conf_limits_I(2) )& (y>Avg_Jump_limits_I(1) & y<Avg_Jump_limits_I(2)));
POP4{t} = find((x>R_Conf_limits_H(1) & x<R_Conf_limits_H(2) )& (y>Avg_Jump_limits_H(1) & y<Avg_Jump_limits_H(2)));
POPSUP{t} = find((x>R_Conf_limits_Upp(1) & x<R_Conf_limits_Upp(2)) |  (y>Avg_Jump_limits_Upp(1) & y<Avg_Jump_limits_Upp(2)));

D1{t} = X{t}(POP1{t},3);
D2{t} = X{t}(POP2{t},3);
D3{t} = X{t}(POP3{t},3);
D4{t} = X{t}(POP4{t},3);
DSUP{t} = X{t}(POPSUP{t},3);

SEM1{t} = std(D1{t})/sqrt(length(D1{t}));
SEM2{t} = std(D2{t})/sqrt(length(D2{t}));
SEM3{t} = std(D3{t})/sqrt(length(D3{t}));
SEM4{t} = std(D4{t})/sqrt(length(D4{t}));
SEMSUP{t} = std(DSUP{t})/sqrt(length(DSUP{t}));

AVE1{t} = -mean(D1{t});
AVE2{t} = -mean(D2{t});
AVE3{t} = -mean(D3{t});
AVE4{t} = -mean(D4{t});
AVESUP{t} = -mean(DSUP{t});

%Plot
data=[AVE1{t} AVE2{t} AVE3{t} AVE4{t}];
errhigh=[SEM1{t} SEM2{t} SEM3{t} SEM4{t}];
errlow=[SEM1{t} SEM2{t} SEM3{t} SEM4{t}];
figure()
bar(1:4,data);                
hold on
er = errorbar(1:4,data,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title(strcat('Differences between H2B & ',strcat(files{t}.dataname(1:end-9)), ' for the different regions'),'FontSize',24,'FontWeight','bold','Interpreter','None'); 
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'Region vL','Region L','Region I','Region H'});
set(gca,'FontSize',20,'FontWeight','bold')

end

%% Save the results in tables

cHeader = {'Radius Conf (nm)' 'Avg Jump (nm)' 'Density difference relative to H2B'};
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

for t=2:ntypes;
mkdir(strcat(files{t}.data{1,2},'Density differences'));
fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences.csv'),X{t},'-append');


fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_vL_Region.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_vL_Region.csv'),X{t}(POP1{t},:),'-append');


fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_L_Region.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_L_Region.csv'),X{t}(POP2{t},:),'-append');


fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_I_Region.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_I_Region.csv'),X{t}(POP3{t},:),'-append');


fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_H_Region.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_H_Region.csv'),X{t}(POP4{t},:),'-append');


fid = fopen(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_SUP_Region.csv'),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(strcat(files{t}.data{1,2},'Density differences'),'Density_Differences_SUP_Region.csv'),X{t}(POPSUP{t},:),'-append');


end








