   %%  %Compare results from several conditions output from the script:
% "SMT_Motion_Classifier_V1.m"
%Pablo Aurelio Gomez Garcia, 2020.

clear all
close all
clc
warning off
set(0,'DefaultFigureColormap',feval('jet'));

%% Inputs
%------------------------------------------------------------------
%Data---------------------------------------------------
% %------------------------------------------------------------------
%Set the data type description
type.t1 = 'DATA';
% type.t2 = 'ALL';

%Set the folder as a cell array. (If only one
% location is included, then it will be used for all data types)
startlocn = {'C:\Users\pgomez.CRG\Dropbox (CRG ADV)\Pablo Gomez\FoxA1 SMT\20200508 New\Dataset\Cell1'};

%------------------------------------------------------------------
%Numerical Inputs---------------------------------------------------
%------------------------------------------------------------------
n_dim = 2; %Dimensionality of the movement (2 for 2D and 3 for 3D).
Frame_interval = 0.010; %Exposure time in seconds

%-----------------------------------------------------------------
%Confined Circle Diffusion Model Fitting---------------------------------------------------
%---------------------------------------------------------------------
num_points = 10; %Number of points for fitting the confined diffussion circle model to the TE-MSD and T-MSD.
offset0 = 0.001; %This is the initial value for the offset of the MSD curve. It depens on the localization precision. offset = 4*(Loc_precition^2) [in µm]. 
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (µm^2/s)
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (µm)

%-----------------------------------------------------------------
%Two-Paramater SMT Analysis and representation ---------------------------------------------------
%--------------------------------------------------------------------
Resnorm_max = 1e-5; %This is used to filter the good fits. The maximum value residuals normalized.
Max_RConf = 300; %Maximum value of radius of confinement (nm)
%Limits for representation of the Two-Parameter scatter plot
x_lim = [0 80]; %Units are in nm
y_lim = [0 80]; %Units are in nm

%% Initialize variables-------------------------------------------------
%----------------------------------------------------------------------
Dcoeffs = {};
datatypes = fieldnames(type);
ntypes = size(datatypes,1);         
Color=lines(ntypes);



%% Load the Data
%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
%Ask the user to select the files
files = cell(ntypes,1);
for t = 1:ntypes;
    files{t}= Select1DataGroup(type.(datatypes{t}),'mat',startlocn{1});
end


%Load the .mat files with the Results of the Tracking Analysis
%You can load multiple files per category.
for t=1:ntypes;

MSD_Results{t}.ma = msdanalyzer(2, 'µm', 's');  %Initialize MSD analizer
  
    for c=1:size(files{t}.data,1);
        
    MSD_Results_prov{t}{c} = load(strcat(files{t}.data{c,2},files{t}.data{c,1}));
    
    MSD_Results{t}.ma = MSD_Results{t}.ma.addAll(MSD_Results_prov{t}{c}.ma.tracks);

    end
    
    MSD_Results{t}.ma = MSD_Results{t}.ma.computeMSD;

end


%% Calculate some numbers
for t=1:ntypes;
    
Num_Tracks(t) = size(MSD_Results{t}.ma.tracks,1);

for j=1:Num_Tracks(t);
Track_Length{t}(j) = size(MSD_Results{t}.ma.tracks{j},1);
end

end


%% - Calculate the geometrical distances of each track
for t=1:ntypes;
for ggg = 1:size(MSD_Results{t}.ma.tracks,1);
points_coord = MSD_Results{t}.ma.tracks{ggg}(:,2:3)*1000; %Units are in nm
[max_dist{t}(ggg), min_dist{t}(ggg), avg_dist{t}(ggg)] = distance_scatter(points_coord);
end
end



%% Fitting the confined diffusion model to each trajectory (T-MSD)-----------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
strB = '';
for t=1:ntypes;

for i=1:size(MSD_Results{t}.ma.msd,1);
    
local_msd = MSD_Results{t}.ma.msd{i};
nonnan = ~isnan(local_msd(:,2));  
local_msd = local_msd(nonnan,1:2);

[Val_fitted_all{t}{i}, residual_all{t}{i},jacobian_all{t}{i},resnorm_all{t}{i}] = confined_diffusion_fit_V2(local_msd,offset0,num_points,D0,R0);
conf_all{t}{i} = nlparci(Val_fitted_all{t}{i},residual_all{t}{i},'jacobian',jacobian_all{t}{i});
f_confined_circ_diff_D_all{t}(i) = Val_fitted_all{t}{i}(2);
f_confined_circ_diff_R_all{t}(i) = abs(Val_fitted_all{t}{i}(1))*1000; %Units are in nm
f_confined_circ_diff_offset_all{t}(i) = Val_fitted_all{t}{i}(3);
Localization_precision_estimation_all{t}(i) = sqrt(f_confined_circ_diff_offset_all{t}(i)/(n_dim^2));


  strOut = sprintf('Completed: % 4.1f',100*i/(size(MSD_Results{t}.ma.msd,1)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
end

fprintf('\n')

%Filter. Just take the good fits. 
iidx{t} = find(cell2mat(resnorm_all{t}) < Resnorm_max);
f_confined_circ_diff_R_Filtered{t} = f_confined_circ_diff_R_all{t}(iidx{t});
avg_dist_Filtered{t} = avg_dist{t}(iidx{t});
Track_Length_Filtered{t} = Track_Length{t}(iidx{t});


end

 
%% Save the results from the two-parameter SMT Analysis
for t=1:ntypes;

RC1 = (f_confined_circ_diff_R_all{t})';
C = (avg_dist{t})';
Tr = (Track_Length{t})';
TABLEFINAL{t} = [RC1 C Tr];

RC1 = (f_confined_circ_diff_R_Filtered{t})';
C = (avg_dist_Filtered{t})';
Tr = (Track_Length_Filtered{t})';
TABLEFINAL_FILTERED{t} = [RC1 C Tr];

cHeader = {'Radius Conf (nm)' 'Avg Jump (nm)' 'TrackLength (frames)'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

fid = fopen(fullfile(files{t}.data{1,2},strcat(files{t}.dataname(1:end-9),'_SMT_TwoParameters.csv')),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(files{t}.data{1,2},strcat(files{t}.dataname(1:end-9),'_SMT_TwoParameters.csv')),TABLEFINAL{t},'-append');

fid = fopen(fullfile(files{t}.data{1,2},strcat(files{t}.dataname(1:end-9),'_SMT_TwoParameters_Filtered.csv')),'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite(fullfile(files{t}.data{1,2},strcat(files{t}.dataname(1:end-9),'_SMT_TwoParameters_Filtered.csv')),TABLEFINAL_FILTERED{t},'-append');

end


%% FIGURES: 

for t=1:ntypes;
sample_sizes(t) = length(f_confined_circ_diff_R_all{t});
end

%SCATTER PLOT
for t=1:ntypes;
figure()
A = datasample(TABLEFINAL{t}, sample_sizes(t), 'Replace',false);
xa = A(:,1);
ya = A(:,2);
filta = find(xa<Max_RConf); %Filtering tracks with Radius of confinement bigger than the input threshold in µm.
x = A(filta,1);
y = A(filta,2);
scatter(x,y);
xlim(x_lim);
ylim(y_lim);
title('Scatter plot of average jump vs radius of confinement');
xlabel('Radius of confinement (nm)');
ylabel('Average jump (nm)');
caxis([0 1]);
% legend(strcat(files{t}.dataname(1:end-9)));
set(gca,'FontSize',20,'FontWeight','bold');
end



%TWO PARAMETER DENSITY SCATTER PLOT
for t=1:ntypes;
figure()
A = datasample(TABLEFINAL{t}, sample_sizes(t), 'Replace',false);
xa = A(:,1);
ya = A(:,2);
filta = find(xa<Max_RConf); %Filtering tracks with Radius of confinement bigger than the input threshold in µm.
x = A(filta,1);
y = A(filta,2);
scatplot(x,y,'circles',3.2,100,5,3,15);
xlim(x_lim);
ylim(y_lim);
title('Scatter density plot of average jump vs radius of confinement');
xlabel('Radius of confinement (nm)');
ylabel('Average jump (nm)');
caxis([0 1]);
% legend(strcat(files{t}.dataname(1:end-9)));
set(gca,'FontSize',20,'FontWeight','bold');

end


