       
%Co  de to analyze tracks for two-parameter single molecule analysis step 1
% Pablo Aurelio Gomez Garcia, Matlab 2015a. 2020* 
clear all
close all
clc

warning off %This is to remove this warning: "NARGCHK will be removed in a future release. Use NARGINCHK or NARGOUTCHK instead". To keep the code compatible with Matlab 2011. 
 
   %% Some inputs:
%%----------------------------------------------------------
%Define starting directory and input units of the tracks----------------------------------------
%--------------------------------------------------------------
% directory_name = uigetdir('C:\Users\pgomez\Desktop\Temporal');
directory_name = 'C:\Users\pgomez.CRG\Dropbox (CRG ADV)\Pablo Gomez\FoxA1 SMT\20200508 New\Dataset';

TrackMate = 0; %If the track file is an .xml file from TrackMate.
SlimFast = 1; %If the track file is a .csv file from slimFAST. 

Input_space_units = 'micron'; %Specify the input spatial units of the Tracks (With Trackmate I usually use microns)
Input_time_units = 'second'; %Specify the input time units of the Tracks ('frames' or 'second') (With TrackMate I usually use Frames)

%----------------------------------------------------------------------------
%Some Numerical Inputs-------------------------------------------------------------
%-------------------------------------------------------------------------
dimension = 2; %Dimensionality of the movement. 2 for 2D SMT and 3 for 3D SMT.

minimum_track_length = 5;  %Minimum track length in frames to take the track into account

Jump_confined_threshold = 0.1; %Value used to separate confined trajectories into two groups [Average Jump (µm)]. 

Frame_interval = 0.01; %Exposure time in seconds
pixel_size = 0.116;  %Effective Pixel size in µm

TMSD_fitting_points = 4;   %Number of points used of the T-MSD for fitting the Diffussion Coefficient. (LINEAR FITTING) The minimum 3 points for being able to calculate the Confidence Intervals
TEMSD_fitting_points = 4;  %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (LINEAR FITTING) The minimum 3 points for being able to calculate the Confidence Intervals
TLOGLOG_fitting_points = 20; %Number of points used of the TE-MSD for fitting the Diffussion Coefficient. (POWER-LAW FITTING).

max_alpha_conf = 0.7; %Maximum alpha value for the power-fitting of the TE-MSD in order to consider Confined Motion
min_alpha_directed = 1; %Minimum alpha value for the power-fitting of the TE-MSD in order to consider Directed Motion
%Everything in the middle will be consider as Brownian Motion

R2LIMIT = 0.7; %Lower limit of the r-squared value for the fitting of the T-MSD and TE-MSD. The Tracks with r-squared lower than that won't be classified.

%------------------------------------------------------------------
%Butterfly Trajectories---------------------------------------------------
%------------------------------------------------------------------
jump_threshold = 1.5; %For identifying butterfly trajectories. One Track MUST jump more than its own (average_jump + jump_threshold*std_jump) to be considered a Butterfly track
minim_dist = 8; %For identifying butterfly trajectories. One Track MUST travell a total distance bigger than its own average_jump*minim_dist to be considered a Butterfly track
Conf2JumpThr = 0; %A Butterfly track need to have a Jump equal or bigger than the average radius of confinement of its confined segments multiplied by Conf2JumpThr.
Out_percentage = 30; % Minimum Percentage of points that a JUMP must have OUTSIDE the previous and posterior polygon (CONVEXHULL) to be considered an OUTER Segment. (Number from 0 to 100). 
N_sliding = 3; %Number of points to check linearity of segments.
P_min_linear = 0.8; %P minimum to consider that a segment is linear. Number between 0 and 1.
angleTH = 45; %Minimum Angle for considering a jump as "directed" in butterfly tracks (degrees) Direct(>angleTH) or non-direct(<angleTH)
Max_Jump_Confined_Butterfly = 0.18; %Maximum jump in (µm) than a confined segment of a butterfly track can have. (Use for track segmentation).
Min_num_points_butt_confined_segment = 4; %Minimum number of points that a confined segment of a butterfly track can have. (The tracks that doesn't fulfill this condition will be discarded)

%------------------------------------------------------------------
%Circle Confined Diffusion---------------------------------------------------
%------------------------------------------------------------------
Offset0 = 0.001; %This is the minimum MSD value. It depens on the localization precision. level = 4*(Loc_precition^2) [in µm]. 
num_points = 12; %Number of points for fitting the confined diffussion circle model to the TE-MSD.
D0 = 0.05; %Starting value for the least squares fitting for the Diffusion Coefficient (µm^2/s) on the Confined Circle Diffusion Model
R0 = 0.05; %Starting value for the least squares fitting for the Radius of Confinement (µm) on the Confined Circle Diffusion Model


  %% Initialize Variables----------------------------------------
%%----------------------------------------------------------
%--------------------------------------------------------------
SPACE_UNITS = 'µm'; %This is just for visualization purposes
TIME_UNITS = 's'; %This is just for visualization purposes
ma = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);  %Initialize MSD analyzer
ma_AllTracks = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS); 
ma_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_High_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_confined_Low_D = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_brownian = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_confined = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);
ma_butterfly_segments_directed = msdanalyzer(dimension, SPACE_UNITS, TIME_UNITS);


%% Load the Data:
%-------------------------------------------------------------------------
%Tracks-------------------------------------------------------------
%--------------------------------------------------------------------------
d = dir(directory_name);
isub = [d(:).isdir]; 
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

if TrackMate == 1;
%Extract the trajectories and other related information from the results of the tracking  
[list_names,directory_name] = uigetfile(fullfile(directory_name,'*.xml'),'select the files with SMT Tracks from TrackMate','MultiSelect','on');
%%Start importing the TrackMate data.
    if iscell(list_names);
        strB = '';
    for g=1:size(list_names,2);
    [trajectory{g}, metadata{g}] = importTrackMateTracks(fullfile(strcat(directory_name,'\'),list_names{g}),'clipZ',1);
    strOut = sprintf('Loading and extracting the Data: % 4.1f',100*g/(size(list_names,2)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
    end
     
    else
    [trajectory{1}, metadata] = importTrackMateTracks(fullfile(strcat(directory_name,'\'),list_names),'clipZ',1); 

    end
end


%Load data from SlimFast
if SlimFast == 1;
[list_names,directory_name] = uigetfile(fullfile(directory_name,'*.csv'),'select the files with SMT Tracks from slimFAST','MultiSelect','on');    
 
    if iscell(list_names);
        strB = '';
        for g=1:size(list_names,2);
        [trajectory{g}, tracklength{g}] = slimFast_ImportTracks(fullfile(strcat(directory_name,'\'),list_names{g}));  
        strOut = sprintf('Loading and extracting the Data: % 4.1f',100*g/(size(list_names,2)));
            fprintf([strB strOut '%%\n']);
            strB = repmat('\b',1,length(strOut)+2);
        end
    else
    [trajectory{1}, tracklength{1}] = slimFast_ImportTracks(fullfile(strcat(directory_name,'\'),list_names));  
    end  
end
 


%% Convert spatial units to microns if input units were in pixels.    
if Input_space_units == 'pixels';
    for i = 1:size(trajectory,2);
        for j=1:size(trajectory{i},1)
            trajectory{i}{j}(:,2:3) = trajectory{i}{j}(:,2:3)*pixel_size;    
        end
    end
end
   
%From now on, all spatial units are in um.


%% Filter The Tracks based on their TrackLength.
%Track Length
if TrackMate == 1;
for ff=1:size(trajectory,2);
   
     
    for u=1:size(trajectory{1,ff},1);   
    Track_length{ff}(u,1) = size(trajectory{1,ff}{u,1},1);   
    end
   
    indices{ff} = find(Track_length{ff}(:,1) >= (minimum_track_length));
    trajectory_filtered{1,ff} = trajectory{1,ff}([indices{ff}]);
        
  
    
    
end
end

%Filter tracks from SlimFast based on their tracklength
if SlimFast == 1;
    for ff=1:size(trajectory,2);
    
    indices_length{ff} = find(tracklength{ff} >= (minimum_track_length));
    trajectory_filtered{1,ff} = trajectory{1,ff}([indices_length{ff}]);
       
    end
end

trajectories = cat(1,trajectory_filtered{:});  %CHOOSE ALL HERE

%% Multiply the frame by the exposure time (in TrackMate, I track them without specifying the frame interval)
if Input_time_units == 'frames'
for ff=1:size(trajectories,1);
    trajectories{ff,1}(:,1)=trajectories{ff,1}(:,1)*Frame_interval;
end

trajectories = trajectories';
else
trajectories = trajectories';
end
%From now on, trajectories are in um and seconds.


%% Data Analysis
fprintf('\n');     
fprintf('Calculating and plotting......\n');     


%% ------ Create a separate variable for All the Trajectories -----------
ma_AllTracks = ma_AllTracks.addAll(trajectories);
ma_AllTracks = ma_AllTracks.computeMSD;
ma_AllTracks = ma_AllTracks.LogTMSD(TLOGLOG_fitting_points);
ma_AllTracks = ma_AllTracks.TMSD(TMSD_fitting_points);



%% Identify Butterfly Tracks --------------------------------------
%--------------------------------------------------------------------------------------------
%%1. Detect "butterfly" motion
[Butterfly_trajectories,Butterfly_track_segments,Butterfly_Reference_Track_segment, motion_type_segment,reference_of_Tracks_butterfly,Num_BigJumps,BigJump_idx]=identify_Butterfly_tracks_V4(trajectories,jump_threshold,minim_dist,TLOGLOG_fitting_points,R2LIMIT,min_alpha_directed,Offset0,num_points,D0,R0,Conf2JumpThr,Out_percentage,N_sliding,P_min_linear,angleTH,Max_Jump_Confined_Butterfly);
trajectories(reference_of_Tracks_butterfly)=[];  



if isempty(Butterfly_trajectories);
else 
Butterfly_trajectories_segments = horzcat(Butterfly_track_segments{:}); %Reshape cell array to colapse all the cells in the top level.   
motion_type_segment = horzcat(motion_type_segment{:}); %Reshape cell array to colapse all the cells in the top level.       
    
%2. Discard those Butterfly tracks with confined segments shorter than N points.
discard_idxxxx = [];
for i=1:length(Butterfly_trajectories);
idxxxx = find(Butterfly_Reference_Track_segment == i);
motion_type_segment_local = motion_type_segment(idxxxx);
idxxxx2 = find(motion_type_segment_local == 1);

if ~isempty(idxxxx2);

    for uu = 1:length(idxxxx2);
        confined_segments_length_local(uu) = length(Butterfly_trajectories_segments{idxxxx(idxxxx2(uu))});
    end
    
    if ~isempty(confined_segments_length_local);
    if ~isempty(find(confined_segments_length_local < Min_num_points_butt_confined_segment));
        discard_idx(i) = 1;
        discard_idxxxx = [discard_idxxxx idxxxx];
    else
        discard_idx(i) = 0;
    end
    end
end  
    clear confined_segments_length_local;
end

Butterfly_trajectories(find(discard_idx==1)) = [];
Butterfly_Reference_Track_segment(discard_idxxxx) = [];
motion_type_segment(discard_idxxxx) = [];
Butterfly_trajectories_segments(discard_idxxxx) = [];


%%3. Add The trajectories to msdanalyzer
ma_butterfly = ma_butterfly.addAll(Butterfly_trajectories);


%FIGURE 1: Plotting of the T-MSD for Butterfly Tracks only
% figure() 
% ma_butterfly.plotMSD;
% title('T-MSD Butterfly Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 2: Plotting of the mean T-MSD (black) and standard deviation (Grey) for Butterfly Tracks only
% figure()
% ma_butterfly.plotMeanMSD(gca, true);
% title('TE-MSD Butterfly Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 3: Visualization of Butterfly Tracks
% figure()
% ma_butterfly.plotTracks;
% ma_butterfly.labelPlotTracks;
% title('Butterfly Tracks');
% set(gca,'Ydir','reverse');

ma_butterfly_segments = ma_butterfly_segments.addAll(Butterfly_trajectories_segments');

%FIGURE 4: Visualization of the confined component of Butterfly Tracks
% figure()
% plotTracks(ma_butterfly_segments,gca,find(motion_type_segment==1));
% ma_butterfly_segments.labelPlotTracks;
% title('Butterfly Tracks Segmentated (Confined Segments)');
% set(gca,'Ydir','reverse');
% 

%FIGURE 5: Visualization of the directed component of Butterfly Tracks
% figure()
% plotTracks(ma_butterfly_segments,gca,find(motion_type_segment==0));
% ma_butterfly_segments.labelPlotTracks;
% title('Butterfly Tracks Segmentated (Directed Segments)');
% set(gca,'Ydir','reverse');

ma_butterfly_segments_confined = ma_butterfly_segments_confined.addAll(Butterfly_trajectories_segments(find(motion_type_segment==1)));
ma_butterfly_segments_directed = ma_butterfly_segments_directed.addAll(Butterfly_trajectories_segments(find(motion_type_segment==0)));

%Compute TMSD and TEMSD
ma_butterfly = ma_butterfly.computeMSD;
ma_butterfly = ma_butterfly.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly = ma_butterfly.TMSD(TMSD_fitting_points);

ma_butterfly_segments_confined = ma_butterfly_segments_confined.computeMSD;
ma_butterfly_segments_confined = ma_butterfly_segments_confined.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments_confined = ma_butterfly_segments_confined.TMSD(TMSD_fitting_points);

ma_butterfly_segments_directed = ma_butterfly_segments_directed.computeMSD;
ma_butterfly_segments_directed = ma_butterfly_segments_directed.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments_directed = ma_butterfly_segments_directed.TMSD(TMSD_fitting_points);


end

%% ------------------------------------------------------------------------------
%Add the trajectories to the msdAnalyzer-----------------------------------
%------------------------------------------------------------------------------   
 if isempty(trajectories)==1;
        error('There are no tracks that fulfill the requirements')
 else
     
ma = ma.addAll(trajectories);
Num_Tracks = size(trajectories,2);

for j=1:Num_Tracks;
Track_Length(j) = size(trajectories{j},1);
Residence_Times(j)=Track_Length(j)*Frame_interval;
end

end
    
%% -------------------------------------------------------------
%Calculate the distance travelled by the Tracks-------------------------
%(This value depends on the Track Length of each Track, but can be informative)--------
%------------------------------------------------------------------
for ggg = 1:size(trajectories,2);
points_coord = trajectories{ggg}(:,2:3);
[max_dist{ggg}, min_dist{ggg}, avg_dist{ggg}] = distance_scatter(points_coord);
end


%% Compute TE-MSD
ma = ma.computeMSD;


%% Preliminary Plotting

% FIGURE 5: Visualization of all the motion tracks
% figure()
% ma_AllTracks.plotTracks;
% ma_AllTracks.labelPlotTracks;
% set(gca,'Ydir','reverse');
% title('All Trajectories')

%FIGURE 6: Plot of the average T-MSD (black) and standard deviation (grey) for all motion tracks
% figure()
% ma.plotMeanMSD(gca, true)
% mmsd = ma.getMeanMSD;
% temps = mmsd(:,1);
% xs = mmsd(:,2);
% dx_plot = mmsd(:,3) ./ sqrt(mmsd(:,4));
% dxs = mmsd(:,3);


%FIGURE 7: Distribution of tracklength for all motion tracks
% figure()
% hist(Track_Length,100);
% title('Histogram of the Length of the Tracks');
% xlabel('Track length (frames)');
% ylabel('Frequency');
% xlim([0 100]);
% set(gca,'FontSize',20,'FontWeight','bold');

%% Compute T-MSD
%This other approach fits every single MSD curve and then the histogram
%of D is very informative
ma = ma.TMSD(TMSD_fitting_points);
good_enough_fit_Ds = find(ma.lfit.r2fit >= R2LIMIT);
Dmean = mean( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit_Ds) ) / 2 / ma.n_dim;
fprintf('**Estimation of the diffusion coefficient from linear fit of the MSD curves (Fitting every MSD curve)**:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, length(good_enough_fit_Ds));
Ds = ma.lfit.a(good_enough_fit_Ds)/ 2 / ma.n_dim;


%% FIGURE 9: Distribution of Diffusion Coefficients for all motion tracks
% figure()
% %Take out negative values from the Diffusion Coefficients list
% idx2=find(Ds > 0);
% histogram(log10(Ds(idx2)),50);
% %xlim([0 0.01]);
% title('Histogram of diffusion coefficients');
% xlabel('Log10(Diffusion Coefficient)')
% ylabel('Frequency');
% set(gca,'FontSize',20,'FontWeight','bold');


%% Motion type analysis from T-MSD curves (fitting each Track MSD)
ma = ma.LogTMSD(TLOGLOG_fitting_points);

r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
good_enough_fit_alpha = find(r2fits >= R2LIMIT);
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas_filtered = alphas(good_enough_fit_alpha);
%Remove NaN from alphas


%% Separate and plot Tracks based on their alpha value from LogLog fit of the T-MSD---------------
%-------------------------------------------------------------------------------------
temporal_tray = trajectories(good_enough_fit_alpha);
conf_filtered = find(alphas_filtered <= max_alpha_conf);
brownian_filtered = find(alphas_filtered>max_alpha_conf & alphas_filtered<min_alpha_directed);
directed_filtered = find(alphas_filtered >= min_alpha_directed);

%---------------------------------------------------------------------------------
%------------------------------------------------------------------------------
Conf_trajectories_filtered = temporal_tray(conf_filtered);
Directed_trajectories_filtered = temporal_tray(directed_filtered);
Brownian_trajectories_filtered = temporal_tray(brownian_filtered);

ma_confined = ma_confined.addAll(Conf_trajectories_filtered);
ma_directed = ma_directed.addAll(Directed_trajectories_filtered);
ma_brownian = ma_brownian.addAll(Brownian_trajectories_filtered);


%Calculate the number of Tracks per type of motion
Percentage_confined = size(ma_confined.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_brownian = size(ma_brownian.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_directed = size(ma_directed.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));
Percentage_butterfly = size(ma_butterfly.tracks,1)*100/(size(ma_brownian.tracks,1) + size(ma_directed.tracks,1) + size(ma_confined.tracks,1) + size(ma_butterfly.tracks,1));


%Compute MSD for all the motion types in separate variables
if isempty(ma_confined.tracks);   
else
ma_confined = ma_confined.computeMSD;
ma_confined = ma_confined.LogTMSD(TLOGLOG_fitting_points);
ma_confined = ma_confined.TMSD(TMSD_fitting_points);
end

if isempty(ma_directed.tracks);   
else
ma_directed = ma_directed.computeMSD;
ma_directed = ma_directed.LogTMSD(TLOGLOG_fitting_points);
ma_directed = ma_directed.TMSD(TMSD_fitting_points);
end

if isempty(ma_brownian.tracks);   
else
ma_brownian= ma_brownian.computeMSD;
ma_brownian = ma_brownian.LogTMSD(TLOGLOG_fitting_points);
ma_brownian = ma_brownian.TMSD(TMSD_fitting_points);
end


%% Separate the Confined Tracks based on their Mean Jump to eliminate "false confined" tracks.
for i=1:length(Conf_trajectories_filtered);  
    tracktemp = Conf_trajectories_filtered{i};
    jd_temp_confined_mean(i) = mean(sqrt((tracktemp(2:end,2) - tracktemp(1:end-1,2)).^2 + (tracktemp(2:end,3) - tracktemp(1:end-1,3)).^2));
end
idx6 = find(jd_temp_confined_mean >= Jump_confined_threshold);
idx7 = find(jd_temp_confined_mean < Jump_confined_threshold);

trajectories_confined_High_D = Conf_trajectories_filtered(idx6);
ma_confined_High_D = ma_confined_High_D.addAll(trajectories_confined_High_D);

trajectories_confined_Low_D = Conf_trajectories_filtered(idx7);
ma_confined_Low_D = ma_confined_Low_D.addAll(trajectories_confined_Low_D);


if isempty(ma_confined_High_D.tracks);  
else
ma_confined_High_D = ma_confined_High_D.computeMSD;
ma_confined_High_D = ma_confined_High_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_High_D = ma_confined_High_D.TMSD(TMSD_fitting_points);

%FIGURE 10: Visualize "false confined" motions
% figure()
% ma_confined_High_D.plotTracks;
% ma_confined_High_D.labelPlotTracks;
% title('Confined Trajectories with Average Jump higher than Jump_threshold','Interpreter', 'none');
% set(gca,'Ydir','reverse');


end

if isempty(ma_confined_Low_D.tracks);  
else
ma_confined_Low_D = ma_confined_Low_D.computeMSD;
ma_confined_Low_D = ma_confined_Low_D.LogTMSD(TLOGLOG_fitting_points);
ma_confined_Low_D = ma_confined_Low_D.TMSD(TMSD_fitting_points);

%FIGURE 11: Visualize "true confined" tracks
% figure()
% ma_confined_Low_D.plotTracks;
% ma_confined_Low_D.labelPlotTracks;
% title('Confined Trajectories with Average Jump lower than Jump_threshold','Interpreter', 'none');
% set(gca,'Ydir','reverse');
end

if isempty(ma_butterfly_segments.tracks);  
else
ma_butterfly_segments = ma_butterfly_segments.computeMSD;
ma_butterfly_segments = ma_butterfly_segments.LogTMSD(TLOGLOG_fitting_points);
ma_butterfly_segments = ma_butterfly_segments.TMSD(TMSD_fitting_points);
end




%% Plotting ---------------------------------------------------
%--------------------------------------------------------------------------------------
%--------------------------------------------
%-------------------------------------------------------------------------------------
%FIGURE 12: Plot T-MSD of unfiltered confined tracks
% figure()
% ma_confined.plotMSD;
% title('T-MSD Confined Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 12: Plot average T-MSD (black) and standard deviation (grey) of unfiltered confined tracks
% figure()
% ma_confined.plotMeanMSD(gca, true);
% title('TE-MSD Confined Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 13: Visualize unfiltered confined tracks
% figure()
% ma_confined.plotTracks;
% ma_confined.labelPlotTracks;
% title('Confined Tracks');
% set(gca,'Ydir','reverse');

%FIGURE 14: Plot T-MSD for directed motion tracks
% figure()
% ma_directed.plotMSD;
% title('T-MSD Directed Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 15: Plot average T-MSD (black) and standard deviation (grey) of
%directed motion tracks
% figure()
% ma_directed.plotMeanMSD(gca, true);
% title('TE-MSD Directed Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 16: Visualize directed motion tracks
% figure()
% ma_directed.plotTracks;
% ma_directed.labelPlotTracks;
% title('Directed Tracks');
% set(gca,'Ydir','reverse');

%FIGURE 17: Plot T-MSD for pure Brownian motion tracks
% figure()
% ma_brownian.plotMSD;
% title('T-MSD Brownian Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 18: Plot average T-MSD (black) and standard deviation (grey) for pure Brownian motion tracks
% figure()
% ma_brownian.plotMeanMSD(gca, true);
% title('TE-MSD Brownian Tracks');
% xlim([0 xmax]);
% ylim([0 ymax]);

%FIGURE 19: Visualize pure Brownian tracks
% figure()
% ma_brownian.plotTracks;
% ma_brownian.labelPlotTracks;
% title('Brownian Tracks');
% set(gca,'Ydir','reverse');


%%   save the results from the MSD Analysis
fprintf('\n');     
fprintf('Saving...\n'); 
clear ma
ma = ma_AllTracks;
mkdir(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'));
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_AllTracks','.mat')),'ma');
clear ma;
ma = ma_brownian;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Brownian','.mat')),'ma');
clear ma;
ma = ma_confined;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Confined','.mat')),'ma');
clear ma;
ma = ma_directed;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Directed','.mat')),'ma');
clear ma
ma = ma_butterfly;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Butterfly','.mat')),'ma');
clear ma
ma = ma_confined_High_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_confined_High_Jump_',num2str(Jump_confined_threshold),'.mat')),'ma');
clear ma
ma = ma_confined_Low_D;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_confined_Low_Jump_',num2str(Jump_confined_threshold),'.mat')),'ma');
clear ma
ma = ma_butterfly_segments_confined;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Butterfly_segments_confined','.mat')),'ma');
clear ma
ma = ma_butterfly_segments_directed;
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Butterfly_segments_directed','.mat')),'ma');
save(fullfile(strcat(directory_name,'\MSD_Results_',num2str(minimum_track_length),'min_frames'),strcat('msd_results_',num2str(Frame_interval),'ms','_Butterfly_segments_reference_ID','.mat')),'Butterfly_Reference_Track_segment');

 
      