clear;close;clc;close all;
cd('D:\Hassett\Matlab\');%Directory where your data are
%addpath('C:\Users\ehassett1\Desktop\Matlab\Code\'); %This should point to wherever you keep all the codes that support this. Assuming they're not in the working directory
% Read data
dat_num='08182022'; %Date vector of when the measurements were taken

%% Read the file.
metadata=ImportHWF_BeaverDamChambers(['.\Data\Trip' dat_num '\flux_metadata_' dat_num '.csv']); %It must be in a folder called "Data" which must then have a folder called "Trip06292022" or something like that. 

%% Makes the calculated flux folder if it does not already exist
f=dir(['.\Data\Trip08182022' dat_num '\Flux' dat_num]);
if isempty(f)
   mkdir(['.\Data\Trip' dat_num '\Flux' dat_num]);
end
%% Set up data structure for output file
fname=['.\Data\Trip' dat_num '\Flux' dat_num '\FluxOutput' dat_num '.csv'];%fid= fopen(fname,'w'); % w= write access
filepath = cd;
file= fullfile(filepath, fname);
fid = fopen(file, 'wt');
fprintf(fid,'Written by Erin Hassett\n'); %n = new line
fprintf(fid,[datestr(now()) '\n']);
fprintf(fid,'Date,Site,CO2_Flux,CO2_NRMSE,CH4_Flux,CH4_NRMSE,Bubble_Flux\n');
fprintf(fid,'Date,Site,umol m2 s-1,%%,umol m2 s-1,%%,umol m2 s-1\n');
fclose(fid);

for i=1:length(metadata.Date)
    %% Chamber analysis for i_th chamber
    [Flux_HM_CO2, NRMSE_CO2,Flux_HM_CH4, NRMSE_CH4,Bubble_CH4]=ProcessPicarroChamber(... 
        metadata.Date(i),metadata.start_time(i),metadata.end_time(i),...
        ['.\Data\Trip' dat_num '\' metadata.data_sheet{i}],metadata.Start_offset(i),metadata.End_offset(i),... 
        metadata.ChamberVolume(i),metadata.ChamberArea(i),1,metadata.Name{i},metadata.Start_Temp(i),metadata.End_Temp(i),...
        metadata.Shift_back(i));
	%% Append data to output file
    fid = fopen(fname,'a+'); % a+ = Open or create new file for reading and writing. Append data to the end of the file.
    fprintf(fid,'%s,%s,%f,%f,%f,%f,%f\n',dat_num,metadata.Name{i},Flux_HM_CO2,NRMSE_CO2,Flux_HM_CH4,NRMSE_CH4,Bubble_CH4);
    fclose(fid);
end
