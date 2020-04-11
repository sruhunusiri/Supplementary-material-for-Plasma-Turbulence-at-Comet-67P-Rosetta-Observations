% Supplementary material to Suranga Ruhunusiri, "Plasma Turbulence at Comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to Journal of Geophysical Resarch Space Physics, 2020

% This program has been tested with MATLAB R2019a.

% This program reads RPCMAG data files, load data  
% into MATLAB Workspace, and save them as .mat files. 

% Before executing the code, the user will need to edit INPUT 1 and INPUT 2
% below.

% INPUT 1: directory: location of the Rosetta MAG data file. 
%          The latest ROSETTA RPCMAG calibrated data can be downloaded from 
%          ftp://psa.esac.esa.int/pub/mirror/INTERNATIONAL-ROSETTA-MISSION/RPCMAG/
directory  = 'C:/Rosetta/MAG/';

% 
% INPUT 2: output_file_directory: location of the MAG data to be saved in .mat format. 
%          This .mat file will contain row matrices for date in seconds,
%          spacecarft position (x,y,z) in CSEQ coordinate system, 
%          magnetic field componenents in CSEQ coordinate system, and
%          eight quality flags
 
output_file_directory = 'C:/Rosetta/MAT_MAG/';

filesAndFolders = dir(directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));
file_i=1;
data_file = filesInDir(file_i).name;                              
FID = fopen(strcat(directory,data_file), 'r');
data = textscan(FID,'%s %f %f %f %f %f %f %f %s');
fclose(FID);
% 
time = cell2mat(data{1,1});
posx =data{1,3}; 
posy =data{1,4};
posz =data{1,5};
Bx = data{1,6};
By = data{1,7};
Bz = data{1,8};
Qflags = cell2mat(data{1,9});

 
 Year = zeros(size(Bx,1),1);
 Month = zeros(size(Bx,1),1);
 Day = zeros(size(Bx,1),1);
 hour = zeros(size(Bx,1),1);
 mint = zeros(size(Bx,1),1);
 sec = zeros(size(Bx,1),1);
 Date_num = zeros(size(Bx,1),1);
 QF1 = zeros(size(Bx,1),1);
 QF2 = zeros(size(Bx,1),1);
 QF3 = zeros(size(Bx,1),1);
 QF4 = zeros(size(Bx,1),1);
 QF5 = zeros(size(Bx,1),1);
 QF6 = zeros(size(Bx,1),1);
 QF7 = zeros(size(Bx,1),1);
 QF8 = zeros(size(Bx,1),1);
% 
for i=1:size(Bx,1)
%   
     Year(i,1) = str2double(time(i,1:4));
     Month(i,1) = str2double(time(i,6:7));
     Day(i,1) = str2double(time(i,9:10));
     hour(i,1) = str2double(time(i,12:13));
     mint(i,1) = str2double(time(i,15:16));
     sec(i,1) = str2double(time(i,18:26));
     Date_num(i,1) = datenum(Year(i,1),Month(i,1),Day(i,1),hour(i,1),mint(i,1),sec(i,1));
     QF1(i,1) = str2double(Qflags(i,1));
     QF2(i,1) = str2double(Qflags(i,2));
     QF3(i,1) = str2double(Qflags(i,3));
     QF4(i,1) = str2double(Qflags(i,4));
     QF5(i,1) = str2double(Qflags(i,5));
     QF6(i,1) = str2double(Qflags(i,6));
     QF7(i,1) = str2double(Qflags(i,7));
     QF8(i,1) = str2double(Qflags(i,8));
 end
% 
 save(strcat(output_file_directory,'RPCMAG_data_',num2str(Year(1,1)),'_',num2str(Month(1,1)),'_',num2str(Day(1,1))),'Date_num','posx','posy','posz','Bx','By','Bz','QF1','QF2','QF3','QF4','QF5','QF6','QF7','QF8','-v7.3');
% 
 clearvars -except Date_num posx posy posz Bx By Bz QF1 QF2 QF3 QF4 QF5 QF6 QF7 QF8;