% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020.

% This program computes the spatial occurrence rates of turbulent processes 
% TP1-TP4 and plots them generating a figures similar to Figure 7 in
% the manuscript.

% Before executing this program, the user may need to change INPUTS 1-25
% below, which are the file locations that contain Data_holder.mat. 
% These files can be generated by executing the RPCMAG_data_analysis.m. 

% The user should execute this program in current form, to plot spatial 
% maps for the intermediately and strongly active phases of the comet.
% The user should change st_time to 1 and en_time to 25 in lines 71-72, 
% and uncomment lines 74 and 213 to plot spatial maps for the weakly active 
% phase of the comet.

clearvars

INPUT1 = 'C:/Rosetta/anlysis_folder/2014/SEP/';

INPUT2 = 'C:/Rosetta/anlysis_folder/2014/OCT/';

INPUT3 = 'C:/Rosetta/anlysis_folder/2014/NOV/';

INPUT4 = 'C:/Rosetta/anlysis_folder/2014/DEC/';

INPUT5 = 'C:/Rosetta/anlysis_folder/2015/JAN/';

INPUT6 = 'C:/Rosetta/anlysis_folder/2015/FEB/';

INPUT7 = 'C:/Rosetta/anlysis_folder/2015/MAR/';

INPUT8 = 'C:/Rosetta/anlysis_folder/2015/APR/';

INPUT9 = 'C:/Rosetta/anlysis_folder/2015/MAY/';

INPUT10 = 'C:/Rosetta/anlysis_folder/2015/JUN/';

INPUT11 = 'C:/Rosetta/anlysis_folder/2015/JUL/';

INPUT12 = 'C:/Rosetta/anlysis_folder/2015/AUG/';

INPUT13 = 'C:/Rosetta/anlysis_folder/2015/SEP/';

INPUT14 = 'C:/Rosetta/anlysis_folder/2015/OCT/';

INPUT15 = 'C:/Rosetta/anlysis_folder/2015/NOV/';

INPUT16 = 'C:/Rosetta/anlysis_folder/2015/DEC/';

INPUT17 = 'C:/Rosetta/anlysis_folder/2016/JAN/';

INPUT18 = 'C:/Rosetta/anlysis_folder/2016/FEB/';

INPUT19 = 'C:/Rosetta/anlysis_folder/2016/MAR/';

INPUT20= 'C:/Rosetta/anlysis_folder/2016/APR/';

INPUT21 = 'C:/Rosetta/anlysis_folder/2016/MAY/';

INPUT22 = 'C:/Rosetta/anlysis_folder/2016/JUN/';

INPUT23 = 'C:/Rosetta/anlysis_folder/2016/JUL/';

INPUT24 = 'C:/Rosetta/anlysis_folder/2016/AUG/';

INPUT25 = 'C:/Rosetta/anlysis_folder/2016/SEP/';

st_time = 1;
en_time = 25;
for inst=st_time:en_time
%if inst<=6 | inst>=18

if inst==1
input_file_directory = INPUT1;
end

if inst==2
input_file_directory = INPUT2;
end

if inst==3
input_file_directory = INPUT3;
end

if inst==4
input_file_directory = INPUT4;
end

if inst==5
input_file_directory = INPUT5;
end

if inst==6
input_file_directory = INPUT6;
end

if inst==7
input_file_directory = INPUT7;
end

if inst==8
input_file_directory = INPUT8;
end

if inst==9
input_file_directory = INPUT9;
end

if inst==10
input_file_directory = INPUT10;
end

if inst==11
input_file_directory = INPUT11;
end

if inst==12
input_file_directory = INPUT12;
end

if inst==13
input_file_directory = INPUT13;
end

if inst==14
input_file_directory = INPUT14;
end

if inst==15
input_file_directory = INPUT15;
end

if inst==16
input_file_directory = INPUT16;
end

if inst==17
input_file_directory = INPUT17;
end

if inst==18
input_file_directory = INPUT18;
end

if inst==19
input_file_directory = INPUT19;
end

if inst==20
input_file_directory = INPUT20;
end

if inst==21
input_file_directory = INPUT21;
end

if inst==22
input_file_directory = INPUT22;
end

if inst==23
input_file_directory = INPUT23;
end

if inst==24
input_file_directory = INPUT24;
end

if inst==25
input_file_directory = INPUT25;
end
                                                                                          
jjj=1;
filesAndFolders = dir(input_file_directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));
num_files = size(filesInDir);
num_files = num_files(1,1);

for file_i=jjj:num_files
    load(strcat(filesInDir(file_i).folder,'/',filesInDir(file_i).name));
    if file_i==jjj & inst==st_time
    
    spec_low = Data_holder(:,1);
    low_err = Data_holder(:,2);
    spec_high = Data_holder(:,3);
    high_err = Data_holder(:,4);
    spec_valid = Data_holder(:,6);
    
     XPOS = Data_holder(:,7);
     YPOS = Data_holder(:,8);
     ZPOS = Data_holder(:,9);
    end
    if file_i>jjj | inst>st_time
       
    spec_low = [spec_low;Data_holder(:,1)];
    low_err = [low_err;Data_holder(:,2)];
    spec_high = [spec_high;Data_holder(:,3)];
    high_err = [high_err;Data_holder(:,4)];
    spec_valid = [spec_valid ;Data_holder(:,6)];
        
        XPOS = [XPOS;Data_holder(:,7)];
        YPOS = [YPOS;Data_holder(:,8)];
        ZPOS = [ZPOS;Data_holder(:,9)];
    end
    
end

end

%end

fin = find(spec_valid==1);

spec_low = spec_low(fin,1);
low_err = low_err(fin,1);
spec_high = spec_high(fin,1);
high_err = high_err(fin,1);
XPOS = XPOS(fin,1);
YPOS = YPOS(fin,1);
ZPOS = ZPOS(fin,1);

for p_num=1:4
find_process = zeros(length(fin),1);
uncertain_minus = zeros(length(fin),1);
uncertain_plus = zeros(length(fin),1);
for ch=1:length(fin)
    
    spec_index_low = spec_low(ch,1);
    spec_index_high = spec_high(ch,1);
    spec_index_low_err = low_err(ch,1);
    spec_index_high_err = high_err(ch,1);
    
    %TP1
    if p_num==1
        if spec_index_low>=-0.5 & spec_index_high<-1.7
            find_process(ch,1) = 1; 
            
            if spec_index_low-spec_index_low_err<-0.5 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,1) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err>=-0.5 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err>=-0.5 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err>=-0.5 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err>=-0.5 & spec_index_high-spec_index_high_err<-1.7)
            
           uncertain_plus(ch,1) = 1; 
            
            end
            
        end
    end
    
    %TP2
    if p_num==2
        
        if spec_index_low<-0.5 & spec_index_low>-1.4 & spec_index_high<-1.7
            find_process(ch,1) = 1; 
            
            if spec_index_low+spec_index_low_err>=-0.5 | spec_index_low-spec_index_low_err<=-1.4 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,1) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<-1.7) | ...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<-1.7)
            uncertain_plus(ch,1) = 1; 
            
            end
            
        end
        
        
        
    end
    
    %TP3
    if p_num==3
        
         if spec_index_low<-1.7 & spec_index_high<-1.7
            find_process(ch,1) = 1; 
            
            if spec_index_low+spec_index_low_err>=-1.7 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,1) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<-1.7 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-1.7 & spec_index_high-spec_index_high_err<-1.7)
               
            uncertain_plus(ch,1) = 1; 
            
            end
            
        end
        
        
    end
        
    %TP4
    if p_num==4
        
        if spec_index_low<=-1.4 & spec_index_low>=-1.7 & spec_index_high<-1.7
            find_process(ch,1) = 1; 
            
            if spec_index_low+spec_index_low_err>-1.4 | spec_index_low-spec_index_low_err<-1.7 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,1) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<-1.7)
            uncertain_plus(ch,1) = 1; 
            
            end
            
        end
        
        
        
    end
        
        
    
end


sel_process = (find(find_process==1));

sel1 = find(XPOS>min(XPOS));
XPOS_sel = XPOS(sel1,1);
RPOS_sel = (YPOS(sel1,1).^2+ZPOS(sel1,1).^2).^(1/2);
theta_sel = atan(RPOS_sel./XPOS_sel);
RPOS_sel = (RPOS_sel.^2+XPOS_sel.^2).^(1/2);
dr = find(theta_sel<0);
theta_sel(dr,1) = theta_sel(dr,1)+pi;

x = [0.5,100,250,500,750,1000,1250,1500,1750];
theta = pi*[0,22.5,45,67.5,90,112.5,135,157.5,180]/180;

tracker = zeros(8,8);
for out=1:length(RPOS_sel)
for i=1:8
    
    for j=1:8

if RPOS_sel(out,1)>=x(1,i) & RPOS_sel(out,1)<x(1,i+1) & theta_sel(out,1)>=theta(1,j) & theta_sel(out,1)<theta(1,j+1)
tracker(i,j) = tracker(i,j)+1;
end
    end
  
end
end

sel1 = sel_process;

XPOS_sel = XPOS(sel1,1);
RPOS_sel = (YPOS(sel1,1).^2+ZPOS(sel1,1).^2).^(1/2);
theta_sel = atan(RPOS_sel./XPOS_sel);
RPOS_sel = (RPOS_sel.^2+XPOS_sel.^2).^(1/2);
dr = find(theta_sel<0);
theta_sel(dr,1) = theta_sel(dr,1)+pi;

x = [0.5,100,250,500,750,1000,1250,1500,1750];
theta = pi*[0,22.5,45,67.5,90,112.5,135,157.5,180]/180;

tracker2 = zeros(8,8);
for out=1:length(RPOS_sel)
for i=1:8
    
    for j=1:8

if RPOS_sel(out,1)>=x(1,i) & RPOS_sel(out,1)<x(1,i+1) & theta_sel(out,1)>=theta(1,j) & theta_sel(out,1)<theta(1,j+1)
tracker2(i,j) = tracker2(i,j)+1;
end
    end
    
end
end
temp_tracker = 100*tracker2./tracker;
temp_trackerp = temp_tracker;
temp_tracker = 63*temp_tracker./100;

if p_num==1
    subplot(4,1,1)

end
if p_num==2
    subplot(4,1,2)

end
if p_num==3
    subplot(4,1,3)

end
if p_num==4
    subplot(4,1,4)

end

colormap jet
cmap = colormap(jet);
x = [0.5,100,250,500,750,1000,1250,1500,1750];
theta = pi*[0,22.5,45,67.5,90,112.5,135,157.5,180]/180;

for i=1:8
    
    for j=1:8
      
           binsx = [x(1,i)*cos(theta(1,j)),x(1,i)*cos((theta(1,j+1)+3.0*theta(1,j))/4.0),x(1,i)*cos((theta(1,j+1)+theta(1,j))/2.0),x(1,i)*cos((3.0*theta(1,j+1)+theta(1,j))/4.0),x(1,i)*cos(theta(1,j+1)),x(1,i+1)*cos(theta(1,j+1)),x(1,i+1)*cos(0.75*theta(1,j+1)+.25*theta(1,j)),x(1,i+1)*cos((theta(1,j+1)+theta(1,j))/2.0),x(1,i+1)*cos((theta(1,j+1)+3.0*theta(1,j))/4.0),x(1,i+1)*cos(theta(1,j)),x(1,i)*cos(theta(1,j))];
           binsy = [x(1,i)*sin(theta(1,j)),x(1,i)*sin((theta(1,j+1)+3.0*theta(1,j))/4.0),x(1,i)*sin((theta(1,j+1)+theta(1,j))/2.0),x(1,i)*sin((3.0*theta(1,j+1)+theta(1,j))/4.0),x(1,i)*sin(theta(1,j+1)),x(1,i+1)*sin(theta(1,j+1)),x(1,i+1)*sin(0.75*theta(1,j+1)+.25*theta(1,j)),x(1,i+1)*sin((theta(1,j+1)+theta(1,j))/2.0),x(1,i+1)*sin((theta(1,j+1)+3.0*theta(1,j))/4.0),x(1,i+1)*sin(theta(1,j)),x(1,i)*sin(theta(1,j))];
      
           
           pgon = polyshape(binsx,binsy);
         
         val = ceil(temp_tracker(i,j));
         
         if val>0
         fill(binsx,binsy,cmap(val,:))
         end
         if val==0
            fill(binsx,binsy,[1,1,1]) 
         end
         
         hold on

    end
    
end

if p_num==1
title('TP1 occurrence rate (%)')
end
if p_num==2
title('TP2 occurrence rate (%)')
end
if p_num==3
title('TP3 occurrence rate (%)')
end
if p_num==4
title('TP4 occurrence rate (%)')
end
xlabel('X-CSEQ (km)')
ylabel('R-CSEQ (km)')

end

clearvars

