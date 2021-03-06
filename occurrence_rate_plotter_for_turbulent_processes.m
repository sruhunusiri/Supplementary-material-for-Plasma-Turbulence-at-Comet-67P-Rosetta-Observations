% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020.

% This program computes monthly occurrence rates and their uncertainties
% for processes TP1-TP8 and plots them generating figures similar to 
% Figures 5 and 6 in the manuscript.

% Before executing this program, the user may need to change INPUTS 1-25
% below, which are the file locations that contain Data_holder.mat
% Data_holder.mat can be generated by executing the RPCMAG_data_analysis.m. 

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

process_p = zeros(8,25);
process_pos_error = zeros(8,25);
process_neg_error = zeros(8,25);
for ijk=1:25

clearvars -except ijk process_p process_pos_error process_neg_error INPUT1...
    INPUT2 INPUT3 INPUT4 INPUT5 INPUT6 INPUT7 INPUT8 INPUT9 INPUT10 INPUT11...
    INPUT12 INPUT13 INPUT14 INPUT15 INPUT16 INPUT17 INPUT18 INPUT19 INPUT20...
    INPUT21 INPUT22 INPUT23 INPUT24 INPUT25 

inst =ijk;
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
    if file_i==jjj
     
    spec_low = Data_holder(:,1);
    low_err = Data_holder(:,2);
    spec_high = Data_holder(:,3);
    high_err = Data_holder(:,4);
    spec_valid = Data_holder(:,6);
    
    end
    
    if file_i>jjj
        
    spec_low = [spec_low;Data_holder(:,1)];
    low_err = [low_err;Data_holder(:,2)];
    spec_high = [spec_high;Data_holder(:,3)];
    high_err = [high_err;Data_holder(:,4)];
    spec_valid = [spec_valid ;Data_holder(:,6)];
     
    end
    
end

fin = find(spec_valid==1);

if length(fin)>0

spec_low = spec_low(fin,1);
low_err = low_err(fin,1);
spec_high = spec_high(fin,1);
high_err = high_err(fin,1);

end

find_process = zeros(length(fin),8);
uncertain_minus = zeros(length(fin),8);
uncertain_plus = zeros(length(fin),8);

for ch=1:length(fin)
    
    spec_index_low = spec_low(ch,1);
    spec_index_high = spec_high(ch,1);
    spec_index_low_err = low_err(ch,1);
    spec_index_high_err = high_err(ch,1);
    
        %dominant energy injection-partial cascade or 
        %dominant energy injection-dispersion 
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
        
       %uniform energy injection-partial cascade or 
       %uniform energy injection-dispersion
       if spec_index_low<-0.5 & spec_index_low>-1.4 & spec_index_high<-1.7
            find_process(ch,2) = 1; 
            
            if spec_index_low+spec_index_low_err>=-0.5 | spec_index_low-spec_index_low_err<=-1.4 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,2) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<-1.7) | ...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<-1.7)
            uncertain_plus(ch,2) = 1; 
            
            end
            
       end
       
        %partial cascade-partial cascade, partial cascade-dispersion,
        %dissipation-dissipation, or dissipation-dispersion  
      if spec_index_low<-1.7 & spec_index_high<-1.7
            find_process(ch,3) = 1; 
            
            if spec_index_low+spec_index_low_err>=-1.7 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,3) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<-1.7 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<-1.7 & spec_index_high-spec_index_high_err<-1.7)
               
            uncertain_plus(ch,3) = 1; 
            
            end
            
        end
      
        %full cascade-dissipation or full cascade-dispersion
        if spec_index_low<=-1.4 & spec_index_low>=-1.7 & spec_index_high<-1.7
            find_process(ch,4) = 1; 
            
            if spec_index_low+spec_index_low_err>-1.4 | spec_index_low-spec_index_low_err<-1.7 | spec_index_high+spec_index_high_err>=-1.7
            uncertain_minus(ch,4) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<-1.7)
            uncertain_plus(ch,4) = 1; 
            
            end
            
        end
          
       %dominant injection-full cascade 
       if spec_index_low>=-0.5 & spec_index_high<=-1.4 & spec_index_high>=-1.7
            find_process(ch,5) = 1; 
            
            if spec_index_low-spec_index_low_err<-0.5 | spec_index_high+spec_index_high_err>-1.4 | spec_index_high-spec_index_high_err<-1.7
            uncertain_minus(ch,5) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err>=-0.5 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err>=-0.5 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low+spec_index_low_err>=-0.5 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err>=-0.5 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7)
            uncertain_plus(ch,5) = 1; 
            
            end
            
       end
    
       %uniform injection-full cascade
       if spec_index_low<-0.5 & spec_index_low>-1.4 & spec_index_high<=-1.4 & spec_index_high>=-1.7
            find_process(ch,6) = 1; 
            
            if spec_index_low+spec_index_low_err>=-0.5 |spec_index_low-spec_index_low_err<-1.4 | spec_index_high-spec_index_high_err<-1.7 | spec_index_high+spec_index_high_err> -1.4
            uncertain_minus(ch,6) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low+spec_index_low_err<-0.5 & spec_index_low+spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err<-0.5 & spec_index_low-spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7)
            uncertain_plus(ch,6) = 1; 
            
            end
            
       end
       
       %uniform or dominant injection-uniform or dominant injection
       if spec_index_low>-1.4 & spec_index_high>-1.4
            find_process(ch,7) = 1; 
            
            if spec_index_low-spec_index_low_err<=-1.4| spec_index_high-spec_index_high_err<=-1.4
            uncertain_minus(ch,7) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err>-1.4) |...
               (spec_index_low-spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err>-1.4) |...
               (spec_index_low+spec_index_low_err>-1.4 & spec_index_high+spec_index_high_err>-1.4) |...
               (spec_index_low-spec_index_low_err>-1.4 & spec_index_high-spec_index_high_err>-1.4)
            uncertain_plus(ch,7) = 1; 
            
            end
            
        end
              
      %full cascade-full cascade 
      if spec_index_low<=-1.4 & spec_index_low>=-1.7 & spec_index_high<=-1.4 & spec_index_high>=-1.7
            find_process(ch,8) = 1; 
            
            if spec_index_low+spec_index_low_err>-1.4 | spec_index_low-spec_index_low_err<-1.7 | spec_index_high+spec_index_high_err>-1.4 | spec_index_high-spec_index_high_err<1.7
            uncertain_minus(ch,8) = 1; 
            
            end
            
        else
            
            if (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low+spec_index_low_err<=-1.4 & spec_index_low+spec_index_low_err>=-1.7 & spec_index_high+spec_index_high_err<=-1.4 & spec_index_high+spec_index_high_err>=-1.7) |...
               (spec_index_low-spec_index_low_err<=-1.4 & spec_index_low-spec_index_low_err>=-1.7 & spec_index_high-spec_index_high_err<=-1.4 & spec_index_high-spec_index_high_err>=-1.7)
            uncertain_plus(ch,8) = 1; 
            
            end
            
        end
      
end


for pind=1:8
process_p(pind,inst) = length(find(find_process(:,pind)==1))*100/length(find_process(:,pind));
process_neg_error(pind,inst)= length(find(uncertain_minus(:,pind)==1))*100/length(find_process(:,pind));
process_pos_error(pind,inst)=length(find(uncertain_plus(:,pind)==1))*100/length(find_process(:,pind));
end
end


x = [1:1:25];
figure
for p_num = 1:4
subplot(4,1,p_num)
errorbar(x,process_p(p_num,:),process_neg_error(p_num,:),process_pos_error(p_num,:))
if p_num ==1
title('TP1-process occurrence rate vs. time')
end
if p_num ==2
title('TP2-process occurrence rate vs. time')
end
if p_num ==3
title('TP3-process occurrence rate vs. time')
end
if p_num ==4
title('TP4-process occurrence rate vs. time')
end
    
xlabel('time (mm/yy)')
ylabel('occurrence rate (%)')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})

end

figure
for p_num = 5:8
subplot(4,1,p_num-4)
errorbar(x,process_p(p_num,:),process_neg_error(p_num,:),process_pos_error(p_num,:))

if p_num ==5
title('TP5-process occurrence rate vs. time')
end
if p_num ==6
title('TP6-process occurrence rate vs. time')
end
if p_num ==7
title('TP7-process occurrence rate vs. time')
end
if p_num ==8
title('TP8-process occurrence rate vs. time')
end
    
xlabel('time (mm/yy)')
ylabel('occurrence rate (%)')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})

end

clearvars
