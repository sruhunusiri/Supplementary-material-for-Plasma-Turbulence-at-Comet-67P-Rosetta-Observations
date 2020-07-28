% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020.

% This program computes monthly medians and uncertainties of spectral break
% frequency, low-frequency spectral index, and high-frequency spectral
% index and plots them generating a figure similar to Figure 4 in the manuscript.

% Before executing this program, the user may need to change INPUTS 1-25
% below, which are the file locations that contain Data_holder.mat.
% Data_holder.mat can be generated by executing RPCMAG_data_analysis.m. 

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


for inst=1:25

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
                spec_high = Data_holder(:,3);
                break_f = Data_holder(:,5);
                spec_valid = Data_holder(:,6);
    
            end
    
            if file_i>jjj 
       
                spec_low = [spec_low;Data_holder(:,1)];
                spec_high = [spec_high;Data_holder(:,3)];
                break_f = [break_f;Data_holder(:,5)];
                spec_valid = [spec_valid;Data_holder(:,6)];
        
            end
    
    end
    
    fin = find(spec_valid==1);
    spec_low = spec_low(fin,1);
    spec_high = spec_high(fin,1);
    spec_break = break_f(fin,1);
    spec_break = 10.^spec_break;
    
    param = quantile(spec_break,[0.25,0.5,0.75]);
    break_median(1,inst) = param(1,2);
    break_pos_error(1,inst) = param(1,3) - param(1,2);
    break_neg_error(1,inst) = param(1,2) - param(1,1);
    
    param = quantile(spec_low,[0.25,0.5,0.75]);
    low_spec_median(1,inst) = param(1,2);
    low_spec_pos_error(1,inst) = param(1,3) - param(1,2);
    low_spec_neg_error(1,inst) = param(1,2) - param(1,1);
    
    param = quantile(spec_high,[0.25,0.5,0.75]);
    high_spec_median(1,inst) = param(1,2);
    high_spec_pos_error(1,inst) = param(1,3) - param(1,2);
    high_spec_neg_error(1,inst) = param(1,2) - param(1,1);

end
figure
x = [1:1:25];

subplot(3,1,1)
errorbar(x,break_median,break_neg_error,break_pos_error)
title('median spectral break frequency vs. time')
xlabel('time (mm/yy)')
ylabel('spectral-break frequency (Hz)')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})

subplot(3,1,2)
errorbar(x,low_spec_median,low_spec_neg_error,low_spec_pos_error)
title('low-frequency spectral index vs. time')
xlabel('time (mm/yy)')
ylabel('spectral-index value')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})
line([1,25],[-5/3,-5/3])

subplot(3,1,3)
errorbar(x,high_spec_median,high_spec_neg_error,high_spec_pos_error)
title('high-frequency spectral index vs. time')
xlabel('time (mm/yy)')
ylabel('spectral-index value')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})
line([1,25],[-5/3,-5/3])

clearvars