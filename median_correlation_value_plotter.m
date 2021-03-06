% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020.

% This program computes the median correlation values and their uncertainties 
% for each month and plots them generating a figure similar to Figure 8c in the
% manuscript 

%The user will need to change INPUT1-INPUT25 below prior to executing
%this program

%The inputs are the locations of files containing autox, autoy, and autoz
%matrices computed for each month.
%These files can be generated by executing 
%magnetic_field_autocorrelation_function_calculator.m.

clearvars 

INPUT1 = 'C:/auto_res/2014/SEP/';

INPUT2 = 'C:/auto_res/2014/OCT/';

INPUT3 = 'C:/auto_res/2014/NOV/';

INPUT4 = 'C:/auto_res/2014/DEC/';

INPUT5 = 'C:/auto_res/2015/JAN/';

INPUT6 = 'C:/auto_res/2015/FEB/';

INPUT7 = 'C:/auto_res/2015/MAR/';

INPUT8 = 'C:/auto_res/2015/APR/';

INPUT9 = 'C:/auto_res/2015/MAY/';

INPUT10 = 'C:/auto_res/2015/JUN/';

INPUT11 = 'C:/auto_res/2015/JUL/';

INPUT12 = 'C:/auto_res/2015/AUG/';

INPUT13 = 'C:/auto_res/2015/SEP/';

INPUT14 = 'C:/auto_res/2015/OCT/';

INPUT15 = 'C:/auto_res/2015/NOV/';

INPUT16 = 'C:/auto_res/2015/DEC/';

INPUT17 = 'C:/auto_res/2016/JAN/';

INPUT18 = 'C:/auto_res/2016/FEB/';

INPUT19 = 'C:/auto_res/2016/MAR/';

INPUT20= 'C:/auto_res/2016/APR/';

INPUT21 = 'C:/auto_res/2016/MAY/';

INPUT22 = 'C:/auto_res/2016/JUN/';

INPUT23 = 'C:/auto_res/2016/JUL/';

INPUT24 = 'C:/auto_res/2016/AUG/';

INPUT25 = 'C:/auto_res/2016/SEP/';

st_time = 1
en_time = 25
for inst=st_time:en_time

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

if isfile(input_file_directory)==1
load(input_file_directory)

A = exist('autox')
if A==1
    
c = find(autox==0);
autox(c) = NaN;
c = find(autoy==0);
autoy(c) = NaN;
c = find(autoz==0);
autoz(c) = NaN;

autov = (autox+autoy+autoz)/3.0;
autov_q = quantile(autov,[0.25,0.5,0.75],1);
hold on

tt=find(autov_q(2,:) > exp(-1));
tau_c = max(incre_t_vec(1,tt));
tt=find(autov_q(1,:) > exp(-1));
tau_c_minus = max(incre_t_vec(1,tt));
tt=find(autov_q(3,:) > exp(-1));
tau_c_plus = max(incre_t_vec(1,tt));
if isempty(tau_c_minus)==1
    tau_c_minus=0;
end

YM(1,inst)= tau_c-tau_c_minus;
Y(1,inst)=tau_c;
YP(1,inst)=tau_c_plus-tau_c;

end
                                                                                          
end
end

X = [1:1:25];
errorbar(X,Y,YM,YP)
title('Median correlation time vs. time')
xlabel('correlation time (s)')
ylabel('time(mm/yy)')
xticks([1 3 5 7 9 11 13 15 17 19 21 23 25])
xticklabels({'9/14','11/14','1/15','3/15','5/15','7/15','9/15','11/15','1/16','3/16','5/16','7/16','9/16'})
