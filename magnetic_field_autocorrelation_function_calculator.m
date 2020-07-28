% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020.

% This program computes the autocorrelation function for magnetic field, 
% which can be used to generate plots similar to Figures 8a and 8b in the
% manuscript

%Prior to executing this program, the user may need to change
%input_file_directory and output_file_directory below.

input_file_directory = 'C:/MAT_MAG/2014/Sep/';
output_file_directory = 'C:/auto_cor_res/2014/Sep/';

filesAndFolders = dir(input_file_directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));
num_files = size(filesInDir);
num_files = num_files(1,1);
auto_index = 0;

for file_i=1:num_files
   
load(strcat(filesInDir(file_i).folder,'/',filesInDir(file_i).name));
[tok,rem] = strtok(filesInDir(file_i).name,'_');


dur = datenum(2014,11,28,4,0,0)-datenum(2014,11,28,0,0,0);
dur2 = datenum(2014,11,28,0,30,0)-datenum(2014,11,28,0,0,0);
row_select = 1;
     
start_date = Date_num(1,1);
size_date_num = size(Date_num);
time_segments = (Date_num(size_date_num(1,1),1)-start_date-dur)/dur2;
time_segments = time_segments-mod(time_segments,1);

for time_incr=1:time_segments

start_t = start_date+dur2*(time_incr-1);
end_t = start_t+dur;
indices = find(Date_num >= start_t & Date_num <= end_t);
if length(indices) >=287000 && length(indices) <=288000+10
    timet_now = Date_num(indices,1); 

    timet_now_new = (timet_now-timet_now(1,1))/(datenum(2014,11,28,0,0,0.05)-datenum(2014,11,28,0,0,0));
    timet_now_new  = floor(timet_now_new)+1;

    time_vector_new = [0:1/20:(max(timet_now_new)-1)/20];
    time_vector_new = time_vector_new';
    Bx_now = Bx(indices,1);
    By_now = By(indices,1);
    Bz_now = Bz(indices,1);
      
    Bx_new = NaN(max(timet_now_new),1);
    Bx_new(timet_now_new,1) = Bx_now;
    Bx_new= resample(Bx_new, time_vector_new);
    By_new = NaN(max(timet_now_new),1);
    By_new(timet_now_new,1) = By_now;
    By_new= resample(By_new, time_vector_new);
    Bz_new = NaN(max(timet_now_new),1);
    Bz_new(timet_now_new,1) = Bz_now;
    Bz_new= resample(Bz_new, time_vector_new);
  
 if length(Bx_new)>=288000
       Bx_now = Bx_new(1:288000,:);
       By_now = By_new(1:288000,:);
       Bz_now = Bz_new(1:288000,:);
 
        Bxx = Bx_now-mean(Bx_now);
        Byy = By_now-mean(By_now);
        Bzz= Bz_now-mean(Bz_now);

        auto_index = auto_index+1;
        
 for kkk=1:500
     incre_t = kkk*9;
     incre_t_vec(1,kkk) = incre_t;
 dur3 = datenum(2014,11,28,4,0,0)-datenum(2014,11,28,0,0,0);
 dur4 = datenum(2014,11,28,0,0,incre_t)-datenum(2014,11,28,0,0,0);

    start_date = Date_num(1,1);
    size_date_num = size(Date_num);
 
 time_incr2=2;
 start_t = start_date+dur4*(time_incr2-1)+dur2*(time_incr-1);
 end_t = start_t+dur3;
 indices = find(Date_num >= start_t & Date_num <= end_t);
    if length(indices) >=287000 && length(indices) <=288000+10

    timet_now = Date_num(indices,1); 

    timet_now_new = (timet_now-timet_now(1,1))/(datenum(2014,11,28,0,0,0.05)-datenum(2014,11,28,0,0,0));
    timet_now_new  = floor(timet_now_new)+1;

    time_vector_new = [0:1/20:(max(timet_now_new)-1)/20];
    time_vector_new = time_vector_new';
    Bx_now = Bx(indices,1);
    By_now = By(indices,1);
    Bz_now = Bz(indices,1);
     
   Bx_new = NaN(max(timet_now_new),1);
   Bx_new(timet_now_new,1) = Bx_now;
   Bx_new= resample(Bx_new, time_vector_new);
   By_new = NaN(max(timet_now_new),1);
   By_new(timet_now_new,1) = By_now;
   By_new= resample(By_new, time_vector_new);
   Bz_new = NaN(max(timet_now_new),1);
   Bz_new(timet_now_new,1) = Bz_now;
   Bz_new= resample(Bz_new, time_vector_new);

if length(Bx_new)>=288000
       Bx_now = Bx_new(1:288000,:);
       By_now = By_new(1:288000,:);
       Bz_now = Bz_new(1:288000,:);
 
        Bxx2 = Bx_now-mean(Bx_now);

        Byy2 = By_now-mean(By_now);

        Bzz2= Bz_now-mean(Bz_now);

 autox( auto_index,kkk) = mean(Bxx.*Bxx2)/mean(Bxx.*Bxx);
 autoy( auto_index,kkk) = mean(Byy.*Byy2)/mean(Byy.*Byy);
 autoz( auto_index,kkk) = mean(Bzz.*Bzz2)/mean(Bzz.*Bzz);
 
end
end
end
end
end
end
end

save(strcat(output_file_directory,'auto_res'),'autox','autoy','autoz','-v7.3');