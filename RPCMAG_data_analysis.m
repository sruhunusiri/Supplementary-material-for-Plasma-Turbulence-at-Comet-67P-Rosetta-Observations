% Supplementary material to Suranga Ruhunusiri, G. G. Howes, & J. S. Halekas' 
% "Plasma Turbulence at comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submitted to JGR Space Physics on 04/11/2020

% This program computes power spectra for magnetic field fluctuations, determine 
% the spectral break frequency, low- and high-frequency spectral
% indices, and the uncertainties of the spectral indices.
% This program will, then organize the data and save an output file for each day of a
% given month provided that the input files are organized in that manner. 
% This form of data processing and organization helps to conveniently
% compute median parameters for each month as depicted in Figures 4, 5, and 6 in the
% manuscript.

% The output files will contain Data_holder.mat, which contains the following 
% parameter values:
% column 1: low-frequency spectral index
% column 2: uncertainty for low-frequency spectral index
% column 3: high-frequency spectral index 
% column 4: uncertainty for high-frequency spectral index
% column 5: spectral-break frequency 
% column 6: data validity indicator 
% column 7: mean X-CSEQ location for the spacecraft
% column 8: mean Y-CSEQ location for the spacecraft
% column 9: mean Z-CSEQ location for the spacecraft

% Before executing this program, the user will need to update the
% 'input_file_directory' and 'output_file_directory' variables below

clearvars

load('back_max.mat')
%Loads the matrix that contains the MAG instrument noise level 

input_file_directory = 'C:/Rosetta/MAT_MAG/2014/Sep/';
%Location of Rosetta RPCMAG data files in mat format

output_file_directory = 'C:/Rosetta/anlysis_folder/2014/Sep/';
%Location to save 'Data_holder' matrix as mat files

filesAndFolders = dir(input_file_directory);
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));
num_files = size(filesInDir);
num_files = num_files(1,1);

for file_i=1:num_files
    
    load(strcat(filesInDir(file_i).folder,'/',filesInDir(file_i).name));
    [tok,rem] = strtok(filesInDir(file_i).name,'_');

    dur = datenum(2014,11,28,0,13,40)-datenum(2014,11,28,0,0,0);
    dur2 = datenum(2014,11,28,0,0,30)-datenum(2014,11,28,0,0,0);
  
    start_date = Date_num(1,1);
    size_date_num = size(Date_num);
    time_segments = (Date_num(size_date_num(1,1),1)-start_date-dur)/dur2;
    time_segments = time_segments-mod(time_segments,1);
    
    Data_holder = zeros(time_segments,9); 

for time_incr=1:time_segments

    start_t = start_date+dur2*(time_incr-1);
    end_t = start_t+dur;%
    indices = find(Date_num >= start_t & Date_num <= end_t);

    if length(indices) >=16384 && length(indices) <=16400
    
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

        
        posx_now = mean(posx(indices,1));
        posy_now = mean(posy(indices,1));
        posz_now = mean(posz(indices,1));
        posx_diff = max(posx(indices,1))-min(posx(indices,1));
  
        Bx_now = Bx_new(1:16384,:);
        By_now = By_new(1:16384,:);
        Bz_now = Bz_new(1:16384,:);
      
        Data_holder(time_incr,:) = [0,0,0,0,0,0,posx_now,posy_now,posz_now];
  
        time_vec = time_vector_new(1:16384,:);
 
        pol_v = polyfit(time_vec,Bx_now,2);
        poly_func = polyval(pol_v,time_vec);
        x = Bx_now-poly_func;
        Fs = 20;
        N = length(x);
        if mod(N,2) ==1
            N = N-1;
        end
  
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx1 = (1/(Fs*N)) * abs(xdft).^2;
        pol_v = polyfit(time_vec,By_now,2);
        poly_func = polyval(pol_v,time_vec);
        x = By_now-poly_func;
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx2 = (1/(Fs*N)) * abs(xdft).^2;

        pol_v = polyfit(time_vec,Bz_now,2);
        poly_func = polyval(pol_v,time_vec);
        x = Bz_now-poly_func;
        xdft = fft(x);
        xdft = xdft(1:N/2+1);
        psdx3 = (1/(Fs*N)) * abs(xdft).^2;

        freq = 0:Fs/length(x):Fs/2;
        psdx = (psdx1+psdx2+psdx3);
        psdx(2:end-1) = 2*psdx(2:end-1);

        psdx = psdx';
        psdx =psdx(1,2:length(freq));
        freq = freq(1,2:length(freq));

        log_period = log10(freq);
        log_power = smooth(log10(psdx),101);
        log_power2 = log10(psdx);
        log_power = log_power';

        spec_index = zeros(1,5);
       
        st_f=-2.5;
        en_f=-2.0;

        for i=1:5

            track_f = find(log_period>st_f & log_period <=en_f);
            result3 = polyfit(log_period(1,track_f),log_power(1,track_f),1);
            spec_index(1,i) = result3(1,1);
  
            st_f = st_f+0.5;
            en_f = st_f+0.5;
  
        end

        if spec_index(1,5)>spec_index(1,4)
            high_range = -0.5;
        else
            high_range = 0.0;
        end

        low_range = -2.7;
        mid_f = -2.0;
        incr = (high_range-mid_f)/4.0;
        coarse_freq =zeros(1,4);
        valid_coarse = zeros(1,4);
        coarse_diff = zeros(1,4);

        for i=1:4
            
            coarse_freq(1,i) = mid_f;
            track_f_temp = find( log_period>-2.7 & log_period <mid_f); 
            track_f_temp1 = find( log_period>mid_f & log_period <high_range); 
            
            if length(track_f_temp)>2 && length(track_f_temp1)>2
                
                valid_coarse(1,i) = 1;
                result31 = polyfit(log_period(1,track_f_temp),log_power(1,track_f_temp),1);
                result32 = polyfit(log_period(1,track_f_temp1),log_power(1,track_f_temp1),1);

                temp_calculated_break = (result31(1,2)-result32(1,2))/(result32(1,1)-result31(1,1));
                coarse_diff(1,i) = abs(mid_f-temp_calculated_break);

            end
            
            mid_f = mid_f+incr;

        end

        valid_track0 = find(valid_coarse==1);
        coarse_diff=coarse_diff(1,valid_track0);

        coarse_freq = coarse_freq(1,valid_track0);
        find_max_cha= find(coarse_diff == min(coarse_diff));
        break_frequency_n = coarse_freq(1,find_max_cha);
        Amincor = islocalmin(coarse_diff);
        fin_amincor = find(Amincor==1);
            
        if length(fin_amincor)>=1
            
            fine_diff_min_cor = coarse_diff(1,fin_amincor);
            fine_freq_min_cor = coarse_freq(1,fin_amincor);
            wc= find(fine_diff_min_cor == min(fine_diff_min_cor));
            break_frequency_n  = fine_freq_min_cor(1,wc);
            
        end

        if spec_index(1,5)>spec_index(1,4)
            
            high_range = -0.5;
            
        else
            
            high_range = 0.0;
            
        end

        mid_range = break_frequency_n-0.5;
        incr = (high_range-mid_range)/20;
        fine_diff = zeros(1,21);
        fine_freq = zeros(1,21);
        valid_fine = zeros(1,21);

        for i=1:21
    
            fine_freq(1,i) = mid_range;
            track_f_temp = find( log_period>-2.7 & log_period <mid_range); 
            track_f_temp1 = find( log_period>mid_range & log_period <high_range); 
  
            if length(track_f_temp)>2 && length(track_f_temp1)>2
      
                valid_fine(1,i) = 1;
        
                result31 = polyfit(log_period(1,track_f_temp),log_power(1,track_f_temp),1);
                result32 = polyfit(log_period(1,track_f_temp1),log_power(1,track_f_temp1),1);

                temp_calculated_break = (result31(1,2)-result32(1,2))/(result32(1,1)-result31(1,1));
                fine_diff(1,i) = abs(mid_range-temp_calculated_break);
                
            end
            
            mid_range = mid_range+0.05;

        end
        
        valid_track = find(valid_fine==1);
        fine_diff=fine_diff(1,valid_track);
        fine_freq = fine_freq(1,valid_track);

        Amin = islocalmin(fine_diff);

        fin_amin = find(Amin==1);

        if length(fin_amin)>=1
    
            fine_diff_min = fine_diff(1,fin_amin);
            fine_freq_min = fine_freq(1,fin_amin);
            w= find(fine_diff_min == min(fine_diff_min));
            break_frequency_n2  = fine_freq_min(1,w);
            
        end
        
        if ~isempty(fin_amin)
        
            break_frequency_n = break_frequency_n2;
        
        end
        
        track_f_low = find(log_period>-2.7 & log_period <break_frequency_n & log_power2 > 2.0*back_max(1,find_max_cha));
        
        if (break_frequency_n >=-0.5)
        
            track_f_high = find( log_period>break_frequency_n & log_period <0.0 & log_power2 > 2.0*back_max(1,5)); 
        
        else
                 
            if(spec_index(1,5)>spec_index(1,4))
            
                track_f_high = find( log_period>break_frequency_n & log_period <-0.5 & log_power2 > 2.0*back_max(1,4)); 
 
            else
                
                track_f_high = find( log_period>break_frequency_n & log_period <0.0 & log_power2 > 2.0*back_max(1,5));

            end
    
        end

        if length(track_f_low)>1 && length(track_f_high)>1
    
            result3 = polyfit(log_period(1,track_f_low),log_power2(1,track_f_low),1);

            m_para = result3(1,1);
            c_para = result3(1,2);
            low_index = m_para;
            scales = log_period(1,track_f_low);
            Pow_spec_2 = log_power2(1,track_f_low);
            x_bar = mean(scales);
            D = sum((scales-x_bar).^2);
            d_i = Pow_spec_2-m_para.*scales-c_para;
            del_m_2 = sum(d_i.^2)/(D*(length(scales)-2));
            del_m = sqrt(del_m_2);
            low_index_err = del_m;
  
            result3 = polyfit(log_period(1,track_f_high),log_power2(1,track_f_high),1);

            m_para = result3(1,1);
            c_para = result3(1,2);
            high_index = m_para;
            scales = log_period(1,track_f_high);
            Pow_spec_2 = log_power2(1,track_f_high);
            x_bar = mean(scales);
            D = sum((scales-x_bar).^2);
            d_i = Pow_spec_2-m_para.*scales-c_para;
            del_m_2 = sum(d_i.^2)/(D*(length(scales)-2));
            del_m = sqrt(del_m_2);
            high_index_err = del_m;

            Data_holder(time_incr,1:6) = [low_index,low_index_err,high_index,high_index_err,break_frequency_n,1];
  
        end


    end
 
    clearvars -except back_max Date_num start_date time_segments Data_holder time_incr dur dur2 rem num_files file_i filesInDir output_file_directory Bx By Bz QF1 QF2 QF3 QF4 QF5 QF6 QF7 QF8 output_file_directory2 posx posy posz

end

save(strcat(output_file_directory,'spectrum_analysis',rem),'Data_holder','-v7.3');

clearvars -except back_max filesInDir file_i num_files output_file_directory 

end

clearvars
 