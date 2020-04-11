% Supplementary material to Suranga Ruhunusiri, "Plasma Turbulence at Comet 67P/Churyumov-Gerasimenko: Rosetta Observations", 
% submittd to Journal of Geophysical Research Space Physics, 2020

% This program has been tested with MATLAB R2019a.

% This program computes a power spectrum for magnetic field fluctuations, determine 
% the spectral-break frequency, low- and high-frequency spectral indices,
% and plot the power spectrum, log10(Power spectral density[nT^2/Hz]) vs. log10(frequency[Hz]). 
% On the plot, the location of the spectral
% break frequency is displayed as a vertical line centered on the spectral
% break frequency. Two lines are also plotted on either side of the
% spectral break frequency whose slope is equal to the low- and
% high-frequency spectral indices respectively.
%

% Prior to executing this code, the user will need to first load RPCMAG
% data into MATLAB workspace by first running RPCMAG_data_reader.m and edit INPUT 1 below.

%INPUT1:SDate
%Start time for computing power spectrum of magnetic field fluctuations in
%Year Month Day Hour Minutes Seconds format

SDate = [2014,11,22,3,10,0];


dur = datenum(2014,11,28,0,17,4)-datenum(2014,11,28,0,0,0);
row_select = 1;
start_date = datenum(SDate(row_select,1),SDate(row_select,2),SDate(row_select,3),SDate(row_select,4),SDate(row_select,5),SDate(row_select,6));

start_t = start_date;
end_t = start_t+dur;

indices = find(Date_num > start_t & Date_num < end_t);

Bx_now = Bx(indices,1);
By_now = By(indices,1);
Bz_now = Bz(indices,1);

x = Bx_now;
Fs = 20;
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx1 = (1/(Fs*N)) * abs(xdft).^2;

x = By_now;
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx2 = (1/(Fs*N)) * abs(xdft).^2;

x = Bz_now;
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx3 = (1/(Fs*N)) * abs(xdft).^2;

freq = 0:Fs/length(x):Fs/2;
psdx = (psdx1+psdx2+psdx3);
psdx(2:end-1) = 2*psdx(2:end-1);
smpsdx=smooth(log10(psdx),41);

figure
hold on;
plot(log10(freq),log10(psdx))
plot(log10(freq),smpsdx)
grid on
title('Magnetic field power spectrum')
xlabel('log(Frequency)')
ylabel('log(PSD)')
freqo= freq;
psdxo=psdx';
psdx = psdx';
psdx =psdx(1,2:length(freq));
freq = freq(1,2:length(freq));

log_period = log10(freq);
log_power = smooth(log10(psdx),41);
log_power = log_power';
n_len = (-0.5 +2.4)/0.2;
n_len = n_len-mod(n_len,1);
spec_index = zeros(1,n_len-1);
spec_index_freq = zeros(1,n_len-1);

for i=1:n_len-1

  st_f = -2.4 + 0.2*i;
  en_f = st_f + 0.2;
  track_f = find(log_period>st_f & log_period < en_f);
  result3 = polyfit(log_period(1,track_f),log_power(1,track_f),1);
  m_para = result3(1,1);
  spec_index(1,i) = m_para;
  spec_index_freq(1,i) = st_f + 0.1;
end

index_change = zeros(1,n_len-1);
index_change_freq = zeros(1,n_len-1);
valid_index = zeros(1,n_len-1);
index_hol = zeros(2,n_len-1);
 
for i=2:n_len-2
if (spec_index(1,i+1) < 0.0 && spec_index(1,i) < 0 && spec_index(1,i-1) < 0.0 &&  spec_index(1,i+1) < spec_index(1,i-1) && spec_index(1,i+1) < spec_index(1,i) && spec_index(1,i) < spec_index(1,i-1))
   
  index_change(1,i) = (abs((spec_index(1,i+1)-spec_index(1,i-1))));
end
  index_change_freq(1,i) =spec_index_freq(1,i); 
  index_hol(1,i) = spec_index(1,i-1);
  index_hol(2,1) =spec_index(1,i+1);
  if (spec_index(1,i+1) < 0.0 && spec_index(1,i) < 0 && spec_index(1,i-1) < 0.0 &&  spec_index(1,i+1) < spec_index(1,i-1) && spec_index(1,i+1) < spec_index(1,i) && spec_index(1,i) < spec_index(1,i-1))
    valid_index(1,i) = 1.0;
  end
  
end

 ind_ch_fre_tracker = find(index_change_freq < 0.0 & index_change_freq > -2.0 & valid_index == 1.0);

 sel_index_change_freq1=index_change_freq(1,ind_ch_fre_tracker);
 sel_index_change = index_change(1,ind_ch_fre_tracker);
 bb = find(sel_index_change == max(sel_index_change));

bb1 = bb;

log_period = log10(freq);
log_power = smooth(log10(psdx),41);
log_power = log_power';
n_len = (-0.5 +2.4)/0.1;
n_len = n_len-mod(n_len,1);
spec_index = zeros(1,n_len-1);
spec_index_freq = zeros(1,n_len-1);

for i=1:n_len-1

  st_f = -2.4 + 0.1*i;
  en_f = st_f + 0.1;
  track_f = find(log_period>st_f & log_period < en_f);
 
  result3 = polyfit(log_period(1,track_f),log_power(1,track_f),1);
  m_para = result3(1,1);
  spec_index(1,i) = m_para;
  spec_index_freq(1,i) = st_f + 0.05;
end

index_change = zeros(1,n_len-1);
index_change_freq = zeros(1,n_len-1);
valid_index = zeros(1,n_len-1);
index_hol = zeros(2,n_len-1);

for i=2:n_len-2
if (spec_index(1,i+1) < 0.0 && spec_index(1,i) < 0 && spec_index(1,i-1) < 0.0 &&  spec_index(1,i+1) < spec_index(1,i-1) && spec_index(1,i+1) < spec_index(1,i) && spec_index(1,i) < spec_index(1,i-1))
   
  index_change(1,i) = (abs((spec_index(1,i+1)-spec_index(1,i-1))));
end
  index_change_freq(1,i) =spec_index_freq(1,i); 
  index_hol(1,i) = spec_index(1,i-1);
  index_hol(2,1) =spec_index(1,i+1);
  if (spec_index(1,i+1) < 0.0 && spec_index(1,i) < 0 && spec_index(1,i-1) < 0.0 &&  spec_index(1,i+1) < spec_index(1,i-1) && spec_index(1,i+1) < spec_index(1,i) && spec_index(1,i) < spec_index(1,i-1))
    valid_index(1,i) = 1.0;
  end
  
end

 ind_ch_fre_tracker = find(index_change_freq < 0.0 & index_change_freq > -2.0 & valid_index == 1.0);

 sel_index_change_freq2=index_change_freq(1,ind_ch_fre_tracker);
 sel_index_change = index_change(1,ind_ch_fre_tracker);
 index_hol_sel = index_hol(1,ind_ch_fre_tracker);
 bb = find(sel_index_change == max(sel_index_change));

 bb2 = bb;

if length(bb1) >0
auto_spec_break1 =  10.0^sel_index_change_freq1(1,bb1);
else
  auto_spec_break1=bb1;  
end
if length(bb2)> 0
auto_spec_break2 =  10.0^sel_index_change_freq2(1,bb2);
else
    auto_spec_break2=bb2; 
end

if auto_spec_break1<=auto_spec_break2 && length(bb1)>0
   auto_spec_break = auto_spec_break1;
     
end
if auto_spec_break2<=auto_spec_break1 && length(bb2)>0
   auto_spec_break = auto_spec_break2;
end
if length(bb1) ==0
   
     auto_spec_break =  auto_spec_break2;
    
end
if length(bb2)==0
   
     auto_spec_break =  auto_spec_break1;
end

if length(bb1)==0 && length(bb2)==0
    
log_period = log10(freq);
log_power = smooth(log10(psdx),41);
log_power = log_power';
n_len = (-0.5 +2.4)/0.2;
n_len = n_len-mod(n_len,1);
spec_index = zeros(1,n_len-1);
spec_index_freq = zeros(1,n_len-1);

for i=1:n_len-1

  st_f = -2.4 + 0.2*i;
  en_f = st_f + 0.2;
  track_f = find(log_period>st_f & log_period < en_f);
 
  result3 = polyfit(log_period(1,track_f),log_power(1,track_f),1);
  m_para = result3(1,1);
  spec_index(1,i) = m_para;
  spec_index_freq(1,i) = st_f + 0.1;
end

index_change = zeros(1,n_len-1);
index_change_freq = zeros(1,n_len-1);
valid_index = zeros(1,n_len-1);
index_hol = zeros(2,n_len-1);
 
for i=2:n_len-2
if (spec_index(1,i+1) < 0.0 && spec_index(1,i-1) > 0)
   
  index_change(1,i) = (abs((spec_index(1,i+1)-spec_index(1,i-1))));
end
  index_change_freq(1,i) =spec_index_freq(1,i); 
  index_hol(1,i) = spec_index(1,i-1);
  index_hol(2,1) =spec_index(1,i+1);
  if (spec_index(1,i+1) < 0.0 && spec_index(1,i-1) > 0)
    valid_index(1,i) = 1.0;
  end
  
end

 ind_ch_fre_tracker = find(index_change_freq < 0.0 & index_change_freq > -2.0 & valid_index == 1.0);

 sel_index_change_freq1=index_change_freq(1,ind_ch_fre_tracker);
 sel_index_change = index_change(1,ind_ch_fre_tracker);
 bb = find(sel_index_change == max(sel_index_change));

bb1 = bb;

log_period = log10(freq);
log_power = smooth(log10(psdx),41);
log_power = log_power';
n_len = (-0.5 +2.4)/0.1;
n_len = n_len-mod(n_len,1);
spec_index = zeros(1,n_len-1);
spec_index_freq = zeros(1,n_len-1);

for i=1:n_len-1

  st_f = -2.4 + 0.1*i;
  en_f = st_f + 0.1;
  track_f = find(log_period>st_f & log_period < en_f);
  
  result3 = polyfit(log_period(1,track_f),log_power(1,track_f),1);
  m_para = result3(1,1);
  spec_index(1,i) = m_para;
  spec_index_freq(1,i) = st_f + 0.05;
end

index_change = zeros(1,n_len-1);
index_change_freq = zeros(1,n_len-1);
valid_index = zeros(1,n_len-1);
index_hol = zeros(2,n_len-1);

for i=2:n_len-2
if (spec_index(1,i+1) < 0.0 && spec_index(1,i-1) > 0)
   
  index_change(1,i) = (abs((spec_index(1,i+1)-spec_index(1,i-1))));
end
  index_change_freq(1,i) =spec_index_freq(1,i); 
  index_hol(1,i) = spec_index(1,i-1);
  index_hol(2,1) =spec_index(1,i+1);
  if (spec_index(1,i+1) < 0.0 && spec_index(1,i-1) > 0)
    valid_index(1,i) = 1.0;
  end
  
end

 ind_ch_fre_tracker = find(index_change_freq < 0.0 & index_change_freq > -2.0 & valid_index == 1.0);

 sel_index_change_freq2=index_change_freq(1,ind_ch_fre_tracker);
 sel_index_change = index_change(1,ind_ch_fre_tracker);
 index_hol_sel = index_hol(1,ind_ch_fre_tracker);
 bb = find(sel_index_change == max(sel_index_change));

 bb2 = bb;

if length(bb1) >0
auto_spec_break1 =  10.0^sel_index_change_freq1(1,bb1);
else
  auto_spec_break1=bb1;  
end
if length(bb2)> 0
auto_spec_break2 =  10.0^sel_index_change_freq2(1,bb2);
else
    auto_spec_break2=bb2; 
end

if auto_spec_break1<=auto_spec_break2 && length(bb1)>0
   auto_spec_break = auto_spec_break1;
     
end
if auto_spec_break2<=auto_spec_break1 && length(bb2)>0
   auto_spec_break = auto_spec_break2;
end
if length(bb1) ==0
   
     auto_spec_break =  auto_spec_break2;
    
end
if length(bb2)==0
   
     auto_spec_break =  auto_spec_break1;
end
    
    
    
end


xx = [log10(auto_spec_break),log10(auto_spec_break)];
yy = [min(log10(psdx)),max(log10(psdx))];

plot(xx,yy)


scales = ((freq));
aa = find(scales <= auto_spec_break & scales >= 0.002);


 scales = log10(scales(1,aa));
 Pow_spec_2 = log_power;

 Pow_spec_2 = (Pow_spec_2(1,aa));

 result3 = polyfit(scales,Pow_spec_2,1);

low_fre_spec_index= result3(1,1);

m_para = result3(1,1);
c_para = result3(1,2);

 xx = scales;

 yy = m_para*xx + c_para + 0.5;
plot(xx,yy)

scales = ((freq));
aa = find(scales <= 1.0 & scales >= auto_spec_break);

 scales = log10(scales(1,aa));
 Pow_spec_2 = log_power;

 Pow_spec_2 = (Pow_spec_2(1,aa));

 result3 = polyfit(scales,Pow_spec_2,1);
 
high_fre_spec_index= result3(1,1);
 
m_para = result3(1,1);
c_para = result3(1,2);

 xx = scales;


yy = m_para*xx + c_para + 0.5;
plot(xx,yy)

Spec_matx = [auto_spec_break,low_fre_spec_index,high_fre_spec_index];

clearvars -except Date_num posx posy posz Bx By Bz QF1 QF2 QF3 QF4 QF5 QF6 QF7 QF8;
