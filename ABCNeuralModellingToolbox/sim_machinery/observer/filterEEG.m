function fdata = filterEEG(Data,srate,locutoff,hicutoff,filtorder)
for n = 1:size(Data,1)
    data = Data(n,:);
t=0:1/srate:length(data)/srate-1/srate;

nyq            = srate*0.5;  % Nyquist frequency
MINFREQ        = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter 
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
m=[0       0                      1            1            0                      0]; 
filtwts = firls(filtorder,f,m); % get FIR filter coefficients

smoothdata = zeros(length(data),1);
fdata(n,:) = filtfilt(filtwts,1,data);
end