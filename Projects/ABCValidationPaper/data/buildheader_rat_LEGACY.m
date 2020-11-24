function R = buildheader_rat()
R.pipestamp = 'rat_020317';
R.subnames = {{'C1','C2','C3','C4','C5','C6','C8'}...
    {'L4','L6','L13','L18','L19','L20','L21','L22','L23'}};
% EXCLUDED RAT C7 - IpsiEEG
R.bbounds = [5 10; 14 24; 25 40]; % Bands of interest
R.bandnames = {'Alpha','Low Beta','High beta'};
R.bandnames_nospace = {'alpha','low_beta','high_beta'};
R.sourcenames = {'fEEG','GP','STR','STN'};
R.condnames = {'control','lesion'};
R.FOI = [4 100]; % Frequencies of interest
R.pp.ds = 250; % sample rate to downsample to
R.pp.hpfilt = [1 2]; % highpass [stop pass]
R.pp.lpfilt = [100 95]; % lowpass [stop pass]
R.montage = 'monopolar';
R.montname = {'monopolar','bipolar','globalaverage','globaldistal','localdistal','bipolar_fixed','bipolarEEG'}; %,'princomp'};
% number of channels per source per subject nreps{cond}(sub srcloc)
R.nreps{1} = [
    1	8	8	4
    1	9	6	4
    1	9	6	5
    1	8	7	5
    1	9	5	4
    1	7	6	3
    1	9	3	4
    1	9	5	4
    ];
R.nreps{2} = [
    1	9	6	4
    1	9	6	4
    1	10	4	2
    1	9	4	4
    1	11	4	2
    1	10	5	5
    1	9	6	2
    1	8	7	3
    1	9	6	3
    ];
if  strcmp(getenv('COMPUTERNAME'), 'FREE') == 1
    R.datapath = 'C:\home\data\ratdata_050816\';
    R.analysispath = 'C:\Users\twest\Documents\Work\PhD\LitvakProject\rat_data\pipeline\';
else
    R.datapath = 'C:\home\data\TimExtracts310516\';
    R.analysispath = 'C:\Users\Tim\Documents\Work\LitvakProject\rat_data\pipeline\';
    
end

