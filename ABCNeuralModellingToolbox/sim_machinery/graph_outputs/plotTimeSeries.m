function plotTimeSeries(TS_data,TS_sim,F,R,bestn,labelna)
if isempty(TS_data)
    TS_data = {zeros(size(TS_sim{1}))};
end
if isempty(TS_sim)
    TS_sim= {zeros(size(TS_data{1}))};
end
if nargin<5
    bestn = 1;
end
if isempty(labelna)
    labelna = 'Magnitude';
end
TS_data_n = TS_data{1};
subplot(2,1,1)
plot(F,TS_data_n{1}(1,:),'b'); hold on
xlabel('Time'); ylabel('X1')
subplot(2,1,2)
plot(F,TS_data_n{1}(2,:),'b'); hold on
xlabel('Time'); ylabel('X2')

for L = 1:length(TS_sim)
    TS_sim_n = TS_sim{L};
    if L == bestn
        lwid = 2;
    else
        lwid = 0.5;
    end
    subplot(2,1,1)
    plot(F,TS_sim_n{1}(1,:),'r:','LineWidth',lwid)
    subplot(2,1,2)
    plot(F,TS_sim_n{1}(2,:),'r:','LineWidth',lwid)
end