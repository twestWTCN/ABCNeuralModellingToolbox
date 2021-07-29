function plotRawTimeSeries(R,data)
tvec_obs = R.IntP.tvec;
tvec_obs(:,2:round(R.obs.brn*(1/R.IntP.dt))) = [];
R.IntP.tvec_obs = tvec_obs;
subplot(2,1,1)
plot(repmat(R.IntP.tvec_obs,size(data,1),1)',data');
xlabel('Time (s)'); ylabel('Amplitude')
subplot(2,1,2)
plot(repmat(R.IntP.tvec_obs,size(data,1),1)',data'); xlim([5 6])
xlabel('Time (s)'); ylabel('Amplitude')
legend(R.chsim_name)
