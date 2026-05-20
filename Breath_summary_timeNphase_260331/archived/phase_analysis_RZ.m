d = designfilt("lowpassfir", ...
    PassbandFrequency=0.2,StopbandFrequency=0.25, ...
    PassbandRipple=1,StopbandAttenuation=60, ...
    DesignMethod="equiripple");
breath_filt = filtfilt(d,breath);

%%
y = hilbert(breath_filt);

%%
figure;
subplot(3,1,1);
plot(t_breath, breath);

subplot(3,1,2);
plot(t_breath, breath_filt);

subplot(3,1,3);
plot(t_breath, angle(y));

%%
phi = angle(y);


%%
%figure;

for i = 1:numel(ca_onsets)
    ca_phi(i) = phi(ca_onsets(i));
end

