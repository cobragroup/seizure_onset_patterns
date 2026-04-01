function featurevector=calc_feat(d,fs)

%frequency
b0power= bandpower(d, fs, [0 0.5]);
b1power= bandpower(d, fs, [0.5 4]);
b2power= bandpower(d, fs, [4 12]);
b3power= bandpower(d, fs, [12 64]);
b4power= bandpower(d, fs, [64 fs/2-1]);

%spikes
q1 = quantile(d,0.05);
q3 = quantile(d,0.95);
lb = quantile(d,0.25) - 1.5*(quantile(d,0.75)-quantile(d,0.25));
ub = quantile(d,0.75) + 1.5*(quantile(d,0.75)-quantile(d,0.25));

alphdiff= q3-q1;
spikeabs = length(find(or(d>ub,d<lb)));

%general mean and variance
sigmean=mean(d,'omitnan');
sigvar=var(d,'omitnan');

%critical slowing
acf = autocorr(d,'NumLags',round(fs/1000*5)); %lag = 5 ms
autocorrel=acf(end);
linelen=sum(sqrt(diff(d).^2+(1/fs)^2),'omitnan');

featurevector=[ b0power, b1power, b2power, b3power, b4power, alphdiff, spikeabs, sigmean, sigvar, autocorrel, linelen];
