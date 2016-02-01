function [period_samples,period_time,i] = Shadysegment(filePath, method, param)

if nargin == 2
    load DefaultSetting.mat
elseif nargin == 1
    load DefaultSetting.mat
    method = 'PfNmf'; %by default, use PfNmf
end

fprintf('Selected method is %s\n', method);

%//load file
[x, fs] = audioread(filePath); 
x = mean(x,2); %down-mixing   
x = resample(x, 44100, fs); %sample rate consistency
fs = 44100;
L = length(x);
p = nextpow2(fs*5);
x = buffer(x,2^p,2^(p-1));
i = round(L/2^p);
period_samples = zeros([1,i]);
period_time = zeros([1,i]);
for k = 1:i
    overlap = param.windowSize - param.hopSize;
    X = spectrogram(x(:,k), param.windowSize, overlap, param.windowSize, fs);    
    X = abs(X);
    [~, HD, ~, ~, ~] = PfNmf(X, param.WD, [], [], [], param.rh, param.sparsity);
    
        [r,~] = xcorr((HD(1,:)+HD(2,:)+HD(3,:)),'coeff');
        [pks,loc] = findpeaks(r(length(HD(2,:))+1:length(r)));
        [~,I] = max(pks);
        period_samples(k)= loc(I)*512;
        period_time(k) = 120/(loc(I)*512/44100);
    
end
    
    