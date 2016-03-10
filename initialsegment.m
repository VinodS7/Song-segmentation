function [period_samples1,period_samples2,period_samples3, period_samples] = initialsegment(filePath, method, param)

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
n = 10;
x = buffer(x,44100*n,44100*n/2);
[~,i] = size(x);
period_samples = zeros([1,i]);
period_time = zeros([1,i]);
for k = 1:i
    overlap = param.windowSize - param.hopSize;
    X = spectrogram(x(:,k), param.windowSize, overlap, param.windowSize, fs);    
    X = abs(X);
    [~, HD, ~, ~, ~] = PfNmf(X, param.WD, [], [], [], param.rh, param.sparsity);
    
        [r1,~] = xcorr(HD(1,:),'coeff');
        [r2,~] = xcorr(HD(2,:),'coeff');
        [r3,~] = xcorr(HD(3,:),'coeff');
        [pks1,loc1] = findpeaks(r1(length(HD(1,:))+35:length(HD(1,:))+104));
        [pks2,loc2] = findpeaks(r2(length(HD(2,:))+35:length(HD(2,:))+104));
        [pks3,loc3] = findpeaks(r3(length(HD(3,:))+35:length(HD(3,:))+104));
        
        [~,I1] = max(pks1);
        period_samples1(k)= (loc1(I1)+34)*512;
        
        [~,I2] = max(pks2);
        period_samples2(k)= (loc2(I2)+34)*512;
        
        [~,I3] = max(pks3);
        period_samples3(k)= (loc3(I3)+34)*512;
        
         [r,~] = xcorr((HD(1,:)+HD(2,:)+HD(3,:)),'coeff');
        [pks,loc] = findpeaks(r(length(HD(1,:))+35:length(HD(1,:))+104));
        [~,I] = max(pks);
        period_samples(k)= (loc(I)+34)*512;
        
        figure(1)
        plot(r1)
        title('Autocorrelation for Hi-hat')
        
        figure(2)
        plot(r2)
        title('Autocorrelation for Bass drum')
        
        figure(3)
        plot(r3)
        title('Autocorrelation for Snare drum')
        
        figure(4)
        plot(r)
        title('Autocorrelation for combined')
    
        
end
    
    