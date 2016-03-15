function [period_samples1,period_samples2,period_samples3, period_samples] = initialsegment(filePath, method, param)
period_samples1 = 0;
period_samples2 = 0;
period_samples3 = 0;
period_samples = 0;
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
n = 10;
x = buffer(x,44100*n);
[~,i] = size(x);
period_samples = zeros([1,i]);
for k = 1:i
    overlap = param.windowSize - param.hopSize;
    X = spectrogram(x(:,k), param.windowSize, overlap, param.windowSize, fs);    
    X = abs(X);
    [~, HD, ~, ~, ~] = PfNmf(X, param.WD, [], [], [], param.rh, param.sparsity);
    hd(1,:) = (HD(1,:) - min(HD(1,:)))/(max(HD(1,:)-min(HD(1,:))));
    hd(2,:) = (HD(2,:) - min(HD(2,:)))/(max(HD(2,:)-min(HD(2,:))));
    hd(3,:) = (HD(3,:) - min(HD(3,:)))/(max(HD(3,:)-min(HD(3,:))));
    for n = 1:length(hd(2,:))
        for m = 1:3
            if hd(m,n) < 0.50
                hd(m,n) = 0;
            end
        end
    end
        [r1,~] = xcorr(hd(1,:));
        [r2,~] = xcorr(hd(2,:));
        [r3,~] = xcorr(hd(3,:));
        hd1 = HD(1,:);
        hd2 = HD(2,:);
        hd3 = HD(3,:);
        blksize = 16;
        hd1 = padarray(hd1,[0 blksize],'pre');
        hd3 = padarray(hd3,[0 blksize],'pre');
        hd2 = padarray(hd2,[0 blksize],'pre');
        for s = blksize+1:length(hd1)
            nf1(s-blksize) = median(hd1(s-blksize:s));
        nf2(s-blksize) = median(hd2(s-blksize:s));
        nf3(s-blksize) = median(hd3(s-blksize:s));
        end
%         [pks1,loc1] = findpeaks(r1(length(HD(1,:))+35:length(HD(1,:))+104));
%         [pks2,loc2] = findpeaks(r2(length(HD(2,:))+35:length(HD(2,:))+104));
%         [pks3,loc3] = findpeaks(r3(length(HD(3,:))+35:length(HD(3,:))+104));
%         
%         [~,I1] = max(pks1);
%         period_samples1(k)= (loc1(I1)+34)*512;
%         
%         [~,I2] = max(pks2);
%         period_samples2(k)= (loc2(I2)+34)*512;
%         
%         [~,I3] = max(pks3);
%         period_samples3(k)= (loc3(I3)+34)*512;
%         
         [r,~] = xcorr((hd(1,:)+hd(2,:)+hd(3,:))/3);
%         [pks,loc] = findpeaks(r(length(HD(1,:))+35:length(HD(1,:))+104));
%         [~,I] = max(pks);
%         period_samples(k)= (loc(I)+34)*512;
%         
        t = (0:10/length(hd(2,:)):10-10/length(hd(2,:)))+10*(k-1);
        
       h1 = figure(1);
       
        subplot(331)
        plot(r1)
        title('Autocorrelation for Hi-hat')
        subplot(332)
        plot(t,HD(1,:))
        title('Novelty function Hi-Hat')
        xlabel('Time in Seconds')
        subplot(333)
        plot(t,nf1)
        title('Median LP novelty')
        xlabel('Time in Seconds')
        
        
        subplot(334)
        plot(r2)
        title('Autocorrelation for Bass drum')
        subplot(335)
        plot(t,HD(2,:))
        title('Novelty function Bass Drum')
        xlabel('Time in Seconds')
        subplot(336)
        plot(t,nf2)
        title('Median LP novelty')
        xlabel('Time in Seconds')
        
        subplot(337)
        plot(r3)
        title('Autocorrelation for Snare drum')
        subplot(338)
        plot(t,HD(3,:))
        title('Novelty function Snare Drum')
        xlabel('Time in Seconds')
        subplot(339)
        plot(t,nf3)
        title('Median LP novelty')
        xlabel('Time in Seconds')
%         subplot(427)
%         plot(r)
%         title('Autocorrelation for combined')
%         subplot(428)
%         plot(t,(hd(1,:)+hd(2,:)+hd(3,:))/3)
%         title('Novelty function Combined')
%         xlabel('Time in Seconds')
        saveas(h1,sprintf('Graph%d.jpg',k));
end
    
    