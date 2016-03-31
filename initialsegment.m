function [period_samples1,period_samples2,period_samples3, period_samples] = initialsegment(filePath,predicted_period)

period_samples1 = 0;
period_samples2 = 0;
period_samples3 = 0;
period_samples = 0;

tm = round(predicted_period*44100/512);

load DefaultSetting.mat

%//load file
[x, fs] = audioread(filePath); 
x = mean(x,2); %down-mixing   
x = resample(x, 44100, fs); %sample rate consistency
fs = 44100;
L = length(x);
n = 5;
x = buffer(x,44100*n,44100*n/2);
[~,i] = size(x);


for k = 1:i
    overlap = param.windowSize - param.hopSize;
    X = spectrogram(x(:,k), param.windowSize, overlap, param.windowSize, fs);    
    X = abs(X);
    [~, HD, ~, ~, ~] = PfNmf(X, param.WD, [], [], [], param.rh, param.sparsity);

        [r1,~] = xcorr(HD(1,:));
        [r2,~] = xcorr(HD(2,:));
        [r3,~] = xcorr(HD(3,:));
        hd1 = HD(1,:);
        hd2 = HD(2,:);
        hd3 = HD(3,:);
        marker1 = zeros([1 length(r1)]);
        marker2 = zeros([1 length(r1)]);
        marker3 = zeros([1 length(r1)]);
        marker4 = zeros([1 length(r1)]);
        marker5 = zeros([1 length(r1)]);
        marker6 = zeros([1 length(r1)]);
        
        blksize = 10;
        hd1 = padarray(hd1,[0 blksize],'pre');
        hd3 = padarray(hd3,[0 blksize],'pre');
        hd2 = padarray(hd2,[0 blksize],'pre');
  
        for s = blksize+1:length(hd1)
            nf1(s-blksize) = median(hd1(s-blksize:s));
        nf2(s-blksize) = median(hd2(s-blksize:s));
        nf3(s-blksize) = median(hd3(s-blksize:s));
        end
        
        nf1 = HD(1,:)-nf1;
        nf2 = HD(2,:)-nf2;
        nf3 = HD(3,:)-nf3;
        nf1(nf1<0) = 0;
        nf2(nf2<0) = 0;
        nf3(nf3<0) = 0;
        [r4,~] = xcorr(nf1);
        [r5,~] = xcorr(nf2);
        [r6,~] = xcorr(nf3);
        [pks1,loc1] = findpeaks(r1(length(HD(1,:))+86:length(HD(1,:))+258));
        [pks2,loc2] = findpeaks(r2(length(HD(2,:))+86:length(HD(2,:))+258));
        [pks3,loc3] = findpeaks(r3(length(HD(3,:))+86:length(HD(3,:))+258));
        
        [~,I1] = max(pks1);
        period_samples1(k)= (loc1(I1)+85)*512;
        
        [~,I2] = max(pks2);
        period_samples2(k)= (loc2(I2)+85)*512;
        
        [~,I3] = max(pks3);
        period_samples3(k)= (loc3(I3)+85)*512;

        marker1(tm+length(hd1)) = max(r1);
        marker2(tm+length(hd1)) = max(r2);
        marker3(tm+length(hd1)) = max(r3);
        marker4(tm+length(hd1)) = max(r4);
        marker5(tm+length(hd1)) = max(r5);
        marker6(tm+length(hd1)) = max(r6);
        t = (-2.5:5/length(HD(2,:)):2.5-5/length(HD(2,:)))+2.5*(k-1);
       close all 
      
       h1 = figure(1);
       subplot(321)
        hold on
        plot(marker1,'r')
        plot(r1)
        hold off
        title('Autocorrelation for Hi-hat')
        subplot(322)
        hold on
        plot(marker4,'r')
        plot(r4)
        title('Autocorrelation for Hi-Hat with median')
        hold off
        subplot(323)
        hold on
        plot(marker2,'r')
        plot(r2)
        title('Autocorrelation for Bass drum')
        hold off
        subplot(324)
        hold on
        plot(marker5,'r')
        plot(r5)
        title('Autocorrelation for Bass Drum with median')
        hold off
        subplot(325)
        hold on
        plot(marker3,'r')
        plot(r3)
        title('Autocorrelation for Snare drum')
        hold off
        subplot(326)
        hold on
        plot(marker6,'r')
        plot(r6)
        title('Autocorrelation for Snare with median')
        hold off
        
        h2 = figure(2);
        subplot(321)
        plot(t,HD(1,:))
        title('Novelty function Hi-Hat')
        xlabel('Time in Seconds')
        subplot(322)
        plot(t,nf1)
        title('Median LP novelty')
        xlabel('Time in Seconds')
        subplot(323)
        plot(t,HD(2,:))
        title('Novelty function Bass Drum')
        xlabel('Time in Seconds')
        subplot(324)
        plot(t,nf2)
        title('Median LP novelty')
        xlabel('Time in Seconds')
        subplot(325)
        plot(t,HD(3,:))
        title('Novelty function Snare Drum')
        xlabel('Time in Seconds')
        subplot(326)
        plot(t,nf3)
        title('Median LP novelty')
        xlabel('Time in Seconds')
        
        saveas(h1,sprintf('Autocorrelation%d.jpg',k));
        saveas(h2,sprintf('Novelty%d.jpg',k));
end
    
    