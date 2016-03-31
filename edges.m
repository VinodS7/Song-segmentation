function [edge,indices] = edges(fileName,bar_length)

[hh,bd,sd] = NmfDrum(fileName);
[y,Fs] = audioread(fileName);

HH = zeros([1 length(y)]);
BD = zeros([1 length(y)]);
SD = zeros([1 length(y)]);
ref1 = zeros([1 length(y)]);
ref2 = zeros([1 length(y)]);
HH(round(hh*Fs)) = 1;
BD(round(bd*Fs)) = 1;
SD(round(sd*Fs)) = 1;

[~,I1] = findpeaks(HH);
[~,I2] = findpeaks(BD);
[~,I3] = findpeaks(SD);
m = min(I1(1),I2(1));
m = min(m,I3(1));
I4 = cat(2,I1,I2);
I4 = cat(2,I4,I3);

p = bar_length*Fs;
k = m:p:length(y);
[~,I] = pdist2(I4',k','euclidean','Smallest',1);
    edge = I4(I);
   ref1(I4) = 1;
  [~,~,indices] = histcounts(ref1,edge);
 
    
    