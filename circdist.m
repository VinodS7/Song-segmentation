function distance = circdist(file1,file2,segments)
d = zeros([1 segments]);
for k = 1:segments
    file2  = circshift(file2,[0,k]);
    for i = 1:segments
    d(k) = d(k) + abs(file1(i)-file2(i));
    end
   
end
distance = min(d);