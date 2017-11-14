function [ BinMean, BinRadius, edges,BinIndex ] = BinAverage(dR, VariableSer, RadiusSer )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
N=int16(max(RadiusSer)/dR); %Number of bins
binrange=0:dR:max(RadiusSer);
[N,edges,BinIndex] = histcounts(RadiusSer,[binrange]);
for i = 1:length(edges)-1
   BinRadius(i)=(edges(i)+edges(i+1))/2; 
end

for n = 1:length(BinRadius)
    k=find(n==BinIndex);
    BinMean(n) = mean(VariableSer(k));
end

end

