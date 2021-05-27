function [radialLengthList,normRadialLengthList,meanNormRadialLength,normRadialLengthStandDev] = findRadialLengths(stateLattice)

C = regionprops(stateLattice,'Centroid');
centroid = C.Centroid;

B = bwboundaries(stateLattice,8,'noholes');
for k = 1:length(B)
   boundary = B{k};
end

[boundaryLength,~] = size(boundary);
radialLengthList = zeros(boundaryLength);
for k = 1:(boundaryLength-1)
    x1 = [boundary(k,1) boundary(k,2)];
    radialLengthList(k) = norm(x1-centroid);
end

maxRadialLength = max(radialLengthList);

normRadialLengthList = radialLengthList/maxRadialLength;

meanNormRadialLength = mean(normRadialLengthList);

normRadialLengthStandDev = std(normRadialLengthList);

