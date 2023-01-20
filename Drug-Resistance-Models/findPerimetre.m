function perimetre = findPerimetre(stateLattice,showPerimetre)

B = bwboundaries(stateLattice,8,'noholes');
for k = 1:length(B)
   boundary = B{k};
   if showPerimetre == true
       plot(boundary(:,1), boundary(:,2), 'bl', 'LineWidth', 1)
   end
end
perimetre = 0;
[boundaryLength,~] = size(boundary);
for k = 1:(boundaryLength-1)
    x1 = [boundary(k,1) boundary(k,2)];
    x2 = [boundary(k+1,1) boundary(k+1,2)];
    d = norm(x1-x2);
    perimetre = perimetre + d;
end