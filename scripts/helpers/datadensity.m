%~~~~~~~~~~ Data Density ~~~~~~~~~~~~~~
function dd = datadensity(x,y,r)
% Jason Buenrostro, Stanford University
% Computes the data density (points/area) of scattered points

Ld = length(x);
dd = zeros(Ld,1);

for k=1:Ld
    dd(k) = sum( sqrt((x-x(k)).^2 + (y-y(k)).^2) < r );
end
area = pi*r^2;
%area = (4/3)*pi*r^3;
dd = dd/area;

return
