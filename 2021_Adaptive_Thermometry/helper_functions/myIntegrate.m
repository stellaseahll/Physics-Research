function Iy = myIntegrate(x,y,Iy0)
if nargin ==2
    Iy0 = 0;
end
Iy = Iy0;
for i = 2:length(x)
    Iy(i) = trapz(x(1:i),y(1:i))+Iy(1);
end