function res = fundistant(p,r,cost)
global RI epsl;
% r = a(1-etsl(a) (cos^2 \theta - 1/3))

ep = interp1(RI,epsl,p,'pchip');
res = r- p*(1-ep*(cost-1/3));

end