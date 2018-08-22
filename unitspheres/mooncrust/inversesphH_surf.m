function F = inversesphH_surf(p)
load LOLA_PA.mat
pcut = 60; lsiz = (pcut+1)*(pcut+2)/2;  
coef = LOLA_PA_310(1:lsiz,3:4);



ltp = length(p(:,1));
gpsiz = min(round(1e7/(pcut+1)^2),ltp);
t0 = 0; t1 = 0; 
ng = ceil(ltp/gpsiz);
F = zeros(ltp,1);
for i = 1:ng
    t1 = min(t1 + gpsiz,ltp);
    F(t0+1:t1) = subgroupinvsphH(p(t0+1:t1,:),coef,pcut);
    t0 = t1;
end


F = F/1e3;
%max(F)