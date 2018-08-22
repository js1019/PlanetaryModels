function F = subgroupinvsphH(p,coef,pcut)

[azimuth,elevation,~] = cart2sph(p(:,1),p(:,2),p(:,3));
%min(elevation)
%dirs = [azimuth(:)';pi/2-elevation(:)']';
dirs(:,1) = pi/2-azimuth(:);
dirs(:,2) = pi/2-elevation(:);
F_N = zeros((pcut+1)^2,1);
if pcut == 0
    F_N = coef(:,1);
else
for i = 1:pcut
    l = i^2;
    F_N(l+1:l+i) = coef(((i+1)*(i+2))/2:-1:((i+1)*i)/2+2,2); 
    F_N(l+i+1:(i+1)^2) = coef((i+1)*i/2+1:(i+2)*(i+1)/2,1);
    %F_N(l+1:l+i) = -coef(((i+1)*(i+2))/2:-1:((i+1)*i)/2+2,1);
    %F_N(l+i+1:(i+1)^2) = coef((i+1)*i/2+1:(i+2)*(i+1)/2,1);
end
%F_N(:)
%F_N(4)
end
F = inverseSHT(F_N, dirs, 'real');
return