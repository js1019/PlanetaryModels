function d = sphdist(p)
   R = 1737.0-1497;
   d = p(:,1).^2+p(:,2).^2+p(:,3).^2;
   d = d - R^2;
return