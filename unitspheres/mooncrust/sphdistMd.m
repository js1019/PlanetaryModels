function d = sphdistMd(p)
   R = 1737.0;
   %F = inversesphH(p);
   %F = inversesphH_cmi(p);
   %F = inversesphH_surf(p);
   %F = max((-34+0.6)/sqrt(4*pi),F);
   %min((F))*sqrt(4*pi)+34
   %max((F))*sqrt(4*pi)+34
   
   Fs = inversesphH_surf(p);
   min(Fs)*sqrt(4*pi)
   max(Fs)*sqrt(4*pi)
   
   %sum(F)/length(F)*1e-3
   r = p(:,1).^2+p(:,2).^2+p(:,3).^2;
   %Fd = R+Fs-real(F)*sqrt(4*pi)-34;
   Fd = R + Fs;
   d = sqrt(r) - Fd;
return