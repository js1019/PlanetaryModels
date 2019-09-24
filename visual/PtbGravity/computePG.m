% compute the values
G = 6.6723*10^-5; % gravitational constant
% prepare the sources
iprec= 5; %5; 
nsource = K+lsrf; 

source = [v0elm(:,1:3)' v5srf(:,1:3)']; 


ifcharge = 1; 
charge = [v0elm(:,4)' v5srf(:,6)']; 

% prepar the targets
ntarget = size(pxyz,1); 
target = pxyz';


% other parameters
ifdipole = 0; dipstr = zeros(1,nsource); 
dipvec = rand(3,nsource); ifpot = 0; iffld = 0;
ifpottarg = 1; iffldtarg = 1;

% create an initial structure
U.pottarg = zeros(1,ntarget); 
U.fldtarg = zeros(3,ntarget); 
% memory control
if mod(nsource,Gpsiz) == 0
   ns = nsource/Gpsiz; 
else
   ns = (nsource - mod(nsource,Gpsiz))/Gpsiz + 1;
end   
   
if mod(ntarget,Gpsiz) == 0
   nt = ntarget/Gpsiz; 
else
   nt = (ntarget - mod(ntarget,Gpsiz))/Gpsiz + 1;
end  


nss = 1;  
for i = 1:ns
    if mod(nsource,Gpsiz)~=0 && i == ns
       ns0 = mod(nsource,Gpsiz);
    else
       ns0 = Gpsiz;
    end
    
    nse = nss + ns0 - 1; 
    nts = 1; 
    for j = 1:nt
        if mod(ntarget,Gpsiz)~=0 && j == nt
           nt0 = mod(ntarget,Gpsiz);
        else
           nt0 = Gpsiz;
        end
        nte = nts + nt0 - 1; 
        
        [U0]=lfmm3dpart(iprec,ns0,source(:,nss:nse),...
             ifcharge,charge(:,nss:nse),...
             ifdipole,[],dipvec(:,nss:nse),ifpot,iffld,...
             nt0,target(:,nts:nte),ifpottarg,iffldtarg);

        U.pottarg(:,nts:nte) = U.pottarg(:,nts:nte)+U0.pottarg;
        U.fldtarg(:,nts:nte) = U.fldtarg(:,nts:nte)+U0.fldtarg;
        clear U0;
        nts = nte + 1;
    end
    nss = nse+1; 
    
end



