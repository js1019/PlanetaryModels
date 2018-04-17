function [x,y,z,link] = construct(meshname,N)

[p,tet,~,~]=read_mesh3d(meshname);


nodal=load('porder2nodal.mat','porder');
r=nodal.porder{N}.r;
s=nodal.porder{N}.s;
t=nodal.porder{N}.t;
pNp=length(r);
Nele=size(tet,1);

x = 0.5*(-(1+r+s+t)*p(tet(:,1),1)'...
             +(1+r)*p(tet(:,2),1)'...
             +(1+s)*p(tet(:,3),1)'...
             +(1+t)*p(tet(:,4),1)');
y = 0.5*(-(1+r+s+t)*p(tet(:,1),2)'...
             +(1+r)*p(tet(:,2),2)'...
             +(1+s)*p(tet(:,3),2)'...
             +(1+t)*p(tet(:,4),2)');
z = 0.5*(-(1+r+s+t)*p(tet(:,1),3)'...
             +(1+r)*p(tet(:,2),3)'...
             +(1+s)*p(tet(:,3),3)'...
             +(1+t)*p(tet(:,4),3)');
%x=x-min(x(:));y=y-min(y(:));z=z-min(z(:));


tt=tet_subelement(N);
Nsele=size(tt,1);
link=zeros(Nele*Nsele,4);
offset=0;
for i=1:Nele
    link((i-1)*Nsele+1:i*Nsele,:)=tt+offset;
    offset=offset+pNp;
end

