function trisurf2poly(fname,p_in,t_in,regs,boxrange)

Nnode=size(p_in,1);
Nface=size(t_in,1);
boxed=0;
Nregs=0;

if nargin>=4
    Nregs=size(regs,1);
end

if nargin>=5 && length(boxrange)>=6
%   boxrange=[xmin xmax ymin ymax z1 z2 ... zN]
    Boxnode=[boxrange(1) boxrange(3)
             boxrange(2) boxrange(3)
             boxrange(2) boxrange(4)
             boxrange(1) boxrange(4)];
    Boxface=[1,2,3,4;
             1,2,6,5;
             2,3,7,6;
             3,4,8,7;
             4,1,5,8];
    Nlayer=length(boxrange)-4;
    boxnode=zeros(4*Nlayer,3);
    boxface=zeros(5*(Nlayer-1)+1,4);
    for i=1:Nlayer-1
        boxnode((i-1)*4+(1:4),1:2)=Boxnode;
        boxnode((i-1)*4+(1:4),3)=boxrange(i+4);
        boxface((i-1)*5+(1:5),:)=Boxface+4*(i-1);
    end
    boxnode((Nlayer-1)*4+(1:4),1:2)=Boxnode;
    boxnode((Nlayer-1)*4+(1:4),3)=boxrange(end);
    boxface(end,:)=Boxface(1,:)+4*(Nlayer-1);
    boxface=boxface+size(p_in,1);
    p_in=[p_in;boxnode];
    Nnode=Nnode+size(boxnode,1);
    Nface=Nface+size(boxface,1);
    boxed=1;
end

fid=fopen([fname,'.poly'],'w');
fprintf(fid,'%d\t%d\t%d\t%d\n',Nnode,3,0,0);
fprintf(fid,'%d\t%f\t%f\t%f\n',[(1:Nnode)',p_in]');
fprintf(fid,'%d\t%d\n',Nface,1);
for i=1:size(t_in,1)
    fprintf(fid,'%d\n',1);
    fprintf(fid,'%d\t%d\t%d\t%d\n',[3,t_in(i,:)]);
end
if boxed==1
    for i=1:size(boxface,1)
        fprintf(fid,'%d\n',1);
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',[4,boxface(i,:)]);
    end
end
fprintf(fid,'%d\n',0);
fprintf(fid,'%d\n',Nregs);
if Nregs>0
    for i=1:Nregs
        fprintf(fid,'%d\t',i);
        fprintf(fid,'%f\t',regs(i,:));
        fprintf(fid,'%d\n',i);
    end
end

fclose(fid);