function [x,it,resvec,trueres] = helsing_gmres_mv(f,b,n,m,tol,debug,grid)
% *** GMRES with low-threshold stagnation control ***
if nargin < 6
    debug = 0;
end
V=zeros(n,m+1);
H=zeros(m);
cs=zeros(m,1);
sn=zeros(m,1);
bnrm2=norm(b);
V(:,1)=b/bnrm2;
s=bnrm2*eye(m+1,1);
for it=1:m  
it1=it+1;                                   
%w=A*V(:,it);
w=f(V(:,it))-V(:,it);
for k=1:it
  H(k,it)=V(:,k)'*w;
  w=w-H(k,it)*V(:,k);
end
H(it,it)=H(it,it)+1;
wnrm2=norm(w);
V(:,it1)=w/wnrm2;
for k=1:it-1                                
  temp     = cs(k)*H(k,it)+sn(k)*H(k+1,it);
  H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
  H(k,it)  = temp;
end
[cs(it),sn(it)]=rotmat(H(it,it),wnrm2);     
H(it,it)= cs(it)*H(it,it)+sn(it)*wnrm2;
s(it1) =-sn(it)*s(it);                      
s(it)  = cs(it)*s(it);                         
myerr=abs(s(it1))/bnrm2;
resvec(it) = myerr;
    if debug
        it
        myerr
        
        y=triu(H(1:it,1:it))\s(1:it);             
        x=fliplr(V(:,1:it))*flipud(y);

       %  res = f(x);
       %  err = res-b;
       % 
       %  err = err(1:end/2)+1i*err(end/2+1:end); 
       % % err_max = [err_max; max(abs(err))];
       % 
       %  figure(58)
       %  clf; 
       %  plot(res)
       %  hold on
       %  plot(b,'+')
       %  %pause(1);
       % 
       %  figure(59)
       %  clf;
       %  scatter(real(grid),imag(grid),10,log10(abs(err)),'filled'); 
       %  colorbar
       %  caxis([-6,1])
       %  str = sprintf('Max err %1.3e',max(abs(err)));
       %  title(str)
        %pause(1);
        
%         figure(60)
%         clf; 
%         plot(x)



    end





if (myerr<=tol)||(it==m)                     
  disp(['predicted residual = ' num2str(myerr)])
  y=triu(H(1:it,1:it))\s(1:it);             
  x=fliplr(V(:,1:it))*flipud(y);
  %trueres=norm(A*x-b)/bnrm2;
  trueres=norm(f(x)-b)/bnrm2;
  disp(['true residual      = ',num2str(trueres)])
  disp(['iterations taken:' ,num2str(it)])
  break
end
end
end


