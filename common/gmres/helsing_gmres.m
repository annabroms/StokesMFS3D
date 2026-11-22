function [x,it,resvec,trueres] = helsing_gmres( ...
    f,b,n,m,tol,debug,verbose,progressDots,reorthogonalize)
% *** GMRES with low-threshold stagnation control taken from Johan Helsing's RCIP demo***
% Optional progressDots prints one dot per iteration, wrapping every 50 dots.
% Optional reorthogonalize adds a second Gram-Schmidt pass and checks
% the recomputed residual before accepting convergence. This is useful
% for strongly coupled image systems. Once the recurrence reaches tol,
% five nonimproving true-residual checks return the best checked iterate;
% trueres then reports the unmet tolerance. Legacy calls retain their behavior.
if nargin < 6
    debug = 0;
    verbose = 0;
elseif nargin < 7 
    verbose = 1;
end
if nargin < 8
    progressDots = false;
end
if nargin < 9
    reorthogonalize = false;
end
V=zeros(n,m+1);
H=zeros(m);
cs=zeros(m,1);
sn=zeros(m,1);
resvec=zeros(1,m);
bestx=zeros(n,1);
bestres=inf;
stagnation=0;
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
if reorthogonalize
  % Recover orthogonality lost in a long, strongly coupled Krylov basis.
  correction=V(:,1:it)'*w;
  w=w-V(:,1:it)*correction;
  H(1:it,it)=H(1:it,it)+correction;
end
H(it,it)=H(it,it)+1;
wnrm2=norm(w);
if wnrm2>0
  V(:,it1)=w/wnrm2;
end
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
    if progressDots
        fprintf('.');
        if mod(it,50) == 0
            fprintf('\n');
        end
    end
    if debug
        disp(it);
        disp(myerr);
    end

if ((myerr<=tol)||(it==m)||(wnrm2==0))
  y=triu(H(1:it,1:it))\s(1:it);             
  x=fliplr(V(:,1:it))*flipud(y);
  %trueres=norm(A*x-b)/bnrm2;
  trueres=norm(f(x)-b)/bnrm2;
  if reorthogonalize
    if trueres<bestres
      bestx=x;
      bestres=trueres;
      stagnation=0;
    else
      stagnation=stagnation+1;
    end
  end
  done=~reorthogonalize || trueres<=tol || it==m || wnrm2==0 || stagnation>=5;
  if done && reorthogonalize
    x=bestx;
    trueres=bestres;
  end
  if done && progressDots && mod(it,50) ~= 0
    fprintf('\n');
  end
  if verbose && done
    disp(['predicted residual = ' num2str(myerr)])
    disp(['true residual      = ',num2str(trueres)])
    disp(['iterations taken:' ,num2str(it)])
  end
  if done
    resvec=resvec(1:it);
    break
  end
end
end
end


