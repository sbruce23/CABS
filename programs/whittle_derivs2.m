function [f,g,h] =whittle_derivs2(param,n,nu_mat,ytemp,tau_temp,nbeta,nbasis,sigmasqalpha)
%Function computes function, gradient, and hessian for Whittle likelihood
   
   ydev=bsxfun(@minus,ytemp,nu_mat*param);
   
   %n is the segment length in the time domain
   n1=floor(n/2);  
   if (mod(n,2)==1) %odd n
       f=sum(sum(bsxfun(@plus,nu_mat(2:n1+1,:)*param,exp(ydev(2:n1+1,:)))))+...
           0.5*(sum(nu_mat(1,:)*param+exp(ydev(1,:))))+...
           0.5*(param(2:nbeta)'*param(2:nbeta)/tau_temp+param(1)^2/sigmasqalpha);
   else
       f=sum(sum(bsxfun(@plus,nu_mat(2:n1,:)*param,exp(ydev(2:n1,:)))))+...
           0.5*(sum(nu_mat(1,:)*param+exp(ydev(1,:))))+...
           0.5*(sum(nu_mat(n1+1,:)*param+exp(ydev(n1+1,:))))+...
           0.5*(param(2:nbeta)'*param(2:nbeta)/tau_temp+param(1)^2/sigmasqalpha);
   end
   
   g=zeros(nbeta,1);
   g(1)=param(1)/sigmasqalpha;
   g(2:nbeta)=param(2:nbeta)/tau_temp;
   
   h=zeros(nbeta,nbeta);
   h(1,1)=1/sigmasqalpha;
   h(2:nbeta,2:nbeta)=1/tau_temp*eye(nbasis);
   
   
   if(mod(n,2)==1)
      for i=1:size(ydev,2) 
          temp_mat=bsxfun(@times,nu_mat(2:end,:),1-exp(ydev(2:end,i)));
          g=g+sum(temp_mat,1)'+0.5*(1-exp(ydev(1,i)))*nu_mat(1,:)';
      end
      jj=1:nbeta;
      big_mat=repmat(nu_mat(2:end,:)',nbeta,1).*nu_mat(2:end,jj(ones(nbeta,1),:)).';
      for i=1:size(ydev,2) 
          coef_mat=repmat(exp(ydev(2:end,i))',nbeta^2,1);
          h=h+sum(reshape(big_mat.*coef_mat,[nbeta nbeta n1]),3)+...
              0.5*exp(ydev(1,i))*nu_mat(1,:)'*nu_mat(1,:);
      end
   else
      for i=1:size(ydev,2)
          temp_mat=bsxfun(@times,nu_mat(2:n1,:),1-exp(ydev(2:n1,i)));
          g=g+sum(temp_mat,1)'+0.5*(1-exp(ydev(1,i)))*nu_mat(1,:)'+...
            0.5*(1-exp(ydev(n1+1,i)))*nu_mat(n1+1,:)';
      end
      jj=1:nbeta;
      big_mat=repmat(nu_mat(2:n1,:)',nbeta,1).*nu_mat(2:n1,jj(ones(nbeta,1),:)).';
      for i=1:size(ydev,2) 
          coef_mat=repmat(exp(ydev(2:n1,i))',nbeta^2,1);
          h=h+sum(reshape(big_mat.*coef_mat,[nbeta nbeta n1-1]),3)+...
              0.5*exp(ydev(1,i))*nu_mat(1,:)'*nu_mat(1,:)+...
              0.5*exp(ydev(n1+1,i))*nu_mat(n1+1,:)'*nu_mat(n1+1,:);
      end
   end
