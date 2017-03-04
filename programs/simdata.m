function [x, u] = simdata(change, nrep, nobs, seed)
%Simulate data for abrupt or slowly varying change

%%%Simulate covariate values(u) and time series(x)
rng(seed);

u=linspace(0,1,nrep);

if strcmp(change,'abrupt')
 
    sd_true(1)=1;
    sd_true(2)=1;
    e1=normrnd(0,sd_true(1),nrep,floor(nobs/2));
    e2=normrnd(0,sd_true(2),nrep,floor(nobs/2));

    %u<=0.5
    phi_true=zeros(2,1);
    phi_true(1)=-0.5;
    phi_true(2)=0.5;
    x1=zeros(nrep,floor(nobs/2));
    x2=zeros(nrep,floor(nobs/2));   
    for i=1:(nrep/2)
        x1(i,:)=filter(1,[1 -phi_true(1)],e1(i,:));
        x2(i,:)=filter(1,[1 -phi_true(2)],e2(i,:));
    end

    %u>0.5
    phi_true=zeros(2,1);
    phi_true(1)=-0.9;
    phi_true(2)=0.9;

    for i=((nrep/2)+1):nrep
        x1(i,:)=filter(1,[1 -phi_true(1)],e1(i,:));
        x2(i,:)=filter(1,[1 -phi_true(2)],e2(i,:));
    end
    x=[x1 x2]';
  
elseif strcmp(change,'slow')

    sd_true=1;
    e=normrnd(0,sd_true,nrep,nobs);
    x=zeros(nrep,nobs);
    for i=1:nrep
        if u(i) <= 0.5
            %slowly varying -0.5 to 0.5
            phi_true=-0.5+(1:nobs)/nobs;
            for j=1:nobs
               if j==1
                   x(i,j)=e(i,j);
               else
                   x(i,j)=phi_true(j)*x(i,j-1)+e(i,j);
               end
            end   
        else
            %slowly varying from -0.9 to 0.9
            phi_true=-0.9+((1:nobs)/nobs)*1.8;
            for j=1:nobs
               if j==1
                   x(i,j)=e(i,j);
               else
                   x(i,j)=phi_true(j)*x(i,j-1)+e(i,j);
               end
            end
        end
    end
    x=x';
        
else
    error('Change input should be either abrupt or slow.');
end

for i=1:nrep
    xmat=[ones(length(x),1) linspace(1,length(x),length(x))'];
    linfit=inv(xmat'*xmat)*xmat'*x(:,i);
    x(:,i)=x(:,i)-xmat*linfit;
end

end

