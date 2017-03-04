function[beta_mean, beta_var, nu_mat, y]=postbeta(i,j,nseg_time_temp,x,u,xi_temp,ui_temp,tau_temp,opt)
%Function returns mean and variance for normal approximation to posterior distribution of beta

nfreq=floor(nseg_time_temp/2);
freq=(0:nfreq)'/(2*nfreq);

% create set of indices for replications in jth covariate segment
if j==1
    uu=find(u<=ui_temp(j));
else
    uu=find(u<=ui_temp(j) & u>ui_temp(j-1));
end

% create log periodograms for replications in ith time seg and jth cov seg
y=zeros(nseg_time_temp,size(uu,2));
yy=zeros(nfreq+1,size(uu,2));
jj=0;
for k=uu
    jj=jj+1;
    if i>1            
            y(:,jj)=log(abs(fft(x(xi_temp(i-1)+1:xi_temp(i),k))).^2/nseg_time_temp);
            yy(:,jj)=y(1:nfreq+1,jj);
    else
            y(:,jj)=log(abs(fft(x(1:xi_temp(i),k))).^2/nseg_time_temp);
            yy(:,jj)=y(1:nfreq+1,jj);
    end
end

% pass log periodograms into optimizer to obtain beta_mean and beta_var
% for normal approximation to beta posterior
nbeta=opt.nbasis+1;
[nu_mat]=lin_basis_func(freq,nbeta);
nn=nseg_time_temp;
ytemp=yy;
param=zeros(nbeta,1);
[beta_mean,~,~,~,~,beta_inv_var]=...
	fminunc(@whittle_derivs2,param,opt.options,nn,nu_mat,ytemp,...
            tau_temp,nbeta,opt.nbasis,opt.sigmasqalpha); 
beta_var=opt.var_inflate*eye(nbeta)/beta_inv_var;
