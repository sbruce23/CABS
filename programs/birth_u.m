function[met_rat,nseg_prop,ui_prop,tau_prop,beta_prop]=birth_u(x,u,...
    nexp_tcurr_temp,nexp_ucurr_temp,nexp_prop,...
    tau_curr_temp,xi_curr_temp,ui_curr_temp,...
    nseg_curr_temp,nseg_cov_curr_temp,beta_curr_temp,log_move_curr,log_move_prop,opt)
%Function to propose birth move in covariate dimension and acceptance probability

nbeta=opt.nbasis+1;
nu_unique=size(unique(u),2);
beta_prop=zeros(nbeta,nexp_tcurr_temp,nexp_prop);
tau_prop=ones(nexp_tcurr_temp,nexp_prop);
nseg_prop=zeros(nexp_prop,1);
ui_prop=zeros(nexp_prop,1);

%Drawing  segment to split
%Available segments (taking duplicate covariate values into consideration)
u_sort=unique(sort(u));
kk=zeros(nexp_ucurr_temp,1);
for i=1:nexp_ucurr_temp
    if (i==1) && (sum(u_sort<=ui_curr_temp(i))>=2*opt.umin)
        kk(i)=i;
    elseif (i>1) && (sum(u_sort<=ui_curr_temp(i)&u_sort>ui_curr_temp(i-1))>=2*opt.umin)
        kk(i)=i;
    end
end
kk=kk(kk~=0)';

nposs_seg=length(kk);
seg_cut=kk(unidrnd(nposs_seg));%Drawing segment to split

%Drawing new cutpoint
if seg_cut==1
    nposs_cut=sum(u_sort<=ui_curr_temp(seg_cut))-2*opt.umin+1;
else
    nposs_cut=sum(u_sort<=ui_curr_temp(seg_cut)&u_sort>ui_curr_temp(seg_cut-1))-...
        2*opt.umin+1;
end
   
for jj=1:nexp_ucurr_temp
    if jj<seg_cut
        ui_prop(jj)=ui_curr_temp(jj);
		tau_prop(:,jj)=tau_curr_temp(:,jj);
		nseg_prop(jj)=nseg_cov_curr_temp(jj);
		beta_prop(:,:,jj)=beta_curr_temp(:,:,jj);
    elseif jj==seg_cut
        index=unidrnd(nposs_cut);
        if (seg_cut==1)
            ui_prop(seg_cut)=u_sort(index+opt.umin-1);
        else
            ui_prop(seg_cut)=u_sort(sum(u_sort<=ui_curr_temp(seg_cut-1))-1+opt.umin+index);
        end
        ui_prop(seg_cut+1)=ui_curr_temp(jj);
        zz=rand(nexp_tcurr_temp,1);%Drawing new tausq
        tau_prop(:,seg_cut)=tau_curr_temp(:,seg_cut).*(zz./(1-zz));
        tau_prop(:,seg_cut+1)=tau_curr_temp(:,seg_cut).*((1-zz)./zz);
       
        if (seg_cut==1)
            nseg_prop(seg_cut)=sum(u_sort<=ui_prop(seg_cut));
        else
            nseg_prop(seg_cut)=sum(u_sort<=ui_prop(seg_cut)&u_sort>ui_prop(seg_cut-1));
        end
        nseg_prop(seg_cut+1)=sum(u_sort<=ui_prop(seg_cut+1)&u_sort>ui_prop(seg_cut));
        
        for k=jj:jj+1
            for i=1:nexp_tcurr_temp
                [beta_mean, beta_var, ~, ~]=postbeta(i,k,nseg_curr_temp(i),x,u,xi_curr_temp,ui_prop,tau_prop(i,k),opt);
                beta_prop(:,i,k)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'))';%Drawing a new value of beta
            end
        end
    else
        ui_prop(jj+1)=ui_curr_temp(jj);
		tau_prop(:,jj+1)=tau_curr_temp(:,jj);
		nseg_prop(jj+1)=nseg_cov_curr_temp(jj);
		beta_prop(:,:,jj+1)=beta_curr_temp(:,:,jj);
    end
end
%Calculating Jacobian(sum of log Jacobian for each covariate seg)
log_jacobian=sum(log(2*tau_curr_temp(:,seg_cut)./(zz.*(1-zz))));
%=======================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
%=======================================
loglike_prop=0;
log_beta_prop=0;
log_beta_prior_prop=0;
log_tau_prior_prop=0;
for jj=seg_cut:seg_cut+1
    for i=1:nexp_tcurr_temp
        [beta_mean, beta_var, nu_mat, y]=postbeta(i,jj,nseg_curr_temp(i),x,u,xi_curr_temp,ui_prop,tau_prop(i,jj),opt);
        log_beta_prop=log_beta_prop+...
            log(mvnpdf(beta_prop(:,i,jj),beta_mean,0.5*(beta_var+beta_var')));
        log_beta_prior_prop=log_beta_prior_prop+...
            log(mvnpdf(beta_prop(:,i,jj),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_prop(i,jj)*ones(opt.nbasis,1)])));%Prior Density of beta
        log_tau_prior_prop=log_tau_prior_prop-log(tau_prop(i,jj));%Prior Density of Tausq
        fhat=nu_mat*beta_prop(:,i,jj);
        [log_prop_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));%Loglikelihood  at proposed values
        loglike_prop=loglike_prop+log_prop_spec_dens;
    end
end

log_seg_prop=-log(nposs_seg);%Proposal for segment choice
log_cut_prop=-log(nposs_cut);%Proposal for cut point choice
%Evaluating prior density for cut points at proposed values
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nu_unique-(nexp_prop-k+1)*opt.umin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nu_unique-...
            sum(u_sort<=ui_prop(k-1))-(nexp_prop-k+1)*opt.umin+1);
	end
end


%Calculating log proposal density at proposed values
log_proposal_prop=log_beta_prop+log_seg_prop+log_move_prop+log_cut_prop;
%Calculating log prior density at proposed values
log_prior_prop= log_beta_prior_prop+log_tau_prior_prop+log_prior_cut_prop;
%Calculating target density at proposed values
log_target_prop=loglike_prop+log_prior_prop;		
%*************************************************************
%CURRENT VALUES		
%*************************************************************
%=======================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Current values
%=======================================
loglike_curr=0;
log_beta_curr=0;
log_beta_prior_curr=0;
log_tau_prior_curr=0;
for i=1:nexp_tcurr_temp
    %Beta proposal and prior
    [beta_mean, beta_var, nu_mat, y]=postbeta(i,seg_cut,nseg_curr_temp(i),x,u,xi_curr_temp,ui_curr_temp,tau_curr_temp(i,seg_cut),opt);
    log_beta_curr=log_beta_curr+log(mvnpdf(beta_curr_temp(:,i,seg_cut),beta_mean,0.5*(beta_var+beta_var')));
    log_beta_prior_curr=log_beta_prior_curr+log(mvnpdf(beta_curr_temp(:,i,seg_cut),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_curr_temp(i,seg_cut)*ones(opt.nbasis,1)])));
    log_tau_prior_curr=log_tau_prior_curr-log(tau_curr_temp(i,seg_cut));
    %Log likelihood  at current values
    fhat=nu_mat*beta_curr_temp(:,i,seg_cut);
    [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));
    loglike_curr=loglike_curr+log_curr_spec_dens;
end
%Calculating log proposal density at current values
log_proposal_curr=log_beta_curr+log_move_curr;
%Evaluating  prior density for cut points at current values
log_prior_cut_curr=0;
for k=1:nexp_ucurr_temp-1
	if k==1
		log_prior_cut_curr=-log(nu_unique-(nexp_ucurr_temp-k+1)*opt.umin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nu_unique-...
            sum(u_sort<=ui_curr_temp(k-1))-(nexp_ucurr_temp-k+1)*opt.umin+1);
	end
end

%Calculating priors at current values
log_prior_curr=log_beta_prior_curr+log_tau_prior_curr+log_prior_cut_curr;
%Evaluating target densities at current values
log_target_curr=loglike_curr+log_prior_curr;

met_rat=min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop+log_jacobian));