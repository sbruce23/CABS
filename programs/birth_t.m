function[met_rat,nseg_prop,xi_prop,tau_prop,beta_prop]=birth_t(x,u,...
    nexp_tcurr_temp,nexp_ucurr_temp,nexp_prop,...
    tau_curr_temp,xi_curr_temp,ui_curr_temp,...
    nseg_curr_temp,beta_curr_temp,log_move_curr,log_move_prop,opt)
%Function to propose birth move in time dimension and acceptance probability

nobs=size(x,1);
nbeta=opt.nbasis+1;
beta_prop=zeros(nbeta,nexp_prop,nexp_ucurr_temp);
tau_prop=ones(nexp_prop,nexp_ucurr_temp);
nseg_prop=zeros(nexp_prop,1);
xi_prop=zeros(nexp_prop,1);

%Drawing  segment to split
kk=find(nseg_curr_temp>=2*opt.tmin); %Number of segments available for splitting
nposs_seg=length(kk);
seg_cut=kk(unidrnd(nposs_seg));%Drawing segment to split
nposs_cut=nseg_curr_temp(seg_cut)-2*opt.tmin+1;%Drawing new cutpoint

for jj=1:nexp_tcurr_temp
    if jj<seg_cut
        xi_prop(jj)=xi_curr_temp(jj);
		tau_prop(jj,:)=tau_curr_temp(jj,:);
		nseg_prop(jj)=nseg_curr_temp(jj);
		beta_prop(:,jj,:)=beta_curr_temp(:,jj,:);
    elseif jj==seg_cut
        index=unidrnd(nposs_cut);
        if (seg_cut==1)
            xi_prop(seg_cut)=index+opt.tmin-1;
        else
            xi_prop(seg_cut)=xi_curr_temp(jj-1)-1+opt.tmin+index;
        end
        xi_prop(seg_cut+1)=xi_curr_temp(jj);
        zz=rand(1,nexp_ucurr_temp);%Drawing new tausq
        tau_prop(seg_cut,:)=tau_curr_temp(seg_cut,:).*(zz./(1-zz));
        tau_prop(seg_cut+1,:)=tau_curr_temp(seg_cut,:).*((1-zz)./zz);
        nseg_prop(seg_cut)=index+opt.tmin-1;
        nseg_prop(seg_cut+1)=nseg_curr_temp(jj)-nseg_prop(seg_cut);                  
        for k=jj:jj+1
            for i=1:nexp_ucurr_temp
                [beta_mean, beta_var, ~, ~]=postbeta(k,i,nseg_prop(k),x,u,xi_prop,ui_curr_temp,tau_prop(k,i),opt);
                beta_prop(:,k,i)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'))';%Drawing a new value of beta
            end
        end
    else
        xi_prop(jj+1)=xi_curr_temp(jj);
		tau_prop(jj+1,:)=tau_curr_temp(jj,:);
		nseg_prop(jj+1)=nseg_curr_temp(jj);
		beta_prop(:,jj+1,:)=beta_curr_temp(:,jj,:);
    end
end

%Calculating Jacobian(sum of log Jacobian for each covariate seg)
log_jacobian=sum(log(2*tau_curr_temp(seg_cut,:)./(zz.*(1-zz))));

%=======================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
%=======================================
loglike_prop=0;
log_beta_prop=0;
log_beta_prior_prop=0;
log_tau_prior_prop=0;
for jj=seg_cut:seg_cut+1
    for i=1:nexp_ucurr_temp
        [beta_mean, beta_var, nu_mat, y]=postbeta(jj,i,nseg_prop(jj),x,u,xi_prop,ui_curr_temp,tau_prop(jj,i),opt);
        log_beta_prop=log_beta_prop+...
            log(mvnpdf(beta_prop(:,jj,i),beta_mean,0.5*(beta_var+beta_var')));
        log_beta_prior_prop=log_beta_prior_prop+...
            log(mvnpdf(beta_prop(:,jj,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_prop(jj,i)*ones(opt.nbasis,1)])));%Prior Density of beta
        log_tau_prior_prop=log_tau_prior_prop-log(tau_prop(jj,i));%Prior density of tausq
        fhat=nu_mat*beta_prop(:,jj,i);
        [log_prop_spec_dens]=whittle_like(y,fhat,nseg_prop(jj));%Log likelihood  at proposed values
        loglike_prop=loglike_prop+log_prop_spec_dens;
    end
end
log_seg_prop=-log(nposs_seg);%Proposal for segment choice
log_cut_prop=-log(nposs_cut);%Proposal for cut point choice
%Evaluating prior density for cut points at proposed values
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*opt.tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*opt.tmin+1);
	end
end
%Calculating Log Proposal density at proposed values
log_proposal_prop=log_beta_prop+log_seg_prop+log_move_prop+log_cut_prop;
%Calculating Log Prior density at proposed values
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
for i=1:nexp_ucurr_temp
    %Beta proposal and prior
    [beta_mean, beta_var, nu_mat, y]=postbeta(seg_cut,i,nseg_curr_temp(seg_cut),x,u,xi_curr_temp,ui_curr_temp,tau_curr_temp(seg_cut,i),opt);
    log_beta_curr=log_beta_curr+log(mvnpdf(beta_curr_temp(:,seg_cut,i),beta_mean,0.5*(beta_var+beta_var')));
    log_beta_prior_curr=log_beta_prior_curr+log(mvnpdf(beta_curr_temp(:,seg_cut,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_curr_temp(seg_cut,i)*ones(opt.nbasis,1)])));
    log_tau_prior_curr=log_tau_prior_curr-log(tau_curr_temp(seg_cut,i));
    %Log likelihood  at current values
    fhat=nu_mat*beta_curr_temp(:,seg_cut,i);
    [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(seg_cut));
    loglike_curr=loglike_curr+log_curr_spec_dens;
end
%Calculating log proposal density at current values
log_proposal_curr=log_beta_curr+log_move_curr;
%Evaluating  prior density for cut points at current values
log_prior_cut_curr=0;
for k=1:nexp_tcurr_temp-1
	if k==1
		log_prior_cut_curr=-log(nobs-(nexp_tcurr_temp-k+1)*opt.tmin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_tcurr_temp-k+1)*opt.tmin+1);
	end
end
%Calculating priors at current values
log_prior_curr=log_beta_prior_curr+log_tau_prior_curr+log_prior_cut_curr;
%Evaluating target densities at current values
log_target_curr=loglike_curr+log_prior_curr;

met_rat=min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop+log_jacobian));