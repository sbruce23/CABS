function[epsilon,ui_prop,beta_prop,nseg_new,size_w]=...
    within_u(x,u,nexp_tcurr_temp,nexp_ucurr_temp,...
    xi_curr_temp,ui_curr_temp,beta_curr_temp,...
    nseg_curr_temp,nseg_cov_curr_temp,tau_temp,opt)      
%Function to propose within model move in covariate dimension and acceptance probability

nbeta=opt.nbasis+1;
nrep=size(x,2);
nu_unique=size(unique(u),2);

ui_prop=ui_curr_temp;
beta_prop=beta_curr_temp;
nseg_new=nseg_cov_curr_temp;

u_sort=unique(sort(u));

if nexp_ucurr_temp>1
	seg_temp=unidrnd(nexp_ucurr_temp-1);%Drawing Segment to cut
	r=rand;    
    cut_poss_curr=sum(u_sort<=ui_curr_temp(seg_temp)); %index of current cut point    
    if (seg_temp==1)
        nposs_prior=sum(u_sort<=ui_curr_temp(seg_temp+1))-...
            2*opt.umin+1;
    else
        nposs_prior=sum(u_sort>ui_curr_temp(seg_temp-1)&u_sort<=ui_curr_temp(seg_temp+1))-...
        2*opt.umin+1;
    end 
    if (seg_temp==1)
        nseg_low=sum(u_sort<=ui_curr_temp(seg_temp));
    else
        nseg_low=sum(u_sort>ui_curr_temp(seg_temp-1)&u_sort<=ui_curr_temp(seg_temp));
    end
    nseg_high=sum(u_sort>ui_curr_temp(seg_temp)&u_sort<=ui_curr_temp(seg_temp+1));

    if r<opt.prob_mm1_u
		%small move
        size_w=1; 
               
        if nseg_low==opt.umin && nseg_high==opt.umin
            nposs=1;% Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new=cut_poss_curr-1+new_index;
        elseif nseg_low==opt.umin
			nposs=2;% Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new=cut_poss_curr-1+new_index;
        elseif nseg_high==opt.umin
			nposs=2;% Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new=cut_poss_curr+1-new_index;
        else
			nposs=3;% Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint 
			cut_poss_new=cut_poss_curr-2+new_index;
        end
    else		
        %big move
        size_w=2;  
        new_index=unidrnd(nposs_prior);       
		%cut_poss_new=sum(nseg_cov_curr_temp(1:seg_temp-1))-1+opt.umin+new_index;
        if (seg_temp==1)
            cut_poss_new=0-1+opt.umin+new_index;
        else
            cut_poss_new=sum(u_sort<=ui_curr_temp(seg_temp-1))-1+opt.umin+new_index;
        end
    end
	ui_prop(seg_temp)=u_sort(cut_poss_new);
    
    if(seg_temp>1)
		nseg_new(seg_temp)=sum(u<=ui_prop(seg_temp)&u>ui_prop(seg_temp-1));%Number of observations in lower part of new cutpoint
	else
		nseg_new(seg_temp)=sum(u<=ui_prop(seg_temp));
    end
	nseg_new(seg_temp+1)=sum(u<=ui_prop(seg_temp+1)&u>ui_prop(seg_temp));%Number of observations in upper part of new cutpoint
    
    %Evaluating the Proposal density for the cut-point at the current and proposed values
    if (seg_temp==1)
        nseg_new_low=sum(u_sort<=ui_prop(seg_temp));
    else
        nseg_new_low=sum(u_sort>ui_prop(seg_temp-1)&u_sort<=ui_prop(seg_temp));
    end
    nseg_new_high=sum(u_sort>ui_prop(seg_temp)&u_sort<=ui_prop(seg_temp+1));
        
    if(abs(cut_poss_new-cut_poss_curr)>1)
        log_prop_cut_prop=log(1-opt.prob_mm1_u)-log(nposs_prior);
        log_prop_cut_curr=log(1-opt.prob_mm1_u)-log(nposs_prior);
    elseif nseg_low==opt.umin && nseg_high==opt.umin
        log_prop_cut_prop=0;
        log_prop_cut_curr=0;
    else
        if nseg_low==opt.umin || nseg_high==opt.umin
           log_prop_cut_prop=log(1-opt.prob_mm1_u)-log(nposs_prior)+log(1/2)+log(opt.prob_mm1_u);
        else
           log_prop_cut_prop=log(1-opt.prob_mm1_u)-log(nposs_prior)+log(1/3)+log(opt.prob_mm1_u); 
        end
        if(nseg_new_low==opt.umin || nseg_new_high==opt.umin)
           log_prop_cut_curr=log(1-opt.prob_mm1_u)-log(nposs_prior)+log(1/2)+log(opt.prob_mm1_u);
        else
           log_prop_cut_curr=log(1-opt.prob_mm1_u)-log(nposs_prior)+log(1/3)+log(opt.prob_mm1_u); 
        end
    end
        
    %Evaluating the Loglikelihood, Priors and Proposals at the current values
	loglike_curr=0;
	log_beta_curr_temp=0;
	log_prior_curr=0;
    for i=seg_temp:seg_temp+1
        for j=1:nexp_tcurr_temp
            [beta_mean, beta_var, nu_mat, y]=postbeta(j,i,nseg_curr_temp(j),x,u,xi_curr_temp,ui_curr_temp,tau_temp(j,i),opt);
            %Compute log proposal density of beta at current  values
            log_beta_curr_temp=log_beta_curr_temp+log(mvnpdf(beta_curr_temp(:,j,i),beta_mean,0.5*(beta_var+beta_var')));
            fhat=nu_mat*beta_curr_temp(:,j,i);
            %Compute Loglike at current values
            [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(j));
            loglike_curr=loglike_curr+log_curr_spec_dens;
            %Compute priors at current values
            log_prior_curr=log_prior_curr+...
            log(mvnpdf(beta_curr_temp(:,j,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_temp(j,i)*ones(opt.nbasis,1)])));
        end
    end
%Evaluating the Loglikelihood, Priors and Proposals at the proposed values
%Likelihood
	loglike_prop=0;
	log_beta_prop=0;
	log_prior_prop=0;
    for i=seg_temp:seg_temp+1
        for j=1:nexp_tcurr_temp
            [beta_mean, beta_var, nu_mat, y]=postbeta(j,i,nseg_curr_temp(j),x,u,xi_curr_temp,ui_prop,tau_temp(j,i),opt);           
            beta_prop(:,j,i)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'))';	
            %Compute log proposal density of beta at proposed  values
            log_beta_prop=log_beta_prop+log(mvnpdf(beta_prop(:,j,i),beta_mean,0.5*(beta_var+beta_var')));
            fhat=nu_mat*beta_prop(:,j,i);
            %Compute Loglike at proposed values
            [log_prop_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(j));
            loglike_prop=loglike_prop+log_prop_spec_dens;
            %Compute priors at proposed values
            log_prior_prop=log_prior_prop+...
            log(mvnpdf(beta_prop(:,j,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_temp(j,i)*ones(opt.nbasis,1)])));
        end
    end
    %Proposal for beta
    log_proposal_curr=log_beta_curr_temp+log_prop_cut_curr;
    log_proposal_prop=log_beta_prop+log_prop_cut_prop;
    log_prior_cut_prop=0;
    log_prior_cut_curr=0;
  
    for k=1:nexp_ucurr_temp-1
        if k==1
			log_prior_cut_prop=-log(nu_unique-(nexp_ucurr_temp-k+1)*opt.umin+1);
			log_prior_cut_curr=-log(nu_unique-(nexp_ucurr_temp-k+1)*opt.umin+1);
		else
			log_prior_cut_prop=log_prior_cut_prop-...
                log(nu_unique-sum(u_sort<=ui_prop(k-1))-(nexp_ucurr_temp-k+1)*opt.umin+1);
			log_prior_cut_curr=log_prior_cut_curr-...
                log(nu_unique-sum(u_sort<=ui_curr_temp(k-1))-(nexp_ucurr_temp-k+1)*opt.umin+1);
        end
    end
 log_target_prop=loglike_prop+log_prior_prop+log_prior_cut_prop;
 log_target_curr=loglike_curr+log_prior_curr+log_prior_cut_curr;
else
	%no move since only one segment
    size_w=0;    
    nseg_new=nrep;
	log_beta_prop=0;
    log_beta_curr_temp=0;
    loglike_prop=0;
    loglike_curr=0;
    log_prior_prop=0;
    log_prior_curr=0;
    for i=1:nexp_tcurr_temp
        [beta_mean, beta_var, nu_mat, y]=postbeta(i,1,nseg_curr_temp(i),x,u,xi_curr_temp,ui_prop,tau_temp(i,1),opt);           
        beta_prop(:,i,1)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'))';
        %Compute log proposal density of beta at proposed  values
        log_beta_prop=log_beta_prop+log(mvnpdf(beta_prop(:,i,1),beta_mean,0.5*(beta_var+beta_var')));
        %Compute log proposal density of beta at current  values
        log_beta_curr_temp=log_beta_curr_temp+log(mvnpdf(beta_curr_temp(:,i,1),beta_mean,0.5*(beta_var+beta_var')));
        %Compute Loglike at proposed values
        fhat=nu_mat*beta_prop(:,i,1);
        [log_prop_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));
        loglike_prop=loglike_prop+log_prop_spec_dens;
        %Compute Loglike at proposed values
        fhat=nu_mat*beta_curr_temp(:,i,1);
        [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));
        loglike_curr=loglike_curr+log_curr_spec_dens;
        %Compute Priors at proposed values
        log_prior_prop=log_prior_prop+log(mvnpdf(beta_prop(:,i,1),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_temp(i,1)*ones(opt.nbasis,1)])));
        %Compute Priors at current values
        log_prior_curr=log_prior_curr+log(mvnpdf(beta_curr_temp(:,i,1),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_temp(i,1)*ones(opt.nbasis,1)])));
    end
        
    log_proposal_curr=log_beta_curr_temp;
	log_proposal_prop=log_beta_prop;
	log_target_prop=loglike_prop+log_prior_prop;
	log_target_curr=loglike_curr+log_prior_curr;
end
epsilon=min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop));