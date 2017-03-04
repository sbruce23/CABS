function[met_rat,nseg_prop,xi_prop,tau_prop,beta_prop]=death_t(x,u,...
    nexp_tcurr_temp,nexp_ucurr_temp,nexp_prop,...
    tau_curr_temp,xi_curr_temp,ui_curr_temp,...
    nseg_curr_temp,beta_curr_temp,log_move_curr,log_move_prop,opt)
%Function to propose death move in time dimension and acceptance probability

nobs=size(x,1);
nbeta=opt.nbasis+1;
beta_prop=zeros(nbeta,nexp_prop,nexp_ucurr_temp);
tau_prop=ones(nexp_prop,nexp_ucurr_temp);
nseg_prop=zeros(nexp_prop,1);
xi_prop=zeros(nexp_prop,1);

%Drawing cut point to delete
cut_del=unidrnd(nexp_tcurr_temp-1);
j=0;
for k=1:nexp_prop
	j=j+1;
	if k==cut_del
		%*************************************************************
		%PROPOSED VALUES		
		%*************************************************************
		xi_prop(k)=xi_curr_temp(j+1);
		tau_prop(k,:)=sqrt(tau_curr_temp(j,:).*tau_curr_temp(j+1,:));% Combining 2 taus into 1
		nseg_prop(k)=nseg_curr_temp(j)+nseg_curr_temp(j+1);% Combining two segments into 1
        %==================================================================
        %Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
        %==================================================================		
        loglike_prop=0;
        log_beta_prop=0;
        log_beta_prior_prop=0;
        log_tau_prior_prop=0;
        %Computing mean and variances for beta proposals
        for i=1:nexp_ucurr_temp
            [beta_mean, beta_var, nu_mat, y]=postbeta(k,i,nseg_prop(k),x,u,xi_prop,ui_curr_temp,tau_prop(k,i),opt);
            beta_prop(:,k,i)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'));%Drawing a new value of beta
            %Loglikelihood  at proposed values
            fhat=nu_mat*beta_prop(:,k,i);
            [log_prop_spec_dens]=whittle_like(y,fhat,nseg_prop(k));
            loglike_prop=loglike_prop+log_prop_spec_dens;
            %==================================================================
            %Evaluating the Prior Densities at the Proposed values for tau and beta
            %==================================================================
            %Beta
            log_beta_prior_prop=log_beta_prior_prop+log(mvnpdf(beta_prop(:,k,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_prop(k,i)*ones(opt.nbasis,1)])));
            %Tau
            log_tau_prior_prop=log_tau_prior_prop-log(tau_prop(k,i));
        
            %==================================================================
            %Evaluating the Proposal Densities at the Proposed values of beta,
            %the cut points
            %==================================================================
            %Beta
            log_beta_prop=log_beta_prop+log(mvnpdf(beta_prop(:,k,i),beta_mean,0.5*(beta_var+beta_var')));
        end   	
		%Segment
		log_seg_prop=-log(nexp_tcurr_temp-1);
        %Calculating Jacobian(sum of log Jacobian for each covariate seg)
		log_jacobian=-sum(log(2*(sqrt(tau_curr_temp(j,:))+sqrt(tau_curr_temp(j+1,:))).^2));
		%Calculating log proposal density at proposed values
		log_proposal_prop=log_beta_prop+log_seg_prop+log_move_prop;
        
		%*************************************************************
		%CURRENT VALUES		
		%*************************************************************
		%=======================================
		%Evaluating the Likelihood, Proposal and Prior Densities at the current values
		%=======================================
		%Beta proposal and prior
		loglike_curr=0;
        log_beta_curr=0;
		log_beta_prior_curr=0;
        log_tau_prior_curr=0;
        for i=1:nexp_ucurr_temp
            for jj=j:j+1               
                [beta_mean, beta_var, nu_mat, y]=postbeta(jj,i,nseg_curr_temp(jj),x,u,xi_curr_temp,ui_curr_temp,tau_curr_temp(jj,i),opt);
                log_beta_curr=log_beta_curr+...
                    log(mvnpdf(beta_curr_temp(:,jj,i),beta_mean,0.5*(beta_var+beta_var')));
                log_beta_prior_curr=log_beta_prior_curr+...
                    log(mvnpdf(beta_curr_temp(:,jj,i),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_curr_temp(jj,i)*ones(opt.nbasis,1)])));
                log_tau_prior_curr=log_tau_prior_curr-log(tau_curr_temp(jj,i));	
                %Log likelihood at current values
                fhat=nu_mat*beta_curr_temp(:,jj,i);
                [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(jj));
                loglike_curr=loglike_curr+log_curr_spec_dens;
            end
        end
        %Calculating log proposal density at current values
        log_proposal_curr=log_move_curr+log_beta_curr;
        %Calculating priors at current values
        log_prior_curr=log_beta_prior_curr+log_tau_prior_curr;
        j=j+1;
	else
		xi_prop(k)=xi_curr_temp(j);
		tau_prop(k,:)=tau_curr_temp(j,:);
		nseg_prop(k)=nseg_curr_temp(j);
		beta_prop(:,k,:)=beta_curr_temp(:,j,:);
	end
end
%Evaluating target density at proposed values
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*opt.tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*opt.tmin+1);
	end
end
log_target_prop=loglike_prop+log_tau_prior_prop+log_beta_prior_prop+log_prior_cut_prop;
%Evaluating target density at current values
log_prior_cut_curr=0;
for k=1:nexp_tcurr_temp-1
	if k==1
		log_prior_cut_curr=-log(nobs-(nexp_tcurr_temp-k+1)*opt.tmin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_tcurr_temp-k+1)*opt.tmin+1);
	end
end
log_target_curr=loglike_curr+log_prior_curr+log_prior_cut_curr;

met_rat=min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop+log_jacobian));
