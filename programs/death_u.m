function[met_rat,nseg_prop,ui_prop,tau_prop,beta_prop]=death_u(x,u,...
    nexp_tcurr_temp,nexp_ucurr_temp,nexp_prop,...
    tau_curr_temp,xi_curr_temp,ui_curr_temp,...
    nseg_curr_temp,nseg_cov_curr_temp,beta_curr_temp,log_move_curr,log_move_prop,opt)
%Function to propose death move in covariate dimension and acceptance probability

nbeta=opt.nbasis+1;
nu_unique=size(unique(u),2);

beta_prop=zeros(nbeta,nexp_tcurr_temp,nexp_prop);
tau_prop=ones(nexp_tcurr_temp,nexp_prop);
nseg_prop=zeros(nexp_prop,1);
ui_prop=zeros(nexp_prop,1);

%Drawing cut point to delete
cut_del=unidrnd(nexp_ucurr_temp-1);
j=0;
for k=1:nexp_prop
	j=j+1;
	if k==cut_del
		%*************************************************************
		%PROPOSED VALUES		
		%*************************************************************
		ui_prop(k)=ui_curr_temp(j+1);
		tau_prop(:,k)=sqrt(tau_curr_temp(:,j).*tau_curr_temp(:,j+1));% Combining 2 taus into 1		
        nseg_prop(k)=nseg_cov_curr_temp(j)+nseg_cov_curr_temp(j+1);% Combining two segments into 1
        %==================================================================
        %Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
        %==================================================================		
        loglike_prop=0;
        log_beta_prop=0;
        log_beta_prior_prop=0;
        log_tau_prior_prop=0;
        %Computing mean and variances for beta proposals
        for i=1:nexp_tcurr_temp
            [beta_mean, beta_var, nu_mat, y]=postbeta(i,k,nseg_curr_temp(i),x,u,xi_curr_temp,ui_prop,tau_prop(i,k),opt);
            beta_prop(:,i,k)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'));%Drawing a new value of beta
            %Loglikelihood  at proposed values
            fhat=nu_mat*beta_prop(:,i,k);
            [log_prop_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));
            loglike_prop=loglike_prop+log_prop_spec_dens;
            %==================================================================
            %Evaluating the Prior Densities at the Proposed values for tau and
            %beta==============================================================
            % Beta
            log_beta_prior_prop=log_beta_prior_prop+log(mvnpdf(beta_prop(:,i,k),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_prop(i,k)*ones(opt.nbasis,1)])));
            %Tau
            log_tau_prior_prop=log_tau_prior_prop-log(tau_prop(i,k));
        
            %==================================================================
            %Evaluating the Proposal Densities at the Proposed values of beta,
            %the cut points
            %==================================================================
            %Beta
            log_beta_prop=log_beta_prop+log(mvnpdf(beta_prop(:,i,k),beta_mean,0.5*(beta_var+beta_var')));
        end   	
		%Segment
		log_seg_prop=-log(nexp_ucurr_temp-1);
        %Calculating Jacobian(sum of log Jacobian for each covariate seg)
		log_jacobian=-sum(log(2*(sqrt(tau_curr_temp(:,j))+sqrt(tau_curr_temp(:,j+1))).^2));
		%Calculating log proposal density at proposed values
		log_proposal_prop=log_beta_prop+log_seg_prop+log_move_prop;
        
		%*************************************************************
		%CURRENT VALUES		
		%*************************************************************
		%=======================================
		%Evaluating the Likelihood, Proposal and Prior Densities at the Current values
		%=======================================
		%Beta proposal and prior
		loglike_curr=0;
        log_beta_curr=0;
		log_beta_prior_curr=0;
        log_tau_prior_curr=0;
        for i=1:nexp_tcurr_temp
            for jj=j:j+1               
                [beta_mean, beta_var, nu_mat, y]=postbeta(i,jj,nseg_curr_temp(i),x,u,xi_curr_temp,ui_curr_temp,tau_curr_temp(i,jj),opt);
                log_beta_curr=log_beta_curr+...
                    log(mvnpdf(beta_curr_temp(:,i,jj),beta_mean,0.5*(beta_var+beta_var')));
                log_beta_prior_curr=log_beta_prior_curr+...
                    log(mvnpdf(beta_curr_temp(:,i,jj),zeros(nbeta,1),diag([opt.sigmasqalpha; tau_curr_temp(i,jj)*ones(opt.nbasis,1)])));
                log_tau_prior_curr=log_tau_prior_curr-log(tau_curr_temp(i,jj));	
                %Loglikelihood  at current values
                fhat=nu_mat*beta_curr_temp(:,i,jj);
                [log_curr_spec_dens]=whittle_like(y,fhat,nseg_curr_temp(i));
                loglike_curr=loglike_curr+log_curr_spec_dens;
            end
        end
        %Calculating log proposal density at current values
        log_proposal_curr=log_move_curr+log_beta_curr;
        %Calculating priors at current values
        log_prior_curr=log_beta_prior_curr+log_tau_prior_curr;
        j=j+1;
	else
		ui_prop(k)=ui_curr_temp(j);
		tau_prop(:,k)=tau_curr_temp(:,j);
		nseg_prop(k)=nseg_cov_curr_temp(j);
		beta_prop(:,:,k)=beta_curr_temp(:,:,j);
	end
end


%Evaluating Target density at proposed values(taking duplicates into
%consideration)
u_sort=unique(sort(u));
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nu_unique-(nexp_prop-k+1)*opt.umin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nu_unique-...
            sum(u_sort<=ui_prop(k-1))-(nexp_prop-k+1)*opt.umin+1);
	end
end
log_target_prop=loglike_prop+log_tau_prior_prop+log_beta_prior_prop+log_prior_cut_prop;
%Evaluating Target density at current values(taking duplicates into
%consideration)
log_prior_cut_curr=0;
for k=1:nexp_ucurr_temp-1
	if k==1
		log_prior_cut_curr=-log(nu_unique-(nexp_ucurr_temp-k+1)*opt.umin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nu_unique-...
            sum(u_sort<=ui_curr_temp(k-1))-(nexp_ucurr_temp-k+1)*opt.umin+1);
	end
end
log_target_curr=loglike_curr+log_prior_curr+log_prior_cut_curr;

met_rat=min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop+log_jacobian));