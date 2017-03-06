function [out] = CABS_sampler(x, u, opt)
%Main function to implement the CABS sampler [Bruce, Hall, Buysse, and Krafty (2017)]
%   Function takes in the set of replicated times series, covariates, and MCMC options
%   and outputs the results of the MCMC sampler in a MATLAB structure.

%%%set seed
rng(opt.seed);

%%%set up initial values
nobs=size(x,1);
nrep=size(x,2);
nu_unique=size(unique(u),2);
nbeta=opt.nbasis+1;
u_sort=sort(unique(u));

%number of segments
if opt.Rev_Jump_t==1
    nexp_tcurr=unidrnd(opt.nexp_tmax);%random draw if RJ allowed
else
    nexp_tcurr=opt.nexp_tmax;%max if RJ not allowed
end

if opt.Rev_Jump_u==1
    nexp_ucurr=unidrnd(opt.nexp_umax);
else
    nexp_ucurr=opt.nexp_umax;
end

%partitions
xi_curr=zeros(nexp_tcurr,1);
nseg_tcurr=zeros(nexp_tcurr,1);
for i=1:nexp_tcurr
    if nexp_tcurr==1
        xi_curr=nobs;
        nseg_tcurr=nobs;
    else
        if i==1
           nposs=nobs-nexp_tcurr*opt.tmin+1;
           xi_curr(i)=opt.tmin+unidrnd(nposs)-1;
           nseg_tcurr(i)=xi_curr(i);
        elseif i>1 && i<nexp_tcurr
           nposs=nobs-xi_curr(i-1)-opt.tmin*(nexp_tcurr-i+1)+1;
           xi_curr(i)=opt.tmin+unidrnd(nposs)+xi_curr(i-1)-1;
           nseg_tcurr(i)=xi_curr(i)-xi_curr(i-1);
        else
           xi_curr(i)=nobs;
           nseg_tcurr(i)=xi_curr(i)-xi_curr(i-1);	
        end
    end
end

ui_curr=zeros(nexp_ucurr,1);
nseg_ucurr=zeros(nexp_ucurr,1);
for j=1:nexp_ucurr
    if nexp_ucurr==1
        ui_curr=max(u);
        nseg_ucurr=nrep;
    else
        if j==1
           nposs=nu_unique-nexp_ucurr*opt.umin+1;
           ui_curr(j)=u_sort(opt.umin+unidrnd(nposs)-1);
           nseg_ucurr(j)=sum(u<=ui_curr(j));      
        elseif j>1 && j<nexp_ucurr
           nposs=nu_unique-find(u_sort==ui_curr(j-1))-opt.umin*(nexp_ucurr-j+1)+1;
           ui_curr(j)=u_sort(opt.umin+unidrnd(nposs)+find(u_sort==ui_curr(j-1))-1);
           nseg_ucurr(j)=sum(u<=ui_curr(j))-sum(u<=ui_curr(j-1));
        else
           ui_curr(j)=max(u);
           nseg_ucurr(j)=nrep-sum(u<=ui_curr(j-1));	
        end
    end
end

%tau sq
tau_up_limit=1000;
tausq_curr=rand(nexp_tcurr,nexp_ucurr)*tau_up_limit;

%nu mat hat
nfreq_hat=50;
freq_hat=(0:nfreq_hat)'/(2*nfreq_hat);
[nu_mat_hat]=lin_basis_func(freq_hat,nbeta);

%beta
beta_curr=zeros(nbeta,nexp_tcurr,nexp_ucurr);
for i=1:nexp_tcurr
    for j=1:nexp_ucurr
        [beta_mean, beta_var, ~, ~]=postbeta(i,j,nseg_tcurr(i),x,u,xi_curr,ui_curr,tausq_curr(i,j),opt);
        beta_curr(:,i,j)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'));
    end
end

%%%initialize .mat file
S.('x') = x;
S.('u') = u;
S.('opt') = opt;
S.('nexp_t') = [];
S.('nexp_u') = [];
S.('nseg_t') = {};
S.('nseg_u') = {};
S.('xi') = {};
S.('ui') = {};
S.('beta') = {};
S.('tausq') = {};
S.('log_spec_hat') = {};
if strcmp(opt.conv_diag,'on') 
    S.('conv_diag') = [];    
end

save(opt.fname, '-struct', 'S', '-v7.3');

%create matfile object to use in MCMC sampler
out = matfile(opt.fname,'Writable',true);

%create data structures to store temporary results in between outputs
beta_tmp = cell(opt.batchsize,1);
log_spec_hat_tmp = cell(opt.batchsize,1);  
nexp_t_tmp = zeros(opt.batchsize,1);
nexp_u_tmp = zeros(opt.batchsize,1);
nseg_t_tmp = cell(opt.batchsize,1);
nseg_u_tmp = cell(opt.batchsize,1);
tausq_tmp = cell(opt.batchsize,1);
ui_tmp = cell(opt.batchsize,1);
xi_tmp = cell(opt.batchsize,1);
if strcmp(opt.conv_diag,'on')
    conv_diag_tmp=zeros(opt.batchsize,nbeta+1);
end

%preallocate for worst case memory use
%will cause out of memory error if load is too much
for i=1:opt.batchsize
    beta_tmp{i} = zeros(nbeta, opt.nexp_tmax, opt.nexp_umax);
    log_spec_hat_tmp{i} = zeros(nfreq_hat+1, opt.nexp_tmax, opt.nexp_umax);  
    nexp_t_tmp(i) = opt.nexp_tmax;
    nexp_u_tmp(i) = opt.nexp_umax;
    nseg_t_tmp{i} = zeros(opt.nexp_tmax,1);
    nseg_u_tmp{i} = zeros(opt.nexp_umax,1);
    tausq_tmp{i} = zeros(opt.nexp_tmax,opt.nexp_umax);
    ui_tmp{i} = zeros(opt.nexp_umax,1);
    xi_tmp{i} = zeros(opt.nexp_tmax,1);
    if strcmp(opt.conv_diag,'on')
        conv_diag_tmp(i,:)=zeros(1,nbeta+1);
    end
end

%set batch index
batch_idx=1;
     
for p=1:opt.nloop
    
    if opt.Rev_Jump_t==1
        %BETWEEN MODEL MOVE
        %Number of available segments
        kk=length(find(nseg_tcurr>=2*opt.tmin));
        %Deciding on birth or death
        [nexp_prop,log_move_curr,log_move_prop]=move(kk,nexp_tcurr,opt.nexp_tmax);
        if nexp_prop<nexp_tcurr
            %btwn=-1;
            %Death
            [met_rat,nseg_prop,xi_prop,tau_prop,beta_prop]=...
                death_t(x,u,nexp_tcurr,nexp_ucurr,nexp_prop,...
                tausq_curr,xi_curr,ui_curr,...
                nseg_tcurr,beta_curr,...
                log_move_curr,log_move_prop,opt);
        elseif nexp_prop>nexp_tcurr
            %btwn=1;
            %Birth,
            [met_rat,nseg_prop,xi_prop,tau_prop,beta_prop]=birth_t(...
                x,u,nexp_tcurr,nexp_ucurr,nexp_prop,...
                tausq_curr,xi_curr,ui_curr,...
                nseg_tcurr,beta_curr,log_move_curr,log_move_prop,opt);
        else
            %btwn=0;
            xi_prop=xi_curr;
            nseg_prop=nseg_tcurr;
            tau_prop=tausq_curr;
            beta_prop=beta_curr;
            met_rat=1;
        end
        %Update current values
        r=rand;
        if r<met_rat
            nexp_tcurr=nexp_prop;
            xi_curr=xi_prop;
            nseg_tcurr=nseg_prop;
            beta_curr=beta_prop;
            tausq_curr=tau_prop;
        end   
    end
    %WITHIN MODEL MOVE
    %Drawing a new cut point and betas simultaneously
    %First draw the size of the move    
    [epsilon,xi_prop,beta_prop,nseg_time_new,~]=...
        within_t(x,u,nexp_tcurr,nexp_ucurr,xi_curr,ui_curr,...
        beta_curr,nseg_tcurr,tausq_curr,opt);    
    
    r=rand;
    if (r<epsilon|| p==1) 
        xi_curr=xi_prop;
        nseg_tcurr=nseg_time_new;
        beta_curr=beta_prop;
    else
        for j=1:nexp_tcurr
            for i=1:nexp_ucurr
                [beta_mean, beta_var, ~, ~]=postbeta(j,i,nseg_tcurr(j),x,u,xi_curr,ui_curr,tausq_curr(j,i),opt);
                beta_curr(:,j,i)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'));
            end
        end   
    end
    %Drawing tau
    for j=1:nexp_tcurr
        for i=1:nexp_ucurr
            tau_a=opt.nbasis/2+opt.tau_prior_a;
            tau_b=sum(beta_curr(2:nbeta,j,i).^2)/2+opt.tau_prior_b;
            tausq_curr(j,i)=1/gamrnd(tau_a,1/tau_b);
        end
    end 
         
    if opt.Rev_Jump_u==1
        %BETWEEN MODEL MOVE
        %Number of available segments (taking duplicate covariate values into
        %consideration)
        kk=0;
        for i=1:nexp_ucurr
            if i==1
                kk=kk+(sum(u_sort<=ui_curr(i))>=2*opt.umin);
            else
                kk=kk+(sum(u_sort<=ui_curr(i)&u_sort>ui_curr(i-1))>=2*opt.umin);
            end
        end
        %Deciding on birth or death
        [nexp_prop,log_move_curr,log_move_prop]=move(kk,nexp_ucurr,opt.nexp_umax);     
        if nexp_prop<nexp_ucurr
            %btwn=-1;
            %Death
            [met_rat,nseg_prop,ui_prop,tau_prop,beta_prop]=death_u(...
                x,u,nexp_tcurr,nexp_ucurr,nexp_prop,...
                tausq_curr,xi_curr,ui_curr,...
                nseg_tcurr,nseg_ucurr,beta_curr,log_move_curr,log_move_prop,opt);
        elseif nexp_prop>nexp_ucurr
            %btwn=1;
            %Birth,
            [met_rat,nseg_prop,ui_prop,tau_prop,beta_prop]=birth_u(...
                x,u,nexp_tcurr,nexp_ucurr,nexp_prop,...
                tausq_curr,xi_curr,ui_curr,...
                nseg_tcurr,nseg_ucurr,beta_curr,log_move_curr,log_move_prop,opt);
        else
            %btwn=0;
            ui_prop=ui_curr;
            nseg_prop=nseg_ucurr;
            tau_prop=tausq_curr;
            beta_prop=beta_curr;
            met_rat=1;
        end        
        r=rand;
        if r<met_rat
            nexp_ucurr=nexp_prop;
            ui_curr=ui_prop;
            nseg_ucurr=nseg_prop;
            tausq_curr=tau_prop;
            beta_curr=beta_prop;
        end
    end
    %WITHIN MODEL MOVE
    %Drawing a new cut point and betas simultaneously
    %First draw the size of the move 
    [epsilon,ui_prop,beta_prop,nseg_cov_new,~]=...
        within_u(x,u,nexp_tcurr,nexp_ucurr,...
        xi_curr,ui_curr,beta_curr,nseg_tcurr,nseg_ucurr,tausq_curr,opt);    
    
    r=rand;   
    if (r<epsilon|| p==1) 
        ui_curr=ui_prop;
        nseg_ucurr=nseg_cov_new;
        beta_curr=beta_prop;
    else
       for j=1:nexp_tcurr
            for i=1:nexp_ucurr
                [beta_mean, beta_var, ~, ~]=postbeta(j,i,nseg_tcurr(j),x,u,xi_curr,ui_curr,tausq_curr(j,i),opt);
                beta_curr(:,j,i)=mvnrnd(beta_mean,0.5*(beta_var+beta_var'));
            end
        end   
    end
    %Drawing tau
    for j=1:nexp_tcurr
        for i=1:nexp_ucurr
            tau_a=opt.nbasis/2+opt.tau_prior_a;
            tau_b=sum(beta_curr(2:nbeta,j,i).^2)/2+opt.tau_prior_b;
            tausq_curr(j,i)=1/gamrnd(tau_a,1/tau_b);
        end
    end
    %Estimating Spectral Density
    log_spec_hat_curr=zeros(nfreq_hat+1,nexp_tcurr,nexp_ucurr);
    for j=1:nexp_tcurr
        for i=1:nexp_ucurr
            log_spec_hat_curr(:,j,i)=nu_mat_hat*beta_curr(:,j,i);
        end
    end
    
    %Create convergence diagnostics
    if strcmp(opt.conv_diag,'on')  
        conv_diag_curr=convdiag(nobs,nu_unique,nbeta,...
            nseg_tcurr,nseg_ucurr,nexp_tcurr, nexp_ucurr,...
            beta_curr, tausq_curr);
    end
    
    %Output to temporary container
    beta_tmp{batch_idx} = beta_curr;
    log_spec_hat_tmp{batch_idx} = log_spec_hat_curr;  
    nexp_t_tmp(batch_idx) = nexp_tcurr;
    nexp_u_tmp(batch_idx) = nexp_ucurr;
    nseg_t_tmp{batch_idx} = nseg_tcurr;
    nseg_u_tmp{batch_idx} = nseg_ucurr;
    tausq_tmp{batch_idx} = tausq_curr;
    ui_tmp{batch_idx} = ui_curr;
    xi_tmp{batch_idx} = xi_curr;
    if strcmp(opt.conv_diag,'on')
        conv_diag_tmp(batch_idx,:)=conv_diag_curr;
    end   
    
    if(mod(p,opt.disp_iter)==0 && p~=opt.nloop)||(p==opt.nloop && opt.disp_iter~=0)
        %display the current iteration
        disp(['Current MCMC iteration for ', opt.fname,': ' num2str(p)]);
    elseif (p==opt.nloop) 
    end
    
    %once batch is complete, output to .mat file
    if(mod(p,opt.batchsize)==0)
                
        %output data to .mat file
        out.beta((p-opt.batchsize+1):p,1) = beta_tmp;
        out.log_spec_hat((p-opt.batchsize+1):p,1) = log_spec_hat_tmp;    
        out.nexp_t((p-opt.batchsize+1):p,1)=nexp_t_tmp;
        out.nexp_u((p-opt.batchsize+1):p,1)=nexp_u_tmp;
        out.nseg_t((p-opt.batchsize+1):p,1)=nseg_t_tmp;
        out.nseg_u((p-opt.batchsize+1):p,1)=nseg_u_tmp;
        out.tausq((p-opt.batchsize+1):p,1) = tausq_tmp;
        out.ui((p-opt.batchsize+1):p,1) = ui_tmp;
        out.xi((p-opt.batchsize+1):p,1) = xi_tmp;
        if strcmp(opt.conv_diag,'on') 
            out.conv_diag((p-opt.batchsize+1):p,1:nbeta+1)=conv_diag_tmp;
        end  
        
        %reset temporary data containers        
        beta_tmp = cell(opt.batchsize,1);
        log_spec_hat_tmp = cell(opt.batchsize,1);  
        nexp_t_tmp = zeros(opt.batchsize,1);
        nexp_u_tmp = zeros(opt.batchsize,1);
        nseg_t_tmp = cell(opt.batchsize,1);
        nseg_u_tmp = cell(opt.batchsize,1);
        tausq_tmp = cell(opt.batchsize,1);
        ui_tmp = cell(opt.batchsize,1);
        xi_tmp = cell(opt.batchsize,1);
        if strcmp(opt.conv_diag,'on')
            conv_diag_tmp=zeros(opt.batchsize,nbeta+1);
        end
       
        %reset batch index
        batch_idx=0;
        
    end
    
    %increment batch index
    batch_idx = batch_idx+1;
    
    %if final iteration is not multiple of batch size, output remaining to .mat file
    if (p==opt.nloop && mod(p,opt.batchsize)~=0)         
       %output data to .mat file
        out.beta((p-mod(p,opt.batchsize)+1):p,1) = beta_tmp(1:mod(p,opt.batchsize));
        out.log_spec_hat((p-mod(p,opt.batchsize)+1):p,1) = log_spec_hat_tmp(1:mod(p,opt.batchsize));    
        out.nexp_t((p-mod(p,opt.batchsize)+1):p,1)=nexp_t_tmp(1:mod(p,opt.batchsize));
        out.nexp_u((p-mod(p,opt.batchsize)+1):p,1)=nexp_u_tmp(1:mod(p,opt.batchsize));
        out.nseg_t((p-mod(p,opt.batchsize)+1):p,1)=nseg_t_tmp(1:mod(p,opt.batchsize));
        out.nseg_u((p-mod(p,opt.batchsize)+1):p,1)=nseg_u_tmp(1:mod(p,opt.batchsize));
        out.tausq((p-mod(p,opt.batchsize)+1):p,1) = tausq_tmp(1:mod(p,opt.batchsize));
        out.ui((p-mod(p,opt.batchsize)+1):p,1) = ui_tmp(1:mod(p,opt.batchsize));
        out.xi((p-mod(p,opt.batchsize)+1):p,1) = xi_tmp(1:mod(p,opt.batchsize));
        if strcmp(opt.conv_diag,'on')  
            out.conv_diag((p-mod(p,opt.batchsize)+1):p,1:nbeta+1)=conv_diag_tmp(1:mod(p,opt.batchsize),:);
        end          
    end
end

disp(['Results of MCMC sampler output to ' opt.fname, '.mat']);

