%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo file for "Adaptive Bayesian Spectral Analysis of Nonstationary 
% Biomedical Time Series" by Bruce, Hall, Buysse, and Krafty (2017)
%
% Author: Scott A. Bruce
% 
% Date: February 27, 2017
% 
% Description:
% Reproduces simulated examples from "Adaptive Bayesian Spectral Analysis 
% of Nonstationary Biomedical Time Series" by Bruce, Hall, Buysse, and
% Krafty (2017)
%
% Contents:
% 1) Piecewise AR Process
%   1a) Simulate data for piecewise AR process
%   1b) Specify options for MCMC sampler
%   1c) Run CABS procedure to obtain spectral estimator
%   1d) Calculate the true log spectrum for piecewise AR process
%   1e) Plot true and estimated time-varying log spectrum for each replicate
%   1f) Estimate and plot point estimates and credible intervals for spectrum and collapsed measures
%
% 2) Slowly Varying AR Process
%   2a) Simulate data for slowly varying AR process
%   2b) Specify options for MCMC sampler
%   2c) Run CABS procedure to obtain spectral estimator
%   2d) Calculate the true log spectrum for slowly varying AR process
%   2e) Plot true and estimated time-varying log spectrum for each replicate
%   2f) Estimate and plot point estimates and credible intervals for spectrum and collapsed measures
%
% 3) Convergence Diagnostics
%   3a) Simulate data for slowly varying AR process
%   3b) Specify options for MCMC sampler (5 chains)
%   3c) Run CABS procedure to obtain spectral estimator (5 chains)
%   3d) Calculate convergence diagnostics
%   3e) Create convergence diagnostic plots
%
% BEFORE YOU BEGIN remember to change the current working folder to the 
% 'programs' folder containing the MATLAB functions need to run the method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Piecewise AR Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% 1a) Simulate data for piecewise AR process
change='abrupt'; %specifies type of change in power spectrum
i=8; %number of time series
j=1000; %length of time series
seed1=89156; %seed for data simulation
seed2=32486; %seed for MCMC sampler
[x,u]=simdata(change,i,j,seed1);

% 1b) Set options for MCMC sampler
opt=setMCMCOptions('nloop',5000,'seed', seed2, 'nwarmup',1000,...
               'nexp_tmax',10,'nexp_umax',8,'tmin',50,'umin',1,...
               'nbasis',7,'fname','CABSout_abrupt','batchsize',1000,...
               'disp_plots','on','disp_iter',1,'parcomp','on');%parcomp on requires parallel computing toolbox

% 1c) Run CABS method to produce sampler results, plots, and spectral estimator
tic         
[spec_est, out]=CABS(x,u,opt);
toc
%Note: To reproduce the exact results from the paper, use the following
%code in place of step 1c:
% out=cell(8,1);
% for k=1:8
%     out{k}=matfile(['C:\CABSdemo\datasets\CABSout_abrupt_',num2str(k),'.mat'],'Writable',true);
% end
% spec_est=CABS_estimator(out);
% CABS_plot(out);

% 1d) True log spectrum
spec_true=zeros(51,j,i);
ar1 =eye(2);
ar2=[0 0;0 0];
sig = [1 0; 0 1];
freq_hat=linspace(0,0.5,51);
for k=1:i
    if strcmp(change,'slow') %Slowly varying spectrum
        if u(k) <= 0.5
            a= -0.5 + 1*(1:j)/j; 
        else
             a= -0.9 + 1.8*(1:j)/j; 
        end   
    elseif strcmp(change,'abrupt') %Piecewise AR spectrum
        if u(k) <= 0.5
            a= [repmat(-0.5,j/2,1); repmat(0.5,j/2,1)];
        else
            a= [repmat(-0.9,j/2,1); repmat(0.9,j/2,1)]; 
        end
    end
    for t=1:j
        spect_true=ARspec([a(t)*ar1 ar2],sig,freq_hat);
        spec_true(:,t,k) = log(real(squeeze(real(spect_true(1,1,:)))));
    end
end

% 1e) Plot true and estimated time-varying log spectrum for each replicate
for k=1:i
    figure;
    surf_true=surf(1:j,freq_hat,spec_true(:,:,k));
    set(surf_true, 'edgecolor','none')
    set(gca,'YDir','reverse')
    title(['True Log Spectrum for Time Series ',num2str(k)])
    xlim([0 1000])
    ylim([0 0.5])
    zlim([-2 5])
    ax=gca;
    ax.XTick = [0 200 400 600 800 1000];
    ax.YTick = [0 0.1 0.2 0.3 0.4 0.5];
    ax.ZTick = [-2 -1 0 1 2 3 4 5];
   
    figure;
    surf_est=surf(1:j,freq_hat,spec_est(:,:,k));
    set(surf_est, 'edgecolor','none')
    set(gca,'YDir','reverse')
    title(['Estimated Log Spectrum for Time Series ',num2str(k)])
    xlim([0 1000])
    ylim([0 0.5])
    zlim([-2 5])
    ax=gca;
    ax.XTick = [0 200 400 600 800 1000];
    ax.YTick = [0 0.1 0.2 0.3 0.4 0.5];
    ax.ZTick = [-2 -1 0 1 2 3 4 5];
end

% 1f) Estimate and plot credible intervals for spectrum and collapsed measures
bands=cell(3,1);
bands{1} = [.15, .40];               %Define the HF band
bands{2} = [.15, .40; .04, .40];     %Define the HF normalized
bands{3} = [.05, .15; .15, .40];     %Define the LF/HF ratio
[collapsed_out] = CABS_collapsed(out,bands);
%collapsed_out contain the pointwise credible intervals and point estimates for
%each collapsed measure indexed by time and subject.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Slowly varying AR Process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% 2a) Simulate data for slowly varying AR process
change='slow';
i=8;
j=1000;
seed1=98165;
seed2=31085;
[x,u]=simdata(change,i,j,seed1);

% 2b) Set options for MCMC sampler
opt=setMCMCOptions('nloop',5000,'seed', seed2, 'nwarmup',1000,...
               'nexp_tmax',10,'nexp_umax',8,'tmin',50,'umin',1,...
               'nbasis',7,'fname','CABSout_smooth','batchsize',1000,...
               'disp_plots','on','disp_iter',1,'parcomp','on');%parcomp on requires parallel computing toolbox

% 2c) Run CABS method to produce sampler results, plots, and spectral estimator
tic         
[spec_est, out]=CABS(x,u,opt);
toc
%Note: To reproduce the exact results from the paper, use the following
%code in place of step 2c:
% out=cell(8,1);
% for k=1:8
%     out{k}=matfile(['C:\CABSdemo\datasets\CABSout_smooth_',num2str(k),'.mat'],'Writable',true);
% end
% spec_est=CABS_estimator(out);
% CABS_plot(out);

% 2d) True log spectrum
spec_true=zeros(51,j,i);
ar1 =eye(2);
ar2=[0 0;0 0];
sig = [1 0; 0 1];
freq_hat=linspace(0,0.5,51);
for k=1:i
    if strcmp(change,'slow') %Slowly varying spectrum
        if u(k) <= 0.5
            a= -0.5 + 1*(1:j)/j; 
        else
             a= -0.9 + 1.8*(1:j)/j; 
        end   
    elseif strcmp(change,'abrupt') %Piecewise AR spectrum
        if u(k) <= 0.5
            a= [repmat(-0.5,j/2,1); repmat(0.5,j/2,1)];
        else
            a= [repmat(-0.9,j/2,1); repmat(0.9,j/2,1)]; 
        end
    end
    for t=1:j
        spect_true=ARspec([a(t)*ar1 ar2],sig,freq_hat);
        spec_true(:,t,k) = log(real(squeeze(real(spect_true(1,1,:)))));
    end
end

% 2e) Plot true and estimated time-varying log spectrum for each replicate
for k=1:i
    figure;
    surf_true=surf(1:j,freq_hat,spec_true(:,:,k));
    set(surf_true, 'edgecolor','none')
    set(gca,'YDir','reverse')
    title(['True Log Spectrum for Time Series ',num2str(k)])
    xlim([0 1000])
    ylim([0 0.5])
    zlim([-2 5])
    ax=gca;
    ax.XTick = [0 200 400 600 800 1000];
    ax.YTick = [0 0.1 0.2 0.3 0.4 0.5];
    ax.ZTick = [-2 -1 0 1 2 3 4 5];
   
    figure;
    surf_est=surf(1:j,freq_hat,spec_est(:,:,k));
    set(surf_est, 'edgecolor','none')
    set(gca,'YDir','reverse')
    title(['Estimated Log Spectrum for Time Series ',num2str(k)])
    xlim([0 1000])
    ylim([0 0.5])
    zlim([-2 5])
    ax=gca;
    ax.XTick = [0 200 400 600 800 1000];
    ax.YTick = [0 0.1 0.2 0.3 0.4 0.5];
    ax.ZTick = [-2 -1 0 1 2 3 4 5];
end

% 2f) Estimate and plot credible intervals for spectrum and collapsed measures
bands=cell(3,1);
bands{1} = [.15, .40];               %Define the HF band
bands{2} = [.15, .40; .04, .40];     %Define the HF normalized
bands{3} = [.05, .15; .15, .40];     %Define the LF/HF ratio
[collapsed_out] = CABS_collapsed(out,bands);
%collapsed_out contain the pointwise credible intervals and point estimates for
%each collapsed measure indexed by time and subject.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Convergence Diagnostics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% 3a) Simulate data for slowly varing AR process
change='slow';
i=8;
j=1000;
seed1=69814;
seed2=[23785 19872 86744 7860 47845];
out=cell(5,1);
[x,u]=simdata(change,i,j,seed1);

% WARNING: Steps 3b and 3c may take a very long time to run 
parfor k=1:5 %change 'parfor' to 'for' for serial computing
    % 3b) Set options for MCMC sampler
    opt=setMCMCOptions('nloop',5000,'seed', seed2(k), 'nwarmup',1000,...
               'nexp_tmax',10,'nexp_umax',8,'tmin',50,'umin',1,...
               'nbasis',7,'fname', ['CABSout_conv_',num2str(k)],'batchsize',1000,...
               'conv_diag','on'); %last option turns on 
                                  %tracking of convergence diagnostic measures
    % 3c) Run CABS method
    tic         
    out{k}=CABS(x,u,opt);
    toc
end
% Alternatively, forego steps 3b and 3c and link to output file already created using:
% for k=1:5
%     out{k} = matfile(['C:\CABSdemo\datasets\CABSout_conv_',num2str(k),'.mat'],'Writable',true);
% end

% 3d) Calculate convergence diagnostics
C=5; %number of chains run
T=1000; %length of chains
P=9; %number of sets of parameters tracking 
     %(equals number of beta coefficients (basis + 1) + 1 for tausq)
b=5; %initial batch size
Q=floor(T/(2*b)); %number of sweeps to use

%initialize data containers
V_theta=zeros(Q,P);
W_C_theta=zeros(Q,P);
W_M_theta=zeros(Q,P);
W_M_W_C_theta=zeros(Q,P);

V_theta_m=zeros(P,P,Q);
W_C_theta_m=zeros(P,P,Q);
W_M_theta_m=zeros(P,P,Q);
W_M_W_C_theta_m=zeros(P,P,Q);

PSRF1=zeros(Q,P);
PSRF2=zeros(Q,P);
MPSRF1=zeros(Q,1);
MPSRF2=zeros(Q,1);

V_theta_maxeig=zeros(Q,1);
W_C_theta_maxeig=zeros(Q,1);
W_M_theta_maxeig=zeros(Q,1);
W_M_W_C_theta_maxeig=zeros(Q,1);

%read needed parts into memory for processing
conv_diag_curr = cell(C,1);
nexp_tcurr = cell(C,1);
nexp_ucurr = cell(C,1);

for k=1:C
    conv_diag_curr{k}=out{k}.conv_diag;
    nexp_tcurr{k}=out{k}.nexp_t;
    nexp_ucurr{k}=out{k}.nexp_u;
end

%calculate convergence diagnostics for each batch
for q=5:Q
    
    disp(['batch=' num2str(q)]);
    
    theta_all=zeros(q*b*C,P);
    mod_all=zeros(q*b*C,3);

    for i=1:C
        theta_all((q*b*(i-1)+1):(q*b*i),:) = conv_diag_curr{i}((q*b+1):(2*q*b),:);
        mod_all((q*b*(i-1)+1):(q*b*i),:) = [repmat(i,q*b,1)...
                                            nexp_tcurr{i}((q*b+1):(2*q*b),:)...
                                            nexp_ucurr{i}((q*b+1):(2*q*b),:)];
    end

    M=size(unique(mod_all(:,2:3),'rows'),1); %number of different models visited by any chain
    
    %grand average
    theta_bar=mean(theta_all,1);
    
    %chain averages
    theta_Cbar=zeros(C,P+1);
    for i=1:C
        idx=find(mod_all(:,1)==i);
        theta_Cbar(i,:) = [i mean(theta_all(idx,:),1)];
    end
    
    %model averages
    tmp=unique(mod_all(:,2:3),'rows');
    theta_Mbar=zeros(M,P+2);
    for i=1:M
        idx=find(ismember(mod_all(:,2:3),tmp(i,:),'rows'));
        theta_Mbar(i,:) = [tmp(i,:) mean(theta_all(idx,:),1)];
    end
    
    %chain x model averages
    tmp=unique(mod_all,'rows');
    theta_CMbar=zeros(size(tmp,1),P+3);
    for i=1:size(tmp,1)
        idx=find(ismember(mod_all,tmp(i,:),'rows'));
        theta_CMbar(i,:)=[tmp(i,:) mean(theta_all(idx,:),1)];
    end
    
    %variance components 
    for i=1:size(theta_all,1)
        tmp1=(theta_all(i,:)-theta_bar).^2/(C*q*b-1);
        tmp2=(theta_all(i,:)-theta_bar)'*(theta_all(i,:)-theta_bar)/(C*q*b-1);
        V_theta(q,:)=V_theta(q,:)+tmp1;
        V_theta_m(:,:,q)=V_theta_m(:,:,q)+tmp2;
        
        idx=find(theta_Cbar(:,1)==mod_all(i,1));
        tmp1=(theta_all(i,:)-theta_Cbar(idx,2:end)).^2/(C*(q*b-1));        
        tmp2=(theta_all(i,:)-theta_Cbar(idx,2:end))'*(theta_all(i,:)-theta_Cbar(idx,2:end))/(C*(q*b-1));
        W_C_theta(q,:)=W_C_theta(q,:)+tmp1;
        W_C_theta_m(:,:,q)=W_C_theta_m(:,:,q)+tmp2;
        
        idx=find(ismember(theta_Mbar(:,1:2),mod_all(i,2:3),'rows'));
        tmp1=(theta_all(i,:)-theta_Mbar(idx,3:end)).^2/(C*q*b-M);
        tmp2=(theta_all(i,:)-theta_Mbar(idx,3:end))'*(theta_all(i,:)-theta_Mbar(idx,3:end))/(C*q*b-M);
        W_M_theta(q,:)=W_M_theta(q,:)+tmp1;
        W_M_theta_m(:,:,q)=W_M_theta_m(:,:,q)+tmp2;
        
        idx=find(ismember(theta_CMbar(:,1:3),mod_all(i,1:3),'rows'));
        tmp1=(theta_all(i,:)-theta_CMbar(idx,4:end)).^2/(C*(q*b-M));
        tmp2=(theta_all(i,:)-theta_CMbar(idx,4:end))'*(theta_all(i,:)-theta_CMbar(idx,4:end))/(C*(q*b-M));
        W_M_W_C_theta(q,:)=W_M_W_C_theta(q,:)+tmp1;
        W_M_W_C_theta_m(:,:,q)=W_M_W_C_theta_m(:,:,q)+tmp2;
    end   
    
    %potential scale reduction factors
    PSRF1(q,:)=V_theta(q,:)./W_C_theta(q,:);
    PSRF2(q,:)=W_M_theta(q,:)./W_M_W_C_theta(q,:);     
    MPSRF1(q)=real(max(eig(W_C_theta_m(:,:,q)\V_theta_m(:,:,q))));
    MPSRF2(q)=real(max(eig(W_M_W_C_theta_m(:,:,q)\W_M_theta_m(:,:,q))));
    
    %max eigenvals
    V_theta_maxeig(q)=real(max(eig(V_theta_m(:,:,q))));
    W_C_theta_maxeig(q)=real(max(eig(W_C_theta_m(:,:,q)))); 
    W_M_theta_maxeig(q)=real(max(eig(W_M_theta_m(:,:,q)))); 
    W_M_W_C_theta_maxeig(q)=real(max(eig(W_M_W_C_theta_m(:,:,q)))); 
   
end

% 3e) Create convergence diagnostic plots
figure; 
plot(5:Q,MPSRF1(5:Q)); 
ylim([1 2]);
hline = refline([0 1.2]);
hline.Color = 'r';
title('MPSRF_1 vs. q');
%MPSRF1 has settled down close to 1 (less than 1.2)

figure; 
plot(5:Q,MPSRF2(5:Q)); 
ylim([1 2]);
hline = refline([0 1.2]);
hline.Color = 'r';
title('MPSRF_2 vs. q');
%MPSRF2 has settled down close to 1 (less than 1.2)

figure; 
hold on; 
    plot(5:Q,V_theta_maxeig(5:Q));
    plot(5:Q,W_C_theta_maxeig(5:Q)); 
hold off; 
title('Max Eigenvalue for Vhat and W_c vs. q');
legend('Vhat','W_c','Location','northeast');

figure; 
hold on; 
    plot(5:Q,W_M_theta_maxeig(5:Q)); 
    plot(5:Q,W_M_W_C_theta_maxeig(5:Q)); 
hold off; 
title('Max Eigenvalue for W_m and W_mW_c vs. q');
legend('W_m','W_mW_c','Location','northeast');
%Max eigenvalue plots have not yet settled to a common value
%consider increasing number of iterations of the sampler