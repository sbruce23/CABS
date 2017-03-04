function [collapsed_out] = CABS_collapsed(out,bands,varargin)
%Function produces pointwise credible intervals for log spectrum and any specified collapsed measures
    %collapsed_out contains the pointwise credible intervals and point estimates for
    %each collapsed measure indexed by time and subject.
    
    %default settings
    param = struct('alpha',0.05,...%95% pointwise credible intervals
                   'plots','on',... %turn on plot output
                   'samp_iter',1000); %max number of samples from iterations after burn-in
    
    %handle changes to settings
    if nargin > 2
        op_names={'alpha','plots','samp_iter'};
        if ischar(varargin{1}) %match and set arguments by names
            if nargin > 8
                error('Too many input arguments.');
            elseif mod(nargin, 2)==0
                for i=1:2:(nargin-2)
                    id = find(strcmp(op_names,varargin{i}));
                    if isempty(id)
                        error('%s is an invalid name. Option names are alpha, plots, or samp_iter.',varargin{i});
                    else
                       param = setVal(param,op_names{id}, varargin{i+1});
                    end
                end
            else
                error('Either option name or its corresponding value is missing.');
            end
        else %set arguments directly
            if nargin <= 5
                for i=1:(nargin-2)
                    param = setVal(param,op_names{i}, varargin{i});
                end
            else
                error('Too many input arguments.');
            end
        end   
    end

    if size(out,1)==1 %only one chain
        %pull in data from matfile for chain
        xi=out.xi;
        ui=out.ui;
        nseg_t=out.nseg_t;
        nexp_t=out.nexp_t;
        log_spec_hat=out.log_spec_hat;
        x=out.x;
        u=out.u;
        opt=out.opt;

        %initialize cells
        collapsed_out=cell(size(bands,1),3);
        collapsed=cell(size(bands,1),1);
        
        %sample down iterations for computational efficiency
        rng(opt.seed);
        p_set=sort(randsample((opt.nwarmup+1):opt.nloop,min(param.samp_iter,opt.nloop-opt.nwarmup)));
        
        %initialize matrices
        for k=1:size(bands,1)
            collapsed{k}=zeros(1,size(x,1),size(x,2),size(p_set,2));
        end
        
        %collect spectral estimates across all iterations after warmup
        freq_hat=linspace(0,0.5,51);
        for k=1:size(x,2)
            for p=1:size(p_set,2)
                xi_curr=xi{p_set(p)};
                nseg_t_curr=nseg_t{p_set(p)};       
                log_spec_hat_curr=log_spec_hat{p_set(p)}(:,:,find(u(1,k)<=ui{p_set(p)},1));
                log_spec_hat_new=zeros(size(freq_hat,2),size(x,1));
                %obtain spectral estimator
                for g=1:nexp_t(p_set(p))
                    if g==1
                        log_spec_hat_new(:,1:xi_curr(g))=...
                            repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));  
                    else
                        log_spec_hat_new(:,(xi_curr(g-1)+1):xi_curr(g))=...
                            repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));
                    end
                end
                %obtain collapsed measures
                if size(bands,1)>0
                    for h=1:size(bands,1)
                        if isequal(size(bands{h}),[1 2])
                            idx_b1=find(ismember(freq_hat,bands{h}));
                            collapsed{h}(:,:,k,p)=...
                                sum(exp(log_spec_hat_new(idx_b1(1):idx_b1(2),:)),1);
                        elseif isequal(size(bands{h}),[2 2])
                            idx_b1=find(ismember(freq_hat,bands{h}(1,:)));
                            idx_b2=find(ismember(freq_hat,bands{h}(2,:)));                           
                            collapsed{h}(:,:,k,p)=...
                                sum(exp(log_spec_hat_new(idx_b1(1):idx_b1(2),:)),1)./...
                                    sum(exp(log_spec_hat_new(idx_b2(1):idx_b2(2),:)),1);
                        else
                            error(['bands{',num2str(h),'} should be a 1x2 or 2x2 double.'])
                        end
                    end
                end
            end
        end
        
        %obtain lower and upper bounds for credible intervals and mean
        for k=1:size(collapsed,1)
            collapsed_out{k,1}=quantile(collapsed{k},param.alpha/2,4);
            collapsed_out{k,2}=mean(collapsed{k},4);
            collapsed_out{k,3}=quantile(collapsed{k},1-(param.alpha/2),4);
        end               
        
    elseif size(out,1) > 1 %multiple chains

        %pull in data from matfile for chain
        x=out{1}.x;
        u=out{1}.u;
        opt=out{1}.opt;
        
        %set seed for down sampling
        rng(opt.seed);
        
        %initialize cells
        collapsed_out=cell(size(bands,1),3);
        collapsed=cell(size(bands,1),size(out,1));

        %collect spectral estimates across all iterations after warmup
        freq_hat=linspace(0,0.5,51);
        for m=1:size(out,1)
            
            %pull in data from matfile for chain m
            xi=out{m}.xi;
            ui=out{m}.ui;
            nseg_t=out{m}.nseg_t;
            nexp_t=out{m}.nexp_t;
            log_spec_hat=out{m}.log_spec_hat;
            opt=out{m}.opt;
            
            %sample down iterations for computational efficiency
            p_set=sort(randsample((opt.nwarmup+1):opt.nloop,min(ceil(param.samp_iter/size(out,1)),opt.nloop-opt.nwarmup)));
             
            %initialize matrices
            for k=1:size(bands,1)
                collapsed{k,m}=zeros(1,size(x,1),size(x,2),size(p_set,2));
            end        

            for k=1:size(x,2)
                for p=1:size(p_set,2)
                    xi_curr=xi{p_set(p)};
                    nseg_t_curr=nseg_t{p_set(p)};       
                    log_spec_hat_curr=log_spec_hat{p_set(p)}(:,:,find(u(k)<=ui{p_set(p)},1));
                    log_spec_hat_new=zeros(size(freq_hat,2),size(x,1));
                    %obtain spectral estimator
                    for g=1:nexp_t(p_set(p))
                        if g==1
                            log_spec_hat_new(:,1:xi_curr(g))=...
                                repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));  
                        else
                            log_spec_hat_new(:,(xi_curr(g-1)+1):xi_curr(g))=...
                                repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));
                        end
                    end
                    %obtain collapsed measures
                    if size(bands,1)>0
                        for h=1:size(bands,1)
                            if isequal(size(bands{h}),[1 2])
                                idx_b1=find(ismember(freq_hat,bands{h}));
                                collapsed{h,m}(:,:,k,p)=...
                                    sum(exp(log_spec_hat_new(idx_b1(1):idx_b1(2),:)),1);
                            elseif isequal(size(bands{h}),[2 2])
                                idx_b1=find(ismember(freq_hat,bands{h}(1,:)));
                                idx_b2=find(ismember(freq_hat,bands{h}(2,:)));                           
                                collapsed{h,m}(:,:,k,p)=...
                                    sum(exp(log_spec_hat_new(idx_b1(1):idx_b1(2),:)),1)./...
                                        sum(exp(log_spec_hat_new(idx_b2(1):idx_b2(2),:)),1);
                            else
                                error(['bands{',num2str(h),'} should be a 1x2 or 2x2 double.'])
                            end
                        end
                    end                   
                end
            end
        end   
        
        %obtain lower and upper bounds for credible intervals
        for k=1:size(collapsed,1)
            collapsed_out{k,1}=quantile(cat(4,collapsed{k,:}),param.alpha/2,4);
            collapsed_out{k,2}=mean(cat(4,collapsed{k,:}),4);
            collapsed_out{k,3}=quantile(cat(4,collapsed{k,:}),1-(param.alpha/2),4);
        end      
    end
    
    %plotting
    if strcmp(param.plots,'on')
        for k=1:size(x,2)
            for h=1:size(bands,1)
                figure
                plot(1:size(x,1),collapsed_out{h,1}(:,:,k),'b:',...
                     1:size(x,1),collapsed_out{h,2}(:,:,k),'b',...
                     1:size(x,1),collapsed_out{h,3}(:,:,k),'b:');
                title([num2str((1-param.alpha)*100),'% Credible Intervals for Collapsed Measure ',num2str(h),...
                       ' for Time Series ',num2str(k)]);
                xlabel('Time')
            end
        end
        
        if isequal(u,sort(u, 'ascend')) %covariate and time series already sorted by scaled covariate
            for h=1:size(bands,1)
                figure;hold on
                surf_LB=surf(1:size(x,1),sort(u,'ascend'),reshape(collapsed_out{h,1},[size(x,1),size(x,2)])');
                set(surf_LB, 'edgecolor','none')
                surf_UB=surf(1:size(x,1),sort(u,'ascend'),reshape(collapsed_out{h,3},[size(x,1),size(x,2)])');
                set(surf_UB, 'edgecolor','none')
                set(gca,'YDir','reverse')
                title([num2str((1-param.alpha)*100),'% Credible Intervals for for Collapsed Measure ',num2str(h)])
                xlabel('Time')
                ylabel('Covariate')
                grid on
                view(-40,30)                  
            end
        else %sort collapsed_out by scaled covariate then plot
            [~,idx_u]=sort(u, 'ascend');
            for h=1:size(bands,1)
                figure;hold on
                surf_LB=surf(1:size(x,1),sort(u,'ascend'),reshape(collapsed_out{h,1}(:,:,idx_u),[size(x,1),size(x,2)])');
                set(surf_LB, 'edgecolor','none')
                surf_UB=surf(1:size(x,1),sort(u,'ascend'),reshape(collapsed_out{h,3}(:,:,idx_u),[size(x,1),size(x,2)])');
                set(surf_UB, 'edgecolor','none')
                set(gca,'YDir','reverse')
                title([num2str((1-param.alpha)*100),'% Credible Intervals for for Collapsed Measure ',num2str(h)])
                xlabel('Time')
                ylabel('Covariate')
                grid on
                view(-40,30)                  
            end
        end
    end        
    
end

