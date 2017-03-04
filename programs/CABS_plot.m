function CABS_plot(out)
%Function produces summary plots for CABS estimation procedure

    if size(out,1)==1 %only one chain
        opt=out.opt;
        xi=out.xi;
        ui=out.ui;
        nexp_t=out.nexp_t;
        nexp_u=out.nexp_u;
        
        post_prob_nexp_t=zeros(opt.nexp_tmax,1);
        for k=1:opt.nexp_tmax
            kk=find(nexp_t(opt.nwarmup+1:opt.nloop)==k);
            post_prob_nexp_t(k)=length(kk)/(opt.nloop-opt.nwarmup);
            if ~isempty(kk) && k>1
                xi_mat=zeros(length(kk),k);
                for g=1:length(kk)
                    xi_mat(g,:)=xi{kk(g)+opt.nwarmup}';
                end
                figure
                hold
                title(['Time Partition Points for a Partition of ',int2str(k)])
                for g=1:k-1
                    plot(xi_mat(:,g))                    
                end
                for g=1:k-1
                    figure
                    histogram(xi_mat(:,g),'FaceAlpha', 1,'Normalization','probability')
                    title(['Distribution of ',int2str(g), ' Time Partition Point for a mixture of ',int2str(k)])
                end
            end
        end
        figure
        bar(post_prob_nexp_t)
        title('Posterior Probability of Number of Time Segments') 

        post_prob_nexp_u=zeros(opt.nexp_umax,1);
        for k=1:opt.nexp_umax
            kk=find(nexp_u(opt.nwarmup+1:opt.nloop)==k);
            post_prob_nexp_u(k)=length(kk)/(opt.nloop-opt.nwarmup);
            if ~isempty(kk) && k>1
                ui_mat=zeros(length(kk),k);
                for g=1:length(kk)
                    ui_mat(g,:)=ui{kk(g)+opt.nwarmup}';
                end
                figure
                hold
                title(['Covariate Partition Points for a Partition of ',int2str(k)])
                for g=1:k-1
                    plot(ui_mat(:,g))                    
                end
                for g=1:k-1
                    figure
                    histogram(ui_mat(:,g),'FaceAlpha', 1,'Normalization','probability')
                    title(['Distribution of ',int2str(g), ' Covariate Partition Point for a mixture of ',int2str(k)])
                end
            end
        end
        figure
        bar(post_prob_nexp_u)
        title('Posterior Probability of Number of Covariate Segments') 

    elseif size(out,1) > 1 %multiple chains

        opt=out{1}.opt;
        
        post_prob_nexp_t=zeros(opt.nexp_tmax,1);
        xi_cell=cell(size(out,1),opt.nexp_tmax);
        
        for k=1:opt.nexp_tmax
            den=0;
            for m=1:size(out,1)
                opt=out{m}.opt;
                den=den+opt.nloop-opt.nwarmup;
                xi=out{m}.xi;
                kk=find(out{m}.nexp_t(opt.nwarmup+1:end,1)==k);
                post_prob_nexp_t(k)=post_prob_nexp_t(k)+length(kk);
                if ~isempty(kk) && k>1
                    xi_mat=zeros(length(kk),k);
                    for g=1:length(kk)
                        xi_mat(g,:)=xi{kk(g)+opt.nwarmup}';
                    end
                    xi_cell{m,k}=xi_mat;
                end
            end
            post_prob_nexp_t(k)=post_prob_nexp_t(k)/den;
        end

        post_prob_nexp_u=zeros(opt.nexp_umax,1);
        ui_cell=cell(size(out,1),opt.nexp_umax);
        
        for k=1:opt.nexp_umax
            den=0;
            for m=1:size(out,1)
                opt=out{m}.opt;
                den=den+opt.nloop-opt.nwarmup;
                ui=out{m}.ui;
                kk=find(out{m}.nexp_u(opt.nwarmup+1:end,1)==k);
                post_prob_nexp_u(k)=post_prob_nexp_u(k)+length(kk);
                if ~isempty(kk) && k>1
                    ui_mat=zeros(length(kk),k);
                    for g=1:length(kk)
                        ui_mat(g,:)=ui{kk(g)+opt.nwarmup}';
                    end
                    ui_cell{m,k}=ui_mat;
                end
            end
            post_prob_nexp_u(k)=post_prob_nexp_u(k)/den;
        end
       
        figure
        bar(post_prob_nexp_t)
        title('Posterior Probability of Number of Time Segments') 

        for k=1:opt.nexp_tmax
            xi_mat=cell2mat(xi_cell(:,k));
            if ~isempty(xi_mat) && k>1
                figure
                hold
                title(['Time Partition Points for a Partition of ',int2str(k)])
                for g=1:k-1
                    plot(xi_mat(:,g))
                end   
                for g=1:k-1
                    figure
                    histogram(xi_mat(:,g),'FaceAlpha', 1,'Normalization','probability')
                    title(['Distribution of ',int2str(g), ' Time Partition Point for a mixture of ',int2str(k)])
                end
            end
        end

        figure
        bar(post_prob_nexp_u)
        title('Posterior Probability of Number of Covariate Segments') 

        for k=1:opt.nexp_umax
            ui_mat=cell2mat(ui_cell(:,k));
            if ~isempty(ui_mat) && k>1
                figure
                hold
                title(['Covariate Partition Points for a Partition of ',int2str(k)])
                for g=1:k-1
                    plot(ui_mat(:,g))                    
                end   
                for g=1:k-1
                    figure
                    histogram(ui_mat(:,g),'FaceAlpha', 1,'Normalization','probability')
                    title(['Distribution of ',int2str(g), ' Covariate Partition Point for a mixture of ',int2str(k)])
                end
            end
        end
        
    end
    
end

