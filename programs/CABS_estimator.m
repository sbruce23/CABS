function [spec_est]=CABS_estimator(out)
%Function produces CABS estimator by averaging over posterior spectrum distribution

    if size(out,1)==1 %only one chain
        %pull in data from matfile for chain
        xi=out.xi;
        ui=out.ui;
        nseg_t=out.nseg_t;
        nexp_t=out.nexp_t;
        nexp_u=out.nexp_u;
        log_spec_hat=out.log_spec_hat;
        x=out.x;
        u=out.u;
        opt=out.opt;

        %initialize matrix
        spec_est=zeros(51,size(x,1),size(x,2));

        %average over posterior distribution
        for k=1:size(x,2)
            for p=(opt.nwarmup+1):size(log_spec_hat,1)
                xi_curr=xi{p};
                nseg_t_curr=nseg_t{p};       
                log_spec_hat_curr=log_spec_hat{p}(:,:,find(u(1,k)<=ui{p},1));
                for g=1:nexp_t(p)
                    if g==1
                        spec_est(:,1:xi_curr(g),k)=...
                            spec_est(:,1:xi_curr(g),k)+...
                            repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g))/...
                            (opt.nloop-opt.nwarmup);
                    else
                        spec_est(:,(xi_curr(g-1)+1):xi_curr(g),k)=...
                            spec_est(:,(xi_curr(g-1)+1):xi_curr(g),k)+...
                            repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g))/...
                            (opt.nloop-opt.nwarmup);
                    end
                end
            end
        end

    elseif size(out,1) > 1 %multiple chains

        %pull in data from matfile for chain
        x=out{1}.x;
        u=out{1}.u;
        
        %initialize matrix
        spec_est=zeros(51,size(x,1),size(x,2));
          
        %average over posterior distribution across chains
        den=0;
        for m=1:size(out,1)
            xi=out{m}.xi;
            ui=out{m}.ui;
            nseg_t=out{m}.nseg_t;
            nexp_t=out{m}.nexp_t;
            nexp_u=out{m}.nexp_u;
            log_spec_hat=out{m}.log_spec_hat;
            opt=out{m}.opt;
            den=den+opt.nloop-opt.nwarmup;
            for k=1:size(x,2)
                for p=(opt.nwarmup+1):size(log_spec_hat,1)
                    xi_curr=xi{p};
                    nseg_t_curr=nseg_t{p};       
                    log_spec_hat_curr=log_spec_hat{p}(:,:,find(u(k)<=ui{p},1));
                    for g=1:nexp_t(p)
                        if g==1
                            spec_est(:,1:xi_curr(g),k)=...
                                spec_est(:,1:xi_curr(g),k)+...
                                repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));
                        else
                            spec_est(:,(xi_curr(g-1)+1):xi_curr(g),k)=...
                                spec_est(:,(xi_curr(g-1)+1):xi_curr(g),k)+...
                                repmat(log_spec_hat_curr(:,g),1,nseg_t_curr(g));
                        end
                    end
                end
            end
        end       
        spec_est=spec_est./den;
    end
end

