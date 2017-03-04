function [spec_est, out] = CABS(x, u, opt)
%Wrapper function to determine parallel resources available and call sampler and produce estimator 

    %%%set seed
    rng(opt.seed);

    if strcmp(opt.parcomp,'on')

        %%%determine how many workers available for computing and
        %how many iterations per worker
        cluster = parcluster('local');
        nloop=zeros(cluster.NumWorkers,1);
        for i=1:cluster.NumWorkers
          if i < cluster.NumWorkers
            nloop(i)=floor((opt.nloop-opt.nwarmup)/cluster.NumWorkers);
          elseif i == cluster.NumWorkers
            nloop(i)=(opt.nloop-opt.nwarmup) - sum(nloop(1:(i-1)));
          end
        end

        %%%use parfor to start parallel processing
        out=cell(cluster.NumWorkers,1);
        parfor i=1:cluster.NumWorkers
            paropt=setVal(opt,'nloop',nloop(i)+opt.nwarmup);
            paropt=setVal(paropt,'fname',[opt.fname,'_',num2str(i)]);
            paropt=setVal(paropt,'seed',unidrnd(1000000));

            tic         
            out{i}=CABS_sampler(x,u,paropt);
            toc
        end

        %%%shut down parallel pool
        p = gcp;
        delete(p)

    else
        out=CABS_sampler(x,u,opt);
    end
    
    %%%produce CABS spectral estimator
    spec_est=CABS_estimator(out);
    
    %%%produce CABS summary plots
    if strcmp(opt.disp_plots,'on')
        CABS_plot(out);
    end

end