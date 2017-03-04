function [param] = setMCMCOptions(varargin)
%Function to set options for MCMC sampler
%   Options
%   seed: random number seed
%   nloop: number of total MCMC iterations 
%   nwarmup: number of iterations to use as burn-in
%	nexp_tmax: maximum number of time-based partition segments(if
%       Rev_Jump_t=0, this will be the exact number of time-
%       based partition segments since no birth and death 
%       moves are allowed)
%   nexp_umax: maximum number of covariate-based partition segments(if
%       Rev_Jump_u=0, this will be the exact number of covariate-
%       based partition segments since no birth and death 
%       moves are allowed)
%   Rev_Jump_t: 1 if want to have reversible jumps in time partition
%       and 0 otherwise (i.e. 1 enables birth and death moves
%       while 0 limits you only to within model moves) 
%   Rev_Jump_u: 1 if want to have reversible jumps in covariate
%       partition and 0 otherwise (i.e. 1 enables birth and 
%       death moves while 0 limits you only to within model moves) 
%   tmin: minimum number of time observations in each time segment
%   umin: minimum number of unique covariate values in each covariate
%       segment
%   prob_mm1: probability of a small move (+/- 1 value) in
%       time-based within move step
%   prob_mm1_u: probability of a small move (+/- 1 value) in
%       covariate-based within move step
%   nbasis: number of basis functions to use for spectral estimation
%   sigmasqalpha: prior known variance for constant term in spectral
%       estimation
%   tau_prior_a: tau prior distribution parameter 
%   tau_prior_b: tau prior distribution parameter
%   var_inflate: inflation factor for variance of beta distribution
%   conv_diag: calculate and output convergence diagnostics or not
%   options: settings for optimization subroutine used to draw from
%       distribution of beta    
%   fname: name of output file saved containing all the results
%   	(include path if you wish to save the output file to a
%   	different folder than the current directory)
%   batchsize: amount of iterations to hold in memory between outputs
%   	to .mat file.  Increasing this parameter will use more
%   	RAM, but require fewer dumps to .mat file thus
%   	decreasing run time.
%   parcomp: indicator for parallel processing (if parcomp='on', use parallel
%       chains to obtain the desired number of iterations; if parcomp='off',
%       used one chain to obtain desired number of iterations).  Requires
%       parallel computing toolbox
%   disp_plots: turn on/off display of summary plots for sampler
%       diagnostics
%   disp_iter: display every nth iteration in console while sampler running
%       (e.g. 0 means don't display any iterations, 1 means display every
%       iteration, 5 means display every 5th iteration, etc.)
        param = struct('seed',unidrnd(1000000),'nloop',10000,'nwarmup',2000,...
                       'nexp_tmax',10,'nexp_umax',10,...
                       'Rev_Jump_t', 1,'Rev_Jump_u', 1,...
                       'tmin',40,'umin',1,...
                       'prob_mm1', 0.8,'prob_mm1_u',0.8,...
                       'nbasis', 7,'sigmasqalpha', 100,...
                       'tau_prior_a', -1,'tau_prior_b', 0,...
                       'var_inflate', 1.0,'conv_diag','off',...
                       'options',optimset('Display','off','GradObj','on',...
                       'Hessian','on','MaxIter',10000,'MaxFunEvals',10000,...
                       'TolFun',0.001,'TolX',0.00001),...
                       'fname','CABSout.mat',...
                       'batchsize',100, 'parcomp','off',...
                       'disp_plots','off', 'disp_iter',1);
               
        if nargin > 0
            op_names = showOptionNames();
            paramLen = length(struct2cell(param)); 
            if ischar(varargin{1})  %match and set arguments by names
                if nargin > 2*paramLen
                      error('Too many input arguments.');     
                elseif mod(nargin, 2) == 0
                      for i = 1:2:nargin
                         id = find(strcmp(op_names,varargin{i}));
                         if isempty(id)
                          error('%s is an invalid name. Type showOptionNames() for more details.',varargin{i});
                         else                            
		                   param = setVal(param, op_names{id}, varargin{i+1});   
                         end
                      end
                else
		           error('Either option name or its corresponding value is missing.');
                end
	        else  %set arguments directly
               if nargin <= paramLen
	               
                    for i = 1:nargin
	                  param = setVal(param, op_names{i}, varargin{i});
                    end
               else
	               error('Too many input arguments.');
               end      
           end

        end
        
%         %check to make sure nloop multiple of batchsize
%         if parcomp==0 && mod(param.nloop,param.batchsize)~=0
%             warning(['Number of iterations not a multiple of batch size. Last ',...
%                 num2str(mod(param.nloop,param.batchsize)),...
%                 ' iterations will not be output to .mat file.']);
%         elseif parcomp==1
%             %%%determine how many workers available for computing and
%             %how many iterations per worker
%             cluster = parcluster('local');
%             nloop=zeros(cluster.NumWorkers,1);
%             for i=1:cluster.NumWorkers
%               if i < cluster.NumWorkers
%                 nloop(i)=floor((param.nloop-param.nwarmup)/cluster.NumWorkers);
%               elseif i == cluster.NumWorkers
%                 nloop(i)=(param.nloop-param.nwarmup) - sum(nloop(1:(i-1)));
%               end
%             end
%             nloop=nloop+param.nwarmup;
%             
%             if sum(mod(nloop,param.batchsize)~=0)>0
%                 
%         end

end
