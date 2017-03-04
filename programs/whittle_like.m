function[log_prop_spec_dens]=whittle_like(y,fhat,nseg)
%Function computes the log Whittle likelihood
    nfreq=floor(nseg/2);
    log_prop_spec_dens=0;

    if(mod(nseg,2)==1)
        for i=1:size(y,2)
            log_prop_spec_dens=log_prop_spec_dens-sum(fhat(2:nfreq+1)+exp(y(2:nfreq+1,i)-fhat(2:nfreq+1)))-...
                0.5*(fhat(1)+exp(y(1,i)-fhat(1)))-...
                0.5*nfreq*log(2*pi);
        end
    else
        for i=1:size(y,2)
            log_prop_spec_dens=log_prop_spec_dens-sum(fhat(2:nfreq)+exp(y(2:nfreq,i)-fhat(2:nfreq)))-...
                0.5*(fhat(1)+exp(y(1,i)-fhat(1)))-...
                0.5*(fhat(nfreq+1)+exp(y(nfreq+1,i)-fhat(nfreq+1)))-...
                0.5*nfreq*log(2*pi);
        end
    end
end