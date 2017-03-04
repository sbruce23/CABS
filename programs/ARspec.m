function [spect]=ARspec(phi,sigma,freq)
%Function to produce true power spectrum for AR process
%   Function takes in AR coefficients, covariance, and frequencies, and 
%   outputs the true power spectrum for given frequencies.

dim = size(phi);
len = dim(2);
phispect = zeros(2,2,length(freq));
spect = zeros(2,2,length(freq));

for k=1:length(freq)
    phispect(:,:,k) = eye(2);
    for j=1:(len/2)
        if j==1
            bigmat = phi(:,1:2).* exp(-2*pi*sqrt(-1)*freq(k));
        else
            bigmat = phi(:,(2*j-1):2*j) .* exp(-2*j*pi*sqrt(-1)*freq(k));
        end
        phispect(:,:,k) = phispect(:,:,k) - bigmat;
    end
    spect(:,:,k) = phispect(:,:,k)\sigma*conj(inv(phispect(:,:,k)));
end    
