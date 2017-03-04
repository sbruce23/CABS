function [xx]= lin_basis_func(freq,nbeta)
%Function calculates basis function values for specified frequencies
n=length(freq);
xx=ones(n,nbeta);
for j=2:nbeta
    xx(:,j)=sqrt(2)*cos((j-1)*pi*freq)/(pi*(j-1));
end