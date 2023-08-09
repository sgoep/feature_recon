function [TG,filt] = DiscreteFilt(G,filterMode,b,s)

[m,n] = size(G);
TG = zeros(m,n);

Ns = length(s)-1;
Nh = ceil(Ns/2);
ell = -Nh:Nh;

even = find(mod(ell,2)==0);
odd = find(mod(ell,2)==1);
zero = find(ell==0);

if strcmp(filterMode,'quadratictype')
    filt(even) = -4./(pi^2*ell(even).^2);
    filt(odd) = 4./(pi^2*ell(odd).^2);
    filt(zero) = -2/3;
    filt = filt*b^3;
    
      
elseif strcmp(filterMode,'ramlaktype')
    filt(even)=2./(pi^2*ell(even).^2);
    filt(odd)=-2./(pi^2*ell(odd).^2)+24./(pi^4*ell(odd).^4);
    filt(zero) = -1/6;
    filt = filt*b^3;
    
elseif strcmp(filterMode,'gaussian')
    b = 1/2;
    a = 12;
    ll = ell;
    filt = (2*sqrt(pi)/a^3)*exp(-ll.^2/a^2).*(2*ll.^2-a^2);
    
elseif strcmp(filterMode,'ramlaknatterer')
    filt(even)=0;
    filt(odd)=-1./(pi^2*ell(odd).^2);
    filt(zero)=1/4;
    filt = filt*b^2/(2*pi^2);
  
end


h =  1/(2*b); % Peroid h
for i = 1:n
    TG(:,i) = h*conv(G(:,i),filt,'same');
end
end

