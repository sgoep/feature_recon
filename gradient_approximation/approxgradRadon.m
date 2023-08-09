function [Ix, Iy] = approxgradRadon(sino,theta,sigma,method,varargin)

nn = size(sino,1);
filterExtent = ceil(4*sigma);
pp = -filterExtent:filterExtent;

g = @(x,k) (-1/(sigma.^3*sqrt(2*pi))).* x.*exp((-x.^2)/(2*sigma.^2));
G = g(pp,sigma);

% Normalize to ensure kernel sums to zero
negVals = G < 0;
posVals = G > 0;
G(posVals) = G(posVals)/sum(G(posVals));
G(negVals) = G(negVals)/abs(sum(G(negVals)));

dy1 = zeros(size(sino));
dy2 = zeros(size(sino));

% Calculate data corresponding to approximate partial derivatives
for k=1:numel(theta)
    tmp = conv(sino(:,k),G(:),'same');
    dy1(:,k) = cosd(theta(k))*tmp;
    dy2(:,k) = sind(theta(k))*tmp;
end

imgSize = 2*floor( size(sino,1)/(2*sqrt(2)) ) - 2;

if strcmpi(method,'fbp')
    disp('FBP');
    Ix = iradon(dy1,theta,imgSize);
    Iy = iradon(dy2,theta,imgSize);
    
elseif strcmpi(method,'ell1')
    disp('ell1');
    if nargin~=5
        error('regularization parameter must be given');
    end
    
    regparam = varargin{1}
    
    A  = @(x) radon(x,theta,nn);
    At = @(x) iradon(x,theta,imgSize,'None');

    Ix  = ista(A,At,dy1,regparam,true,false);
    Iy  = ista(A,At,dy2,regparam,true,false);
    
else
    error('Method not known');
end

