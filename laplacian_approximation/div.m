function fd = div(Px,Py )

% div - divergence operator
%
%   Note that the -div and grad operator are adjoint
%   of each other such that 
%       <grad(f),g>=<f,-div(g)>
%
%   See also: grad.
%

fx = Px-Px([1 1:end-1],:);         
fx(1,:)   = Px(1,:);        % boundary
fx(end,:) = -Px(end-1,:);        

fy = Py-Py(:,[1 1:end-1]);
fy(:,1)   = Py(:,1);    % boundary
fy(:,end) = -Py(:,end-1);
        
fd = fx+fy;
