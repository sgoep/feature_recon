function x = ista(A,At,y,lambda,varargin)
% function x = ista(A,At,y,lambda,varargin)
%
% INPUT:
% A         (function handle) forward operator
% At        (function handle) adjoint/transpose of A
% y         (1d array) data
% lambda    (real number) regularizaiton parameter
% x0        (array) initial guess
% showlog   (boolean) flag to print interation information
% dispiter  (boolean) flag; if true, current iterate is displayed

maxit = 50;

% Current iterates and gradients
x = At(y);
xOld = zeros(size(x));
gOld = dminfun(A,At,xOld,y);
g = dminfun(A,At,x,y);

showlog = false;

if nargin==5
    showlog = varargin{1};
elseif nargin==6
    showlog  = varargin{1};
    dispiter = varargin{2};
    if dispiter
        ff = figure;
    end
end

SoftThresh = @(x,lambda) (sign(x).* max(abs(x)-lambda,0));
figure, colormap gray
for it=1:maxit
    % Barzilai Borwein Stepsize
    dX = x - xOld;
    dG = g - gOld;
    
    stepBB = dX(:)'*dG(:);
    if abs(stepBB)<eps
        return;
    end
    
    stepBB = (dX(:)'*dX(:))/stepBB;

    gOld = g;
    xOld = x;
    
    % Gradient Step + Soft thresholding
    x = SoftThresh(x - stepBB*g, lambda*stepBB);
    
    if mod(it,1)==0
       subplot(121), imagesc(x)
       drawnow
    end

    % Gradient of the least squares error evaluatet at the current iterate x
    g = dminfun(A,At,x,y);
    
    if showlog
       fprintf(1,'it = %d; \t res = %f \t fval = %f \t step = %f \n',it,norm(A(x)-y),minfun(A,x,y),stepBB);
    end
    if dispiter
       figure(ff); imagesc(abs(x)); colormap gray; drawnow;
       pause;
    end
end

%%%%
function val = minfun(A,x,y)
res = A(x)-y;
val = 0.5*(res(:)'*res(:));

function grad = dminfun(A,At,x,y)
grad = At( A(x) - y );
% grad = grad./norm(grad(:));