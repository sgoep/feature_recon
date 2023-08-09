function xout = fista_H1(x0, L, TG, TF0, lambda, mu, niter,Nx,Nal,Ns,Nt,x,al,s,t)

% ########################################################################
% Solves xout = argmin {||Rf-g||_2^2 + ||f||_H^1^2 + lambda*||f||_1} for Radon Transform R,
% sparse f and g, using FISTA
% x0       - startingpoint
% L        - Lipschitz constant for R
% TG       - sparse data
% TF0      - \Delta F
% lambda   - regularization parameter
% mu       - tuning parameter
% niter    - number of Iterations
% Nx,...,s - variables for Radon transform and backprojection
% ########################################################################

xout = x0;
y = x0;
t_step = 1;

for k=1:niter
        
    tprev = t_step;
    xprev = xout;
      
    TGiter= my_radon(y,Nal,Ns,Nt,x,al,s,t);
    Fgrad = my_backprojection(TGiter-TG,x,al,s);
           
    [g1,g2] = grad(y);
    
    xout = soft_thresholding(y - (Fgrad - mu*div(g1,g2))/L, lambda/L);
    t_step = (1+sqrt(1+4*tprev^2))/2;
    y = xout+(tprev-1)/t_step*(xout-xprev);

    if mod(k,10)==0    
        fprintf('Iteration %d of %d completed. \n', k, niter)
    end
    
end

end