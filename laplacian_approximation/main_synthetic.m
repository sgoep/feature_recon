clear
clc
close all

addpath('functions')

% Specify parameters for regularization
noise       = 0.01;
filter_type = 'ramlak'; % 'ramlak', 'quadratic' or 'gaussian'
lambda      = 0.001;    % regularization parameter ell_1 norm
mu          = 0.0;      % regularization parameter H_1 norm, if mu=0, only ell1 regularization is used
niter       = 100;     % number of iterations

%% Feature reconstruction for synthetic phantom

% Setting discretization
Nx = 200; 
Nal= 40; 
Ns = 300; 
Nt = 300;

x =linspace(-1,1,Nx+1);
al= linspace(0,pi*(1-1/Nal),Nal);
s = linspace(-1.5,1.5,Ns+1);
t = linspace(-1.5,1.5,Nt+1);

dx = x(2)-x(1);
ds = s(2)-s(1);

% create Phantom
f = zeros(Nx+1);
[X,Y]=ndgrid(x,x);
f( (X-0.2).^2+(Y-0.1).^2<=0.15^2  ) = 1;
f( (X-0.1).^2+(Y+0.25).^2<=0.25^2 ) = 1;
f( (X+0.55).^2+(Y-0.35).^2<=0.21^2) = 1;
f = f - 0.3*imrotate(ell2d( Nx, 0.1, 0.15, 0.3, -0.1, 1, 0),-50,'crop');
f = f - 0.3*imrotate(ell2d( Nx, 0.1, 0.15, 0.2, 0.1, 1, 0),-100,'crop');

% numerical data and noise
b = pi/ds;
data = my_radon(f,Nal,Ns,Nt,x,al,s,t);
rng('default');
rng(1);
data = data + randn(size(data))* max(abs(data(:))) * noise;

% Standard filtering and reconstruction
[filtered_data, filter] = discrete_filtering(data,'ramlak_classical',b,s);
backprojection          = my_backprojection(filtered_data,x,al,s);

% sparsify data accordingly
[feature_filtered_data, feature_filter] = discrete_filtering(data,filter_type,b,s);

feature_filtered_data = feature_filtered_data*ds^2;
laplacian = 4*del2(f);

% initiate FISTA
x0 = zeros(size(f));
L = 2*pi;

% Start l1 recon
fprintf('Starting Fista with Lambda=%f, Mu=%f, L=%f and %d Iterations \n', lambda, mu, L, niter)
xout = fista_H1(x0,L,feature_filtered_data,laplacian,lambda,mu,niter,Nx,Nal,Ns,Nt,x,al,s,t);


%% Thresholding

thresh_log = 0.005;
thresh_l1 = 0.05;
sigma = 1.3;%0.75;


fsize          = ceil(sigma*3) * 2 + 1;  % choose an odd fsize > 6*sigma;
op             = fspecial('log',fsize,sigma);
log_fbp        = imfilter(backprojection, op, 'replicate');    
log_fbp_thresh = laplacian_of_gaussian_thresholding(log_fbp, thresh_log);
xout_thresh    = laplacian_of_gaussian_thresholding(xout, thresh_l1);

figure
colormap gray
subplot(221), imagesc(log_fbp), title('LoG F'), colorbar
subplot(222), imagesc(log_fbp_thresh), title('Thresh + LoG F'), colorbar
subplot(223), imagesc(xout), title('xout'), colorbar
subplot(224), imagesc(xout_thresh), title('Feature Extraction + Thresh'), colorbar
