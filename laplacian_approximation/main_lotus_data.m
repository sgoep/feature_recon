clear
clc
%close all

addpath('functions')

% Specify parameters for regularization
filter_type = 'ramlak'; % 'ramlak', 'quadratic' or 'gaussian'
lambda      = 0.01;    % regularization parameter ell_1 norm
mu          = 0.0;      % regularization parameter H_1 norm, if mu=0, only ell1 regularization is used
niter       = 10;     % number of iterations


%%
load('fips_lotus_parallelbeam.mat')
data = sino(1:3:end,1:5:end);
thetas = theta(1:5:end);
sloc = p(1:3:end);

n = -(length(sloc)-1)/2:(length(sloc)-1)/2;
stdev = 8;
y = exp(-1/2*(n/stdev).^2);

Nx = 590-1;
Nal = length(thetas); 
Ns = length(sloc)-1; 
Nt = Ns; 

% Take x in -1.2 to 1.2 otherwise Lotus will not be visible completely
x = linspace(-1.2,1.2,Nx+1);
al = linspace(0,pi*(1-1/Nal),Nal);
s = linspace(-1.5,1.5,Ns+1);
t = s;

dx = x(2)-x(1);
ds = s(2)-s(1);

b = pi/ds;

% Standard filtering and reconstruction
[filtered_data, filter] = discrete_filtering(data,'ramlak_classical',b,s);
backprojection          = my_backprojection(filtered_data,x,al,s);

% sparsify data accordingly
[feature_filtered_data, feature_filter] = discrete_filtering(data,filter_type,b,s);
if ~strcmp(filter_type, 'gaussian')
    feature_filtered_data = feature_filtered_data*ds^2;
end

% initiate FISTA
x0 = zeros(size(backprojection));
L = 2*pi;

% Start l1 recon
fprintf('Starting Fista with Lambda=%f, Mu=%f, L=%f and %d Iterations \n', lambda, mu, L, niter)
xout = fista_H1(x0,L,feature_filtered_data,[],lambda,mu,niter,Nx,Nal,Ns,Nt,x,al,s,t);


%% Thresholding

thresh_log = 0.018;
thresh_l1  = thresh_log; %0.005;
sigma      = 1; %0.75;

fsize          = ceil(sigma*3) * 2 + 1;  % choose an odd fsize > 6*sigma;
op             = fspecial('log',fsize,sigma);
log_fbp        = imfilter(backprojection,op,'replicate');    
log_fbp_trhesh = laplacian_of_gaussian_thresholding(log_fbp, thresh_log);

xout_filt        = imfilter(xout, op, 'replicate');
xout_filt_thresh = laplacian_of_gaussian_thresholding(xoutfilt, thresh_l1);

figure
colormap gray
subplot(221), imagesc(log_fbp), title('LoG F'), colorbar
subplot(222), imagesc(log_fbp_trhesh), title('Thresh + LoG F'), colorbar
subplot(223), imagesc(xout_filt), title('xout'), colorbar
subplot(224), imagesc(xout_filt_thresh), title('Feature Extraction + Thresh'), colorbar
