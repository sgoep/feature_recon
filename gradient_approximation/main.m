clc
close all
clear

Nx = 256;
x =linspace(-1,1,Nx+1);
% % create Phantom
f = zeros(Nx+1);
[X,Y]=ndgrid(x,x);
% f(X.^2+Y.^2<= 1/2) = 0.2;
f( (X-0.2).^2+(Y-0.1).^2<=0.15^2  ) = 1;
f( (X-0.1).^2+(Y+0.25).^2<=0.25^2 ) = 1;
f( (X+0.55).^2+(Y-0.35).^2<=0.21^2) = 1;
f = f - 0.3*imrotate(ell2d( Nx, 0.1, 0.15, 0.3, -0.1, 1, 0),-50,'crop');
f = f - 0.3*imrotate(ell2d( Nx, 0.1, 0.15, 0.2, 0.1, 1, 0),-100,'crop');

% numerical data and noise
theta = 0:5:179;
data = radon(f,theta);
sino = data + randn(size(data))* max(abs(data(:))) * 0.01;

sigma = 4.5;
regparam = 0.01;

% ell1-regularized Gradient
[Ix_ell1,Iy_ell1] = approxgradRadon(sino,theta,sigma,'ell1',regparam);
Gmag_ell1 = sqrt(Ix_ell1.^2+Iy_ell1.^2);

E_ell1 = mycanny(Ix_ell1,-Iy_ell1,[]);
E_ell2 = mycanny(imgaussfilt(Ix_ell1,1),imgaussfilt(-Iy_ell1,1),[]);

figure; 
subplot(221); imagesc(Gmag_ell1); colormap gray; title('ell1: Gmag');
subplot(223); imagesc(E_ell1); colormap gray; title('ell1: canny');
subplot(224); imagesc(E_ell2); colormap gray; title('ell1: canny + gaussian');



%% 2-step and direct FBP gradient
sigma = 4;
[Ix_fbp,Iy_fbp] = approxgradRadon(sino,theta,sigma,'fbp');
Gmag_fbp = sqrt(Ix_fbp.^2+Iy_fbp.^2);

I = iradon(sino,theta,Nx);

[Ix2step,Iy2step] = gradient(imgaussfilt(I,sigma));
Gmag2step = sqrt(Ix2step.^2+Iy2step.^2);

E_fbp = mycanny(Ix_fbp,-Iy_fbp,[]);
E_fbp2 = mycanny(imgaussfilt(Ix_fbp,1),imgaussfilt(-Iy_fbp,1),[]);
E2step = edge(I,'canny',[],sigma);

figure; 
subplot(231);imagesc(I); colormap gray;  title('FBP');
subplot(232);imagesc(Gmag2step); colormap gray;  title('2step: Gmag');
subplot(233);imagesc(Gmag_fbp); colormap gray; title('FBP: Gmag');
subplot(235);imagesc(E2step);colormap gray; title('2-step');
subplot(236);imagesc(E_fbp);colormap gray; title('FBP');
