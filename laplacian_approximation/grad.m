function [fx,fy ] = grad(M)

fx = M([2:end end],:)-M;
fy = M(:,[2:end end])-M;
