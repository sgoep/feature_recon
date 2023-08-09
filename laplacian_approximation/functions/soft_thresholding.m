function y = SoftThresh(x,alpha)
    y = sign(x).*(max(0,abs(x)-alpha));
end