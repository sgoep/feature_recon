function e = laplacian_of_gaussian_thresholding(b, thresh)

    [m, n] = size(b);
    rr = 2:m-1; cc=2:n-1;
    e = false(m,n);
    
    [rx,cx] = find( b(rr,cc) < 0 & b(rr,cc+1) > 0 ...
        & abs( b(rr,cc)-b(rr,cc+1) ) > thresh );   % [- +]
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr,cc-1) > 0 & b(rr,cc) < 0 ...
        & abs( b(rr,cc-1)-b(rr,cc) ) > thresh );   % [+ -]
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc) > 0 ...
        & abs( b(rr,cc)-b(rr+1,cc) ) > thresh);   % [- +]'
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr-1,cc) > 0 & b(rr,cc) < 0 ...
        & abs( b(rr-1,cc)-b(rr,cc) ) > thresh);   % [+ -]'
    e((rx+1) + cx*m) = 1;
    [rz,cz] = find( b(rr,cc)==0 );
    if ~isempty(rz)
        zero = (rz+1) + cz*m; 
        zz = (b(zero-1) < 0 & b(zero+1) > 0 ...
            & abs( b(zero-1)-b(zero+1) ) > 2*thresh);     % [- 0 +]'
        e(zero(zz)) = 1;
        zz = (b(zero-1) > 0 & b(zero+1) < 0 ...
            & abs( b(zero-1)-b(zero+1) ) > 2*thresh);     % [+ 0 -]'
        e(zero(zz)) = 1;
        zz = (b(zero-m) < 0 & b(zero+m) > 0 ...
            & abs( b(zero-m)-b(zero+m) ) > 2*thresh);     % [- 0 +]
        e(zero(zz)) = 1;
        zz = (b(zero-m) > 0 & b(zero+m) < 0 ...
            & abs( b(zero-m)-b(zero+m) ) > 2*thresh);     % [+ 0 -]
        e(zero(zz)) = 1;
    end
end