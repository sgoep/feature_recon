function F = ell2d( Nx, ax, ay, mx, my, X0, angle);

x = linspace( -X0, X0, Nx+1) ;
[ X,Y ] = ndgrid( x, x );
R = sqrt( (X - mx).^2 / ax^2 + (Y - my).^2 / ay^2  );

F = zeros(Nx+1);
F( R < 1) = 1;