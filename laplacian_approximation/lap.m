function L = lap( F )
[F1,F2 ] = grad(F);
L = div(F1,F2 );