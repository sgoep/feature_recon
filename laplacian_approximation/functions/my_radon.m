function Rf = my_radon(f,Nal,Ns,Nt,x,al,s,t)

dt = 3/length(t);

[X,Y] = meshgrid(x,x);
[S,T] = meshgrid(s,t);
Rf = zeros(Ns+1,Nal);

for i=1:Nal        
    Trot = S*cos(al(i))-T*sin(al(i));
    Srot = S*sin(al(i))+T*cos(al(i));
    frot = interp2(X,Y,f,Srot,Trot,'linear',0);
    fsum = sum(frot,1);
    Rf(:,i) = fsum*dt;
end

