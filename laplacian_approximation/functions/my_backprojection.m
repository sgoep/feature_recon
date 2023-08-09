function bf = my_backprojection(g,x,al,s)

sum = zeros(length(x));
[X,Y] = ndgrid(x,x);
theta=[cos(al);sin(al)]';

for i = 1:length(al)
    sp = X.*theta(i,1)+ Y.*theta(i,2);
    proj = g(:,i);
    aux = interp1(s,proj,sp,'linear',0);
    sum = sum + aux;
end
bf = pi/length(al).*sum;
