function [ IntePart ] = Integral1( h1,h2,dmax,D )

S = sqrt(-(D-(h2*dmax)-(h1*dmax))*(D+(h2*dmax)-(h1*dmax))*(D-(h2*dmax)+(h1*dmax))*(D+(h2*dmax)+(h1*dmax)))-2*((h2*dmax)^2)*asec((2*D*(h2*dmax))/(D^2 + (h2*dmax)^2 -(h1*dmax)^2))-2*(((h1*dmax)^2)*asec((2*D*(h1*dmax))/(D^2 - (h2*dmax)^2 +(h1*dmax)^2)));
f = @(x)(2*((h2*dmax)^2-(h1*dmax)^2-D*sqrt(D^2+((h2*dmax)^2-D^2)*sec(x).^2)+D*cos(2*x).*(D-sqrt(D^2+((h2*dmax)^2-D^2)*sec(x).^2))))./(cos(x));
xitam = acos(((h1*dmax)^2+D^2-(h2*dmax)^2)/(2*h1*dmax*D));
IntePart = quadgk(f,0,xitam);
IntePart = IntePart/S;
