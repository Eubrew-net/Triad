function snoon=solar_noon(date,long)
long_min=long;
diajul=diaj(date);
phi=2*pi*(diajul-1.5)/365; % angulo diario (fractional year)
eqtime=222.18*(0.000075+0.001868*cos(phi)-0.038077*sin(phi)...
               -0.014615*cos(2*phi)-0.040849*sin(2*phi));
snoon= 720-4*long_min-eqtime;