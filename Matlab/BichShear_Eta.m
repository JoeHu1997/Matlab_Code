clc;
clear all;

t=0;
xleft=0;
xright=100;
x=xleft:0.1:xright;
g=9.8066;
h=5;
a1=0.45;
a2=0.2;
f1=0.4;
f2=0.4;
k1=0.5096;
k2=0.5124;
k11=0.5215;
k12=0.3297;
Um=0.35;
Om=0.1;
stokes1=1;
stokes2=0;
stokes3=1;

Us=Um+0.5*Om*h;
sigma1=2*pi*f1;
sigma2=2*pi*f2;

sigma01=g*k1*tanh(k1*h)/(sigma1-Us*k1+Om*tanh(k1*h));
sigma02=g*k2*tanh(k2*h)/(sigma2-Us*k2+Om*tanh(k2*h));

sigma101=g*k11*tanh(k11*h)/(sigma1-Us*k11+Om*tanh(k11*h));
sigma102=g*k12*tanh(k12*h)/(sigma2-Us*k12+Om*tanh(k12*h));

[G10, G01, G20, G02, G11p, G11m, L20, L02, L11p, L11m,...
    G30, G03, G21p, G21m, G12p, G12m,...
    L30, L03, L21p, L21m, L12p, L12m, L130, L103]...
    = BSCproperties(a1, a2, Us, Om, k1, k2, sigma01, sigma02, g, h);


sol10=a1*cos(k1*x-sigma1*t);
sol01=a2*cos(k2*x-sigma2*t);

sol20=(a1^2/2)*L20*cos(2*k1*x-2*sigma1*t);
sol02=(a2^2/2)*L02*cos(2*k2*x-2*sigma2*t);
sol11p=(a1*a2)*L11p*cos((k1+k2)*x-(sigma1+sigma2)*t);
sol11m=(a1*a2)*L11m*cos((k1-k2)*x-(sigma1-sigma2)*t);

sol30=(a1^3/2)*L30*cos(3*k1*x-3*sigma1*t);
sol03=(a2^3/2)*L03*cos(3*k2*x-3*sigma2*t);
sol21p=(a1^2*a2/2)*L21p*cos((2*k1+k2)*x-(2*sigma1+sigma2)*t);
sol21m=(a1^2*a2/2)*L21m*cos((2*k1-k2)*x-(2*sigma1-sigma2)*t);
sol12p=(a1*a2^2/2)*L12p*cos((k1+2*k2)*x-(sigma1+2*sigma2)*t);
sol12m=(a1*a2^2/2)*L12m*cos((k1-2*k2)*x-(sigma1-2*sigma2)*t);
sol130=a1*L130*cos(k1*x-sigma1*t);
sol103=a2*L103*cos(k2*x-sigma2*t);

sol1=(sol10+sol01);
sol2=(sol20+sol02+sol11p+sol11m);
sol3=(sol30+sol03+sol21p+sol21m+sol12p+sol12m+sol130+sol103);

sol11=(sol110+sol101);
sol12=(sol120+sol102+sol111p+sol111m);

plot(x,stokes1*sol11,'--b','linewidth',1); %µù¡G¤@¶¥ªi®ö¸Ñ
  hold on;
  plot(x,stokes3*(sol1+sol2+sol3),'-b','linewidth',1); %µù¡G3¶¥ªi®ö¸Ñ
  hold on;
  grid on;
  axis([xleft xright -1.2 1.2]);
  xlabel('x (m)','fontSize',16);
  ylabel('surface Elevation (m)','fontSize',16);
  title('Surface Elevation','fontSize',16);
  legend('first-order','third-order');
  