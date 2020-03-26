
% x0ㄢ跑计癬﹍,funcMat琌ㄢよ祘,varㄢよ祘ㄢ跑计,eps北弘
% 箉猭秆じ獶絬┦よ祘舱

a1=0.45;
a2=0.2;
sigma1=2*pi()*0.4;
sigma2=2*pi()*0.3;
Um=0.35;
Om=0.1;
h=5;
g=9.8066;

Us=Um+0.5*Om*h;

syms k1 k2
sigma01=g*k1*tanh(k1*h)/(sigma1-Us*k1+Om*tanh(k1*h));
sigma02=g*k2*tanh(k2*h)/(sigma2-Us*k2+Om*tanh(k2*h));
G20 = Gamma2(k1, k1, sigma01, sigma01, g, h, Om);
G02 = Gamma2(k2, k2, sigma02, sigma02, g, h, Om);
G11p = Gamma2(k1, k2, sigma01, sigma02, g, h, Om);
G11m = Gamma2(k1, -k2, sigma01, -sigma02, g, h, Om);
L20 = Lambda2(k1, k1, sigma01, sigma01, g, h, Om);
L02 = Lambda2(k2, k2, sigma02, sigma02, g, h, Om);
L11p = Lambda2(k1, k2, sigma01, sigma02, g, h, Om);
L11m = Lambda2(k1, -k2, sigma01, -sigma02, g, h, Om);

DR1=-sigma1+Us*k1+sigma01;
    
DR2=-sigma2+Us*k2+sigma02;

 x0 = [0.1 0.1];
 funcMat=[DR1 DR2];
      
 var=[sym('k1') sym('k2')];
 eps=1.0e-4;

n_Var = size(var,2);%跑计计
n_Func = size(funcMat,2);%ㄧΑ计
n_X = size(x0,2);%跑计计

if n_X ~= n_Var && n_X ~= n_Func
    fprintf('Expression Error!\n');
    exit(0);
end

myfd=myf(x0, funcMat, var);
dmyfd=dmyf(x0, funcMat, var);

r=x0-myfd*inv(dmyfd);
n=0;
tol=1;
while tol>=eps
    x0=r;
    myfd=myf(x0, funcMat, var);
    dmyfd=dmyf(x0, funcMat, var);
    r=x0-myfd*inv(dmyfd);
    tol=norm(r-x0);
    n=n+1;
    if(n>100000)
        disp('˙计びよ祘ぃΜ滥');
        return;
    end
end
[r,n]