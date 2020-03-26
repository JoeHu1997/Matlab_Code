function L3 = Lambda3(n, m, k1, k2, sigma01, sigma02, g, h, Om);

Knm = Kappa3(n, k1, m, k2);
K11 = Kappa3(1, k1, 1, k2);
Snm = Sigma3(n, sigma01, m, sigma02);
S11 = Sigma3(1, sigma01, 1, sigma02);
Gnm = Gamma3(n, m, k1, k2, sigma01, sigma02, g, h, Om);
G11 = Gamma2(k1, k2, sigma01, sigma02, g ,h, Om);
L11 = Lambda2(k1, k2, sigma01, sigma02, g, h, Om);
if n>m
    G2 = Gamma2(k1, k1, sigma01, sigma01, g, h, Om);
    L2 = Lambda2(k1, k1, sigma01, sigma01, g, h, Om);
elseif n<m %Only for waves that propagates +x
    G2 = Gamma2(abs(k2), abs(k2), sigma02, sigma02, g, h, Om);
    L2 = Lambda2(abs(k2), abs(k2), sigma02, sigma02, g, h, Om);
else
    G2 = 0;
    L2 = 0;
end

L3 = (1/Snm)*(Gnm*(Knm*sinh(Knm*h))...
    +G11*K11*(Knm*cosh(K11*h))...
    +G2*(Knm-K11)*(Knm*cosh(2*(Knm-K11)*h))...
    +L11*(Snm-S11)*(Knm*coth((Knm-K11)*h))...
    +L2*(2*S11-Snm)*(0.5*Knm*coth((2*K11-Knm)*h))...
    +Om*Knm*(L11+0.5*L2)...
    +((1/4)*Knm*(2*(Knm-K11)*(Snm-S11)+(2*K11-Knm)*(2*S11-Snm))));
end