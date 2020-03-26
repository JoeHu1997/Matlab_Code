function f=myf(x,funcMat, varMat)
% 块ま计xㄢ计,func1*2才腹跑计痻皚,var1*2才腹跑计痻皚い跑计
% 1*2痻皚,ず甧计

n_X = size(x,2);%跑计计
f_Val = zeros(1,n_X);
for i=1:n_X
    tmp_Var = cell(1,n_X);
    tmp_X = cell(1,n_X);
    for j=1:n_X
        tmp_Var{j} = varMat(1,j);
        tmp_X{j} = x(1,j);
    end
    f_Val(i) = subs(funcMat(1, i), tmp_Var, tmp_X);
end
f=f_Val;
end % end myf