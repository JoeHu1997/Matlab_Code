function df_val=dmyf(x, funcMat, varMat)
% 返回值為2*2矩陣,內容為數值
%df=[df1/x1, df1/x2;
%    df2/x1. df2/x2];
n_X = size(x,2);%變數的個數
df =cell(n_X, n_X);
for i=1:n_X
    for j=1:n_X
        df{i,j} = diff(funcMat(1, i), varMat(1, j));
    end
end

df_val=zeros(n_X, n_X);

for i=1:n_X
    for j=1:n_X
        tmp_Var = cell(1,n_X);
        tmp_X = cell(1,n_X);
        for k=1:n_X
            tmp_Var{k} = varMat(1,k);
            tmp_X{k} = x(1,k);
        end
        df_val(i,j) = subs(df{i,j}, tmp_Var, tmp_X);
    end
end
end % end dmyf