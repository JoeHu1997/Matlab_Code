function f=myf(x,funcMat, varMat)
% ��J�޼�x����Ӽƭ�,func��1*2�Ÿ��ܼƯx�},var��1*2�Ÿ��ܼƯx�}�����ܼ�
% ��^�Ȭ�1*2�x�},���e���ƭ�

n_X = size(x,2);%�ܼƪ��Ӽ�
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