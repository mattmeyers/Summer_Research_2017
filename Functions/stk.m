function x = stk(mat)
    a = mat.domain(1); b = mat.domain(end); n = b-a;
    x = mat(:,1);
    for i = 2:size(mat,2)
       x = join(x,newDomain(mat(:,i),[(i-1)*n,i*n])); 
    end
end