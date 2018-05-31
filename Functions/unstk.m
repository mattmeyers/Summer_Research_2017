function x = unstk(vec ,dom)
    n = dom(2)-dom(1); 
    x = restrict(vec,dom);
    for i = 2:((vec.domain(end)-vec.domain(1))/n)
       x(:,i) = newDomain(restrict(vec,[(i-1)*n,i*n]),dom); 
    end
end