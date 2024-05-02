function random_Laplacian(n,potential=true)
    Q = 10*rand()*rand(n,n); D = -rand(n)
    for i=1:n
     Q[i,i] = 0; Q[i,i] = -sum(Q[i,:])
    end
    # -L = Q + D
    if potential==true
     return -(Q + diagm(D))
    end
    return -Q
 end