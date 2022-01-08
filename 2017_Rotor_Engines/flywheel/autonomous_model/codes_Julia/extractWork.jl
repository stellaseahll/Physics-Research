function extractWork(rho,H)
    (r,v) = eig(rho)
    (E,v) = eig(H)
    r = sort(real(r))
    r = r[end:-1:1]
    E = real(E)
    # println(r)
    # println(E)
    tmp = H.*rho
    S = real(sum(tmp[:])) - sum(r.*E)
    return S
end