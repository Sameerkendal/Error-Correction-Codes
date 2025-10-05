function [G_sys,G]=get_systematic_generator_matrix(k, n, m)

G = gf(zeros(k, n), m);

for ri=1:k
    for ci=1:n
        G(ri, ci)=gf(ci-1, m)^(ri-1);
    end
end

G_sys = G(1:k,1:k)\G;
end