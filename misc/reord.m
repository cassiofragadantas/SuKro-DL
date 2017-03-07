% This function actually gives the transpose of the rearrangement operator
% defined in [1], but it is equally valid to the proposed technique.
%
%  [1] C.F. Dantas, M. N. da Costa, R.R. Lopes, "Learning Dictionaries as a
%      sum of Kronecker products", To appear.

function [reordD] =  reord(D,idx)
[N,M] = size(D);
Ni = sqrt(N);
Mi = sqrt(M);

if nargin > 1
    reordD = reshape(D(idx),Mi*Ni,Mi*Ni);
else
    vec_D = reshape(D,Mi*Mi*Ni*Ni,1);
    Pnm = permut(Ni,Mi);
    P_rank = kron(eye(Mi),kron(Pnm,eye(Ni)));
    vec_reordD = P_rank*vec_D;
    reordD = reshape(vec_reordD,Mi*Ni,Mi*Ni);
end

end