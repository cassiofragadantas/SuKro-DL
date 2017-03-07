function [D] =  reord_inv(reordD,size_D,idx)
Ni = sqrt(size_D(1));
Mi = sqrt(size_D(2));

if nargin > 2
    D = reshape(reordD(idx),Ni*Ni,Mi*Mi);
else
    vec_D = reshape(reordD,Mi*Mi*Ni*Ni,1);
    Pnm = permut(Ni,Mi);
    P_rank = kron(eye(Mi),kron(Pnm,eye(Ni)));
    vec_reordD = P_rank.'*vec_D;
    D = reshape(vec_reordD,Ni*Ni,Mi*Mi);
end
end