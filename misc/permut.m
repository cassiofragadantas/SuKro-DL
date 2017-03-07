function P = permut(A,B)
% Construct a permutation matrix
P = zeros(A*B,B*A);
for a=1:A
    vec_uniA = zeros(A,1);
    vec_uniA(a)=1;
    for b=1:B
        vec_uniB = zeros(B,1);
        vec_uniB(b)=1;
        P_ab = kron(vec_uniA*vec_uniB',vec_uniB*vec_uniA');
        P = P_ab + P;
    end
end
% X = [111 211 311 411 511 611; 
%      112 212 312 412 512 612;
%      113 213 313 413 513 613;
%      121 221 321 421 521 621;
%      122 222 322 422 522 622;
%      123 223 323 423 523 623]
% 
% % Xmod = Prow*X
% 
% Xmod = X'*Prow

