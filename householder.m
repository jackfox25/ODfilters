% =========================================================================
% Function: householder
% Author:   Jack Fox
% Date:     03/10/25
%     Householder algorithm. Implementing notation from Born et al.
% =========================================================================
% Inputs:
%       A   (m+nxn+1) Augmented a priori information matrix
%       k   (1x1) Size of measurement vector
% Outputs:
%      hh   (m+nxn+1)       
% =========================================================================
function hh = householder(A,k)

    n = size(A,2)-1;
    m = size(A,1)-n;

    for k=1:n

        sumterm = 0;
        for i=k:m+n
            sumterm = sumterm + A(i,k)^2;
        end
        
        sigma = sign(A(k,k)) * sqrt(sumterm);
        uk = A(k,k) + sigma;
        A(k,k) = -sigma;
        ui = A(k+1:m+n,k);
        beta = 1/(sigma*uk);

        for j=k+1:n+1
            gamma = beta * 
        end

    end

end