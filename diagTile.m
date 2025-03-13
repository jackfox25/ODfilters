% =================================================================================================
% Function: diagTile
% Replicates a square matrix along the diagonal of a larger square matrix
% =================================================================================================
% Inputs:
%        A   square matrix to tile
%        n   number of copies along diagonal
% Outputs:
%    diagA   i.e. diagTile(A,5) = blkdiag(A,A,A,A,A)
% =================================================================================================
function diagA = diagTile(A,n)

    sa = size(A,1);

    diagA = zeros(n*size(A));
    for i=1:n
        diagA(sa*(i-1)+1:sa*i,sa*(i-1)+1:sa*i) = A;
    end

end