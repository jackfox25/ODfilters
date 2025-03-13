function [Q, R] = householderQR(A)
    [m, n] = size(A);
    Q = eye(m); % Initialize Q as the identity matrix
    R = A;
    
    for k = 1:n
        % Extract the vector to be reflected
        x = R(k:m, k);
        
        % Compute the Householder vector
        e1 = zeros(length(x), 1);
        e1(1) = norm(x) * sign(x(1));
        v = x + e1;
        v = v / norm(v);
        
        % Compute Householder matrix
        Hk = eye(m-k+1) - 2 * (v * v');
        
        % Embed Hk into a full-sized identity matrix
        H = eye(m);
        H(k:m, k:m) = Hk;
        
        % Apply the transformation
        R = H * R;
        Q = Q * H'; % Accumulate Q as product of H'
    end
end