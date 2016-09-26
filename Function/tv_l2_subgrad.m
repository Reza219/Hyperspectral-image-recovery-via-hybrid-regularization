function G = tv_l2_subgrad(X)
% Returns a subgradient G at X

[m,n] = size(X);

% The intensity differences across rows
row_diff = X(1:m-1,2:n)-X(1:m-1,1:n-1);
% The intensity differences across columns
col_diff = X(2:m,1:n-1)-X(1:m-1,1:n-1);

% Compute the total variation
grad_norms = sqrt(row_diff.^2+col_diff.^2);

% Find a subgradient
G = zeros(m,n);
% When non-differentiable, set to 0
grad_norms(grad_norms==0)=inf;
G(1:m-1,1:n-1) = G(1:m-1,1:n-1) - row_diff./grad_norms;
G(1:m-1,2:n)   = G(1:m-1,2:n)   + row_diff./grad_norms;
G(1:m-1,1:n-1) = G(1:m-1,1:n-1) - col_diff./grad_norms;
G(2:m,1:n-1)   = G(2:m,1:n-1)   + col_diff./grad_norms;

end