function v = symmat2vec(M)

% This function computes a vectorization of SPD matrices using Mandel notation.

if ndims(M) == 2
    N = size(M,1);
    v = [];
    v = diag(M);
    for n = 1:N-1
        v = [v; sqrt(2).*diag(M,n)];
    end
else
    [D, ~, N] = size(M);
    v = [];
    for n = 1:N
        vn = [];
        vn = diag(M(:,:,n));
        for d = 1:D-1
            vn = [vn; sqrt(2).*diag(M(:,:,n),d)];
        end
        v = [v vn];
    end
    
end
end
