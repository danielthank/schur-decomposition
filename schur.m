clear;
A = [14 -6 -4 0; 6 2 0 -4; -2 2 4 0; -2 2 0 4];
U = zeros(4); % W_perp

for z = 1:4
    idx = z;
    for i = 1:4
        if idx == 5
            break
        end
        tmp = A(:,i);
        for j = 1:idx-1
            tmp = tmp - dot(A(:,i),U(:, j))/(U(:,j).' * U(:,j)) * U(:, j); % Gram-Schmitz
        end
        if norm(tmp) > 0.00001
            U(:, idx) = tmp;
            idx = idx + 1;
        end
    end
    [a, b, c] = eig(U(:,z:4) \ (A * U(:,z:4))); % A|W_perp * x = /lambda x
    U(:,z) = U(:,z:4) * a(:,1); % add x to W_perp
end
U = U  ./ sqrt(sum(conj(U) .* U)); % normalize
display(U); % basis
display(U^-1 * A * U); % upper triangular
