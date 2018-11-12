function g = G(P, n, deltaQ, sigma)
%%
if nargin < 4, sigma = 1; end

G0 = G_zero(n, deltaQ, sigma);
kernelSize = 2 * P + 1;
g = zeros(kernelSize);

q = deltaQ:deltaQ:n*deltaQ;
qSquare = q.*q;

for l=-P:1:P
    for k=-P:1:P
        r = sqrt(l.^2 + k.^2);
        b = besselj(0, r.*q);
        powers = exp(-sigma * qSquare);
        scalars = qSquare .* powers .* b;
        g(l + P + 1, k + P + 1) = sum(scalars) / G0;
    end
end
end