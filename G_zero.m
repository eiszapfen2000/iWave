function G_0 = G_zero(n, deltaQ, sigma)
%%
if nargin < 3, sigma = 1; end
q = deltaQ:deltaQ:n*deltaQ;
G_0 = sum(q.*q.*exp(-sigma.*q.*q));
end