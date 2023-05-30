function [V,N] = UniformVector(N,M,seed)
    rng('shuffle');
    mu = zeros(1,M);
    sigma = eye(M,M);
    R = mvnrnd(mu,sigma,N);
    V = abs(R./sqrt(sum(R.^2,2)));
end