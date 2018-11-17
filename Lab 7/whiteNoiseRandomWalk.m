function [s,rw] = whiteNoiseRandomWalk(n,seed)
rng(seed);
s = randn(n,1);
rw = cumsum(s);
end