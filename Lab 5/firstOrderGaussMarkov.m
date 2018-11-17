function [fogm] = firstOrderGaussMarkov(t_correlation, t_sample, s)
    % computes the first order guass-markov process given: 
    % correlation time, sample time, and random sequence
    Z = zeros(length(s),1);
    Z(1) = s(1);
    Beta = 1/t_correlation;
    for i=2:length(s)
        Z(i) = exp(-Beta*t_sample)*Z(i-1)+s(i);
    end
    fogm = Z;
end

