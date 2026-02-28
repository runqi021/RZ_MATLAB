function [phi_hat, phi_std, phi_CI] = phase_jackknife(C12_vec)
% phase_jackknife_from_C
%   Jackknife estimate of mean phase + std/CI from complex coherency samples
%
% INPUT
%   C12_vec : K x 1 complex coherency samples (e.g. from tapers/trials)
%
% OUTPUT
%   phi_hat : circular mean phase (radians), in [-pi, pi)
%   phi_std : jackknife std of phase (radians)
%   phi_CI  : 95% CI half-width (radians), so CI ≈ phi_hat ± phi_CI

    C12_vec = C12_vec(:);
    K = numel(C12_vec);

    if K < 2
        error('Need at least 2 samples for jackknife.');
    end

    % phases of each sample
    phi_k = angle(C12_vec);             % K x 1

    % full-sample circular mean (point estimate)
    phi_hat = angle(mean(exp(1i*phi_k)));

    % ---- jackknife: leave-one-out circular mean ----
    phi_jk = zeros(K, 1);
    for i = 1:K
        mask = true(K,1);
        mask(i) = false;

        % circular mean of remaining samples
        phi_jk(i) = angle(mean(exp(1i*phi_k(mask))));
    end

    % need to wrap differences to avoid 2π issues
    diffs = wrapToPi(phi_jk - angle(mean(exp(1i*phi_jk))));

    % jackknife variance
    var_jk = (K-1)/K * sum(diffs.^2);
    phi_std = sqrt(var_jk);

    % 95% CI (Gaussian approx on phase)
    phi_CI = 1.96 * phi_std;
end
