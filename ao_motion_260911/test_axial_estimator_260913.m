function test_axial_estimator_260913()
%TEST_AXIAL_ESTIMATOR_260913  Synthetic ground-truth tests for the axial estimator.
%
% No rig, no data. Builds a synthetic volume whose true z displacement is known
% exactly, then asks whether the estimator recovers it. The point is to catch
% sign errors, off-by-ones and dilution bias BEFORE they are hiding inside a
% real measurement, where a wrong slope looks like biology.
%
% Runqi Zhang / 2026-09-13.

rng(7);
fprintf('\n================ axial estimator, synthetic tests ================\n');

%% --------- a synthetic volume: blobs at random 3-D positions ---------
ny = 128; nx = 128;
zTrue = -10:1.25:10;                 % same sampling as block B
K = numel(zTrue);
nBlob = 60;
bx = rand(1,nBlob)*nx;  by = rand(1,nBlob)*ny;  bz = (rand(1,nBlob)-0.5)*40;
amp = 200 + 800*rand(1,nBlob);
sxy = 2.5;  sz = 3.0;                % lateral / axial PSF sigma, px and um

    function I = render(zPlane, dxPx, dyPx)
        [X, Y] = meshgrid(1:nx, 1:ny);
        I = zeros(ny,nx);
        for b = 1:nBlob
            az = amp(b) * exp(-((zPlane - bz(b))^2) / (2*sz^2));
            if az < 1, continue; end
            I = I + az * exp(-(((X-(bx(b)+dxPx)).^2 + (Y-(by(b)+dyPx)).^2)) / (2*sxy^2));
        end
        I = I + 50;                  % PMT baseline
    end

refVol = zeros(ny,nx,K,'single');
for k = 1:K, refVol(:,:,k) = render(zTrue(k), 0, 0); end

E = dm_axial_estimator_260913('prep', refVol, zTrue, struct('nStrips',4));
fprintf('volume %dx%dx%d, %d blobs, planes %.2f to %.2f um\n', ny,nx,K,nBlob,zTrue(1),zTrue(end));

%% --------- TEST 1: pure axial displacement, noiseless ---------
fprintf('\n--- 1. pure axial displacement, no noise ---\n');
test = -4:0.5:4;
est = nan(size(test));
for i = 1:numel(test)
    I = single(render(test(i), 0, 0));
    o = dm_axial_estimator_260913('run', E, I);
    est(i) = median(o.dz(:));
end
p = polyfit(test, est, 1);
res = est - polyval(p, test);
fprintf('  estimated = %.4f * true %+.3f um,  residual RMS %.3f um\n', p(1), p(2), sqrt(mean(res.^2)));
report('slope within 5%% of 1', abs(p(1)-1) < 0.05);
report('bias below 0.2 um',     abs(p(2)) < 0.2);
report('residual below 0.2 um', sqrt(mean(res.^2)) < 0.2);

%% --------- TEST 2: SIGN. this is the one that runs a loop away ---------
fprintf('\n--- 2. sign convention ---\n');
I = single(render(+3, 0, 0));
o = dm_axial_estimator_260913('run', E, I);
z = median(o.dz(:));
fprintf('  a frame rendered at TRUE z = +3 um is estimated at %+.2f um\n', z);
report('sign is correct (positive)', z > 1);
fprintf('  => the correction to send is -dz = %+.2f um. Send +dz and the loop runs away.\n', -z);

%% --------- TEST 3: lateral motion must NOT leak into dz ---------
fprintf('\n--- 3. lateral rejection: 6 px of pure lateral, zero true dz ---\n');
leak = nan(1,7);  lat = -6:2:6;
for i = 1:numel(lat)
    I = single(render(0, lat(i), lat(i)/2));
    o = dm_axial_estimator_260913('run', E, I);
    leak(i) = median(o.dz(:));
end
fprintf('  dz reported for lateral shifts %s:\n    %s um\n', mat2str(lat), num2str(leak,'%+.2f '));
report('lateral leakage below 0.3 um', max(abs(leak)) < 0.3);

%% --------- TEST 4: brightness change must NOT read as defocus ---------
fprintf('\n--- 4. gain rejection: +-30%% brightness, zero true dz ---\n');
gains = [0.7 0.85 1.0 1.15 1.3];
gz = nan(size(gains));
for i = 1:numel(gains)
    I = single(gains(i) * render(0,0,0));
    o = dm_axial_estimator_260913('run', E, I);
    gz(i) = median(o.dz(:));
end
fprintf('  dz reported for gains %s:\n    %s um\n', mat2str(gains), num2str(gz,'%+.2f '));
report('brightness leakage below 0.3 um', max(abs(gz)) < 0.3);
fprintf('  (drop the gain column from the design matrix and this test fails --\n');
fprintf('   that is the whole reason it is there.)\n');

%% --------- TEST 5: shot noise ---------
fprintf('\n--- 5. Poisson noise at a realistic photon count ---\n');
for scale = [0.05 0.2 1.0]
    err = nan(1,20);
    for r = 1:20
        zt = (rand-0.5)*6;
        I = render(zt,0,0) * scale;
        I = single(poissrnd_local(max(I,0)) / scale);
        o = dm_axial_estimator_260913('run', E, I);
        err(r) = median(o.dz(:)) - zt;
    end
    fprintf('  photon scale %.2f : bias %+.3f um, scatter %.3f um\n', ...
            scale, mean(err), std(err));
end

%% --------- TEST 6: THE DILUTION BIAS a blurred reference causes ---------
fprintf('\n--- 6. what an UNGATED (breath-blurred) reference costs you ---\n');
fprintf('  This is the reason the reference is breath-gated. A reference built\n');
fprintf('  from frames spread over the breath cycle is axially blurred, dR/dz\n');
fprintf('  shrinks, and every dz estimate is pulled toward zero.\n');
for blurUm = [0 1 2 4]
    volB = zeros(ny,nx,K,'single');
    for k = 1:K
        acc = zeros(ny,nx);
        offs = linspace(-blurUm/2, blurUm/2, 9);
        for q = offs, acc = acc + render(zTrue(k)+q, 0, 0); end
        volB(:,:,k) = acc / numel(offs);
    end
    Eb = dm_axial_estimator_260913('prep', volB, zTrue, struct('nStrips',4));
    tt = -3:1:3;  ee = nan(size(tt));
    for i = 1:numel(tt)
        I = single(render(tt(i),0,0));
        o = dm_axial_estimator_260913('run', Eb, I);
        ee(i) = median(o.dz(:));
    end
    pb = polyfit(tt, ee, 1);
    fprintf('  reference blurred over %.0f um -> recovered slope %.3f', blurUm, pb(1));
    if blurUm == 0, fprintf('   <- gated\n'); else, fprintf('   (%.0f%% under-read)\n', 100*(1-pb(1))); end
end

%% --------- TEST 6b: the mechanism that ACTUALLY dilutes ---------
fprintf('\n--- 6b. a NOISY reference (few frames per plane) ---\n');
fprintf('  Test 6 says smoothing alone does not bias the slope -- correct, because\n');
fprintf('  least squares with a smoothed but NOISE-FREE predictor is consistent.\n');
fprintf('  Regression dilution needs NOISE in the predictor. That is what a gate\n');
fprintf('  costs you: it throws frames away, so the reference gets noisier.\n');
fprintf('  This is the real trade-off, and it runs the OTHER way from the gate.\n');
photonScale = 0.1;
for nAvg = [40 20 10 5 2]
    volN = zeros(ny,nx,K,'single');
    for k = 1:K
        acc = zeros(ny,nx);
        for q = 1:nAvg
            acc = acc + poissrnd_local(max(render(zTrue(k),0,0)*photonScale,0))/photonScale;
        end
        volN(:,:,k) = acc / nAvg;
    end
    En = dm_axial_estimator_260913('prep', volN, zTrue, struct('nStrips',4));
    tt = -3:1:3;  ee = nan(size(tt));
    for i = 1:numel(tt)
        I = single(poissrnd_local(max(render(tt(i),0,0)*photonScale,0))/photonScale);
        o = dm_axial_estimator_260913('run', En, I);
        ee(i) = median(o.dz(:));
    end
    pn = polyfit(tt, ee, 1);
    fprintf('  %2d frames/plane -> slope %.3f  (%+.0f%% error), resid RMS %.3f um\n', ...
            nAvg, pn(1), 100*(pn(1)-1), sqrt(mean((ee-polyval(pn,tt)).^2)));
end

%% --------- TEST 7: phase_shift sign and off-by-one ---------
fprintf('\n--- 7. phase_shift: the off-by-one and the sign ---\n');
A = single(render(0,0,0));
B = circshift(A, [5 -3]);            % B is A displaced by (+5, -3)
d = dm_axial_estimator_260913('shift', double(A), double(B));
fprintf('  B is A displaced by [+5 -3]; phase_shift returns [%+d %+d]\n', d(1), d(2));
report('returns the CORRECTION, i.e. [-5 +3]', isequal(d, [-5 3]));
report('circshift(B,d) restores A', isequal(circshift(B,d), A));

%% --------- TEST 8: DYNAMIC FLUORESCENCE (the real-prep problem) ---------
fprintf('\n--- 8. dynamic fluorescence: does activity leak into dz? ---\n');
fprintf('  Frames at TRUE dz = 0, with a fraction of blobs firing a transient.\n');
fprintf('  Any non-zero dz here is activity being misread as axial motion.\n');
    function I = render_active(zPlane, actIdx, actAmp)
        [X, Y] = meshgrid(1:nx, 1:ny);
        I = zeros(ny,nx);
        for b = 1:nBlob
            az = amp(b) * exp(-((zPlane - bz(b))^2) / (2*sz^2));
            if ismember(b, actIdx), az = az * actAmp; end
            if az < 1, continue; end
            I = I + az * exp(-(((X-bx(b)).^2 + (Y-by(b)).^2)) / (2*sxy^2));
        end
        I = I + 50;
    end

for fracActive = [0.05 0.15 0.40]
    nA = max(1, round(fracActive*nBlob));
    lk = nan(2,15);
    for r = 1:15
        actIdx = randperm(nBlob, nA);
        I = single(render_active(0, actIdx, 3.0));      % 3x transient
        for rb = 1:2
            Er = dm_axial_estimator_260913('prep', refVol, zTrue, ...
                      struct('nStrips',4,'robust', rb==2));
            o = dm_axial_estimator_260913('run', Er, I);
            lk(rb,r) = median(o.dz(:));
        end
    end
    fprintf('  %3.0f%% of cells firing 3x: leak %+.3f um plain, %+.3f um robust\n', ...
            100*fracActive, mean(abs(lk(1,:))), mean(abs(lk(2,:))));
end

%% --------- TEST 9: AXIAL DEGENERACY -- which vessels carry z information ---
fprintf('\n--- 9. penetrating vs in-plane vessels: the z aperture problem ---\n');
fprintf('  A vessel running straight down the z axis looks IDENTICAL at every\n');
fprintf('  plane, so it carries NO axial information no matter how bright it is.\n');
    function V = vessel_vol(mode)
        V = zeros(ny,nx,K,'single');
        [X, Y] = meshgrid(1:nx, 1:ny);
        for k = 1:K
            im = zeros(ny,nx);
            switch mode
                case 'penetrating'      % cylinders along z: z-invariant
                    for v = 1:8
                        cxv = 15 + 14*v;  cyv = 30 + mod(v*37, 70);
                        im = im + 900*exp(-(((X-cxv).^2 + (Y-cyv).^2))/(2*3^2));
                    end
                case 'inplane'          % tubes lying in the focal plane at various depths
                    for v = 1:8
                        zc = -9 + 2.4*v;
                        yv = 12 + 14*v;
                        prof = exp(-((zTrue(k)-zc)^2)/(2*sz^2));
                        im = im + 900*prof*exp(-((Y-yv).^2)/(2*3^2));
                    end
            end
            V(:,:,k) = im + 50;
        end
    end

for mode = {'penetrating','inplane'}
    V = vessel_vol(mode{1});
    Ev = dm_axial_estimator_260913('prep', V, zTrue, struct('nStrips',4));
    % axial information = how fast the image changes per um, relative to noise
    g = double(Ev.Rz(:,:,Ev.kMid));
    r0 = double(V(:,:,Ev.kMid));
    snr = sqrt(mean(g(:).^2)) / sqrt(mean(max(r0(:),1)));
    tt = -2:1:2;  ee = nan(size(tt));
    for i = 1:numel(tt)
        zi = tt(i);
        kk = interp1(zTrue, 1:K, zi);
        I = single((1-mod(kk,1))*V(:,:,floor(kk)) + mod(kk,1)*V(:,:,ceil(kk)));
        o = dm_axial_estimator_260913('run', Ev, I);
        ee(i) = median(o.dz(:));
    end
    pv = polyfit(tt, ee, 1);
    fprintf('  %-12s : axial SNR %6.3f, recovered slope %.3f\n', mode{1}, snr, pv(1));
end
fprintf('  => pick a FOV with vessels running ACROSS the plane, or branch points.\n');
fprintf('     A field of pure penetrating vessels is axially blind.\n');

fprintf('\n================================================================\n');
end

function report(name, pass)
if pass, fprintf('  [PASS] %s\n', name);
else,    fprintf('  [FAIL] %s   <<<<<<<<<<\n', name);
end
end

function y = poissrnd_local(lam)
% Poisson draw without the Statistics Toolbox. Knuth below 30, normal above.
y = zeros(size(lam));
small = lam < 30;
L = exp(-lam(small));
k = zeros(nnz(small),1);  p = ones(nnz(small),1);
done = false(size(k));
for it = 1:200
    p = p .* rand(size(p));
    adv = ~done & p > L(:);
    k(adv) = k(adv) + 1;
    done = done | (p <= L(:));
    if all(done), break; end
end
y(small) = k;
y(~small) = round(lam(~small) + sqrt(lam(~small)).*randn(nnz(~small),1));
y = max(y, 0);
end
