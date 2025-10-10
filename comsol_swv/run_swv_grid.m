% run_swv_grid.m
% Batch sweep of frequency and concentration using COMSOL LiveLink.

import com.comsol.model.*
import com.comsol.model.util.*

% 1) COMSOL connection:
% Either MATLAB via the "COMSOL with MATLAB" shortcut (the following script auto-connects)
% ModelUtil.connect;                 % if a COMSOL server is already running
% mphstart;                          % starts a local COMSOL server (easiest)

% 2) Loading model
model = mphload('aptamer_swv.mph');  % Note: must be in the current folder or provide the full path

% 3) Defining sweep vectors
f_vec = logspace(1,3,25);            % 10 Hz to 1 kHz
C_vec = logspace(-7, -4, 15);        % 0.1 µM to 100 µM, in mol/m^3

% 4) Preallocate results: rows = concentrations, cols = frequencies
I = zeros(numel(C_vec), numel(f_vec));

% 5) Parallel pool (8 workers if you like)
p = gcp('nocreate'); if isempty(p), parpool(8); end

% 6) Current expression
% Example: currentExpr = '-intI(j_s)';  % if j_s is the normal current density
% If the operator tag or variable differs, will be edited:
currentExpr = '-intI(j_s)';

parfor ic = 1:numel(C_vec)
    Ci = C_vec(ic);
    Imat = zeros(1, numel(f_vec));
    for jf = 1:numel(f_vec)
        fj = f_vec(jf);

        % Set parameters
        mphparam(model, 'Ctarget', Ci);
        mphparam(model, 'f', fj);

        % Solve Stationary study (tag 'std1' by default; needs to be adjusted if differs)
        mphsolve(model, 'std1');

        % Evaluate global current
        Imat(jf) = mphglobal(model, currentExpr);
    end
    I(ic,:) = Imat;
end

% 7) Save data
save('swv_grid.mat','f_vec','C_vec','I');

% 8) Quick plot
figure;
imagesc(log10(f_vec), log10(C_vec*1e6), I); axis xy;
xlabel('log10 f [Hz]'); ylabel('log10 C [\muM]'); title('I(f, C)');
colorbar;
