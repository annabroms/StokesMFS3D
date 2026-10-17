function approved = mfs_confirm_system_sizes(systemSizes,systemNames,confirmationRequired)
%MFS_CONFIRM_SYSTEM_SIZES Shared Laplace/Stokes size preview and consent gate.
% systemSizes has one [equations,unknowns] row per solve.  Report the matrix
% alone, not peak memory: kernel temporaries and SVD factors need more space.
% No kernel assembly or solve is performed here.  Default is explicit consent.

if nargin < 3
    confirmationRequired = true;
end
systemNames = string(systemNames);
fprintf('\nDense-system preview (both particles together):\n');
for k = 1:size(systemSizes,1)
    matrixMiB = 8*prod(double(systemSizes(k,:)))/2^20;
    fprintf('  %s: %d equations x %d unknowns; matrix alone %.1f MiB\n', ...
        systemNames(k),systemSizes(k,1),systemSizes(k,2),matrixMiB);
end
fprintf('  Kernel temporaries, SVD factors and workspace require additional memory.\n');
if any(systemSizes(:,1) <= systemSizes(:,2))
    error('mfs_confirm_system_sizes:TooFewEquations', ...
        ['Each system must have more equations than unknowns. Increase basic ', ...
         'collocation count/resolution or patch refinement, or reduce sources/families.']);
end
if ~confirmationRequired
    approved = true; % explicit opt-out for an already-approved automated run
    return
end
try
    answer = input('Proceed with the displayed systems, SVDs and evaluations? [y/N]: ','s');
    approved = any(strcmpi(strtrim(answer),{'y','yes'}));
catch
    % Batch/closed input is not consent.  Stop rather than starting a solve.
    approved = false;
    fprintf('No interactive approval is available.\n');
end
if ~approved
    fprintf('Cancelled: no kernel assembly, solve or evaluation will start for these systems.\n');
end
end
