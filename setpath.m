repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root)
    repo_root = pwd;
end

addpath(genpath(repo_root)); % add this repository

fmm3d_root = fullfile(repo_root,'..','FMM3D');
if isfolder(fmm3d_root)
    addpath(genpath(fmm3d_root)); % to build, see https://fmm3d.readthedocs.io/en/latest/install.html
end

% The half-space package is independent; add its MATLAB package directory only.
halfspace_root = fullfile(repo_root,'..','StokesHalfspaceFMM','matlab');
if isfolder(halfspace_root)
    addpath(halfspace_root);
end
