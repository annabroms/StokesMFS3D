repo_root = fileparts(mfilename('fullpath'));
if isempty(repo_root)
    repo_root = pwd;
end

addpath(genpath(repo_root)); % add this repository

stokes_direct_root = fullfile(repo_root,'..','Stokes_Direct');
if isfolder(stokes_direct_root)
    addpath(genpath(stokes_direct_root)); % to build, see https://github.com/annabroms/Stokes_Direct
end

fmm3d_root = fullfile(repo_root,'..','FMM3D');
if isfolder(fmm3d_root)
    addpath(genpath(fmm3d_root)); % to build, see https://fmm3d.readthedocs.io/en/latest/install.html
end
