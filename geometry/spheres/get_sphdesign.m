function [X,w] = get_sphdesign(Nmax)
% """
%     [X, w] = get_sphdesign(Nmax)
% 
% Load then returns the largest spherical design (quadrature points on S^2) with
% no more than `Nmax` points. `X` is the (N,3) coordinate array in R^3
% of the N points, and `w`, a (N,) column vector of their corresponding
% weights, w.r.t. surface measure on S^2. (The weights are equal, thus easily calculated by the user anyway.)
% N<=Nmax, with N also <= 180, the largest design available in Wormseley's files.
% 
% Also see: [`getavailablesphdesigns`](@ref)

% The function is adopted from Alex's function in Julia Nov 13 2023
% """
    

    dirName = "spherical_designs";

    Ns = getavailablesphdesigns(dirName);
    t = find(Ns<=Nmax,1,'last');% degree
    assert(size(t,1)==1,"no available N are <= requested Nmax!");
    N = Ns(t);     %# the largest N not exceeding Nmax -- check indexing here!
    fnam = sprintf('sf%.3d.%.5d', t, N);     % reverse-engineer filename
    

    absfnam = sprintf('%s/%s',dirName,fnam);
    X = load(absfnam);    
    % X = readdlm(absfnam)
    assert(size(X,1)==N,"read wrong number of lines from file! Please see sphdesigns/README");
    w = ones(N,1)*(4*pi/double(N));

end

function res = getavailablesphdesigns(dirName)
    % file list with leading "sf" chars removed (dot is separator)...
    
    % (note linux, OSX, not Windows)
    fileName = sprintf('%s/filelist.txt',dirName);
    %fileName = "~/NOBACKUP/qbx_with_brownian/spherical_designs/filelist.txt";
    fileID = fopen(fileName, 'r');
    dataCell = textscan(fileID, '%*[^.].%d');
    
    % Close the file
    fclose(fileID);
    res = dataCell{1};

end

function hostname = get_hostname()
  [status, hostname] = system('hostname -s');
  if status ~= 0
    [status, hostname] = system('hostname');
  end
  if status ~= 0
    warning('UnknownHostName', 'Host name set to UNKNOWN');
    hostname = 'UNKNOWN';
  end
end


    