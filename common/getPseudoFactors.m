function [Y, U] = getPseudoFactors(N, tol, visualise)
%GETPSEUDOFACTORS Computes factors that provide the matrix psuedoinverse from 
% a truncated SVD 
%
% Syntax:
%   [Y, U] = getPseudoFactors(N, tol, visualise)
%
% Inputs:
%   N         - Input matrix for which to compute a pseudoinverse
%   tol       - Relative truncation threshold; singular values σ are kept if σ > max(σ) * tol
%   visualise - Logical flag: plot singular values (true/false)
%
% Outputs:
%   U - Matrix of left singular vectors corresponding to retained singular values
%   Y - Product VS⁺, where:
%         - S⁺ is a diagonal matrix with entries 1/σ for retained singular values
%         - V contains the corresponding right singular vectors
%
% Description:
%   Computes a truncated SVD of matrix N and returns factor matrices U and Y such that:
%       N⁺ ≈ Y * U'
%   This allows efficient and backward-stable application of the pseudoinverse without explicitly forming it.
%
% Notes:
%   - Only singular values greater than max(σ) * tol are retained
%   - Intended for use in stable pseudoinverse application (e.g., solving least-squares problems)
%
% Anna Broms 4 April 2025

% - With N the discretized Stokeslets, 
%with tol = eps, we throw away only the last singval, corresponding to the
%non-trivial null-space (normal direction for single layer...). Does that make sense?
%yes! testing ra = size(N,2) leads to a horribly irregular 1-body basis. 


if nargin < 3
    visualise = 0; 
end

[UU,S,V] = svd(N);
S = diag(S);

%use relative tolerance 
ra = sum(S>max(S)*tol); 


if visualise
   
    figure(57)
   % clf;
    semilogy(S,'o-');
    hold on
    semilogy(ra*ones(1,2),logspace(-15,5,2),'r--')
    c = max(S)/min(S);  %Condition number
    str = sprintf('Self condition number %1.3e, max sing %1.3e',c,max(S));
    title(str,'interpreter','latex');
    grid on
    xlabel('$j$','interpreter','latex')
    ylabel('$\sigma_j$','interpreter','latex')

    % figure(56)
    % SV = SS*V';
    % semilogy(abs(SV(:,end-50:end)'))
    % 
    % 
    % figure(56)
    % semilogy(abs(UU(end-3:end,:)'))

    % figure(57)
    % semilogy(abs(diff(S)));
    % hold on
    % title('Decay rate of sing vals','interpreter','latex')
    % 
    % 
    % figure(58)
    % semilogy(abs(diff(S)./S(2:end)));
    % hold on
    % title('Relative decay rate of sing vals','interpreter','latex')
end

S = S(1:ra);  %get pseudoinverse of S
iS = 1./S; 
Y = V(:,1:ra)*diag(iS); 
U = UU(:,1:ra); 





end