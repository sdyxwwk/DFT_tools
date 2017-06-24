function [G00, Gb] = greensurf1985(E, eta, H00, H01, H10)
    %
    % A fast iterative algorithm to calculate G00.
    % G00 is used to calculate LDOS: rho = -1./pi .* trace(imag(G00)
    %
    % Input:
    %    E (1,1) : energy
    %    eta (1,1) : imaginary part of E, (Emax-Emin)./num_E
    %    H00 (N,N) : intralayer Hamiltonian
    %    H01 (N,N) : interlayer coupling
    %    H10 (N,N) : H01'
    % Output:
    %    G00 (N,N) : surface Green function
    %    Gb (N,N) : Green function for a bulk layer
    %
    % Reference: M. L. Sancho, et al. J. Phys. F: Metal Physics, 15(4), 851 (1985).
    %
    
    %eta = 1e-12 ;
    accuracy = 1e-12 ;
    maxiter = 1000 ;
    
    alpha = H01 ;
    beta  = H10 ;
    epsilon  = H00 ;
    epsilons = H00 ;
    epsilont = H00 ;
    E_p = E + j.*eta ;
    Nsize = size(H00, 1) ;
    
    niter = 0 ;
    E_inv = inv(E_p.*eye(Nsize)-epsilon) ;
    % E_alpha = norm(E_inv * alpha) ;
    % E_beta  = norm(E_inv * beta) ;
    while((niter<maxiter) && (norm(alpha)>=accuracy || norm(beta)>=accuracy))
        epsilon  = epsilon  + alpha * E_inv * beta + beta * E_inv * alpha ;
        epsilons = epsilons + alpha * E_inv * beta ;
        epsilont = epsilont + beta * E_inv * alpha ;
        
        alpha = alpha * E_inv * alpha ;
        beta = beta * E_inv * beta ;
        
        E_inv = inv(E_p.*eye(Nsize) - epsilon) ;
        % E_alpha = norm(E_inv * alpha) ;
        % E_beta  = norm(E_inv * beta) ;
        niter = niter + 1 ;
    end
    G00 = inv(E_p.*eye(Nsize) - epsilons) ;
    Gb = inv(E_p.*eye(Nsize) - epsilon) ;
end
