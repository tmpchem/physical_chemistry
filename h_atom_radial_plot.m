% prepares a plot of hydrogen atom radial wavefunction (R_nl(r))
% and radial distribution function (r^2 * |R_nl(r)|^2)
function [] = h_atom_radial_plot(n, l, Z)
    ao = 1; % functions plotted in units of Bohr radius
    r_min = 0;
    r_max = 25*ao;
    r_range = r_max - r_min;
    r_step = 0.005; % resolution of spatial grid
    r = linspace(0, r_range, r_range/r_step + 1);
    zeroline = zeros(r_range/r_step + 1, 1);
    psi_max = 2 * (Z / (2*ao))^(3/2); % maximum y value
    rho = Z * r / ao; % scaled length rho
    
    % associated Laguerre polynomials
    if (n == 1 && l == 0)
        L = 1;
        N = 2 * (Z/(ao))^(3/2);
    elseif (n == 2 && l == 0)
        L = 2 - rho;
        N = (Z/(2*ao))^(3/2);
    elseif (n == 2 && l == 1)
        L = rho;
        N = 1/sqrt(3) * (Z/(2*ao))^(3/2);
    elseif (n == 3 && l == 0)
        L = 27 - 18*rho + 2*rho.^2;
        N = 2/27 * (Z/(3*ao))^(3/2);
    elseif (n == 3 && l == 1)
        L = rho .* (6 - rho);
        N = 1/27 * (2*Z/(3*ao))^(3/2);
    elseif (n == 3 && l == 2)
        L = rho.^2;
        N = 4/(27*sqrt(10)) * (Z/(3*ao))^(3/2);
    end
    
    % prepare plot
    psi = N * L .* exp(-rho/n); % psi_n_l(r)
    psi2 = r.^2 .* psi.^2; % psi*_nl(r)*psi_nl(r)
    plot(r, r.^2, r, psi, r, psi2, r, zeroline); % plot psi and r^2*|psi|^2
    axis([r_min r_max -psi_max psi_max]);
    xlabel('Radius r (units of a_o)');
    ylabel('Function Value');
    legend(sprintf('r^2'), ...
        sprintf('R_%i_%i(r)', n, l), ...
        sprintf('r^2[R*_%i_%i(r)][R_%i_%i(r)]', n, l, n, l), ...
        'Location','SouthEast');
    text(0.18*r_range + r_min, 0.9*psi_max, ... % title
        'Hydrogen Atom Radial Wavefunction and R.D.F.');
    text(0.12*r_range + r_min, -0.8*psi_max, ... % energy
        sprintf('E_%i = -(1/%i)*(m_ee^4) / (8E_o^2h^2)', n, n^2));
end