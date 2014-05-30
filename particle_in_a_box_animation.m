% prepares a movie of particle in a box wavefunctions (psi_n(x))
% and probability distribution functions (|psi_n(x)|^2)
function [] = particle_in_a_box(max_n, framerate)
  L = 1; % functions plotted in units of L
  x_step = 0.005; % resolution of spatial grid
  x = linspace(0, L, L/x_step + 1);
  % calculate figure plot for each frame
  for n = 1:max_n
      psi = sqrt(2/L) * sin(n*pi*x/L); % psi_n(x)
      psi_star_psi = conj(psi) .* psi; % psi*_n(x)*psi_n(x)
      plot(x, psi, x, psi_star_psi); % plot psi and |psi|^2
      axis([0 L -2/L 2/L]);
      xlabel('Position x (units of L)');
      ylabel('Function Value (units of 1/L)');
      legend(sprintf('psi_%i(x) = sqrt(2/L)*sin(%ipi*x/L)',n,n), ...
          sprintf('|psi_%i(x)|^2 = (2/L)*sin^2(%ipi*x/L)',n,n), ...
          'Location','SouthEast');
      text(0.18 * L, 0.9 * 2/L, ... % title
          '1D Particle in a Box Wavefunction and P.D.F.');
      text(0.12 * L, -0.8 * 2/L, ... % energy
          sprintf('E_%i = %ih^2 / 8mL^2', n, n^2));
      movie_frames(n) = getframe; % add frame to movie
  end
  % play resulting movie
  n_cycles = 2;
  movie(movie_frames, n_cycles, framerate);
end