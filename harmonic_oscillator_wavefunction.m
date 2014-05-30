% prepares a movie of Psi(x,t) over time with specified conditions
function [] = harmonic_oscillator_wavefunction(C, k, mu, e_scale)
  x_step = 0.01; % resolution of spatial grid
  t_step = 0.02; % time between movie frames
  width = 4;
  x = linspace(-width/2, width/2, width/x_step + 1);
  V = 0.02 * k * x.^2;
  N = 0:(length(C)-1); % array of consecutive integers
  % harmonic oscillator functions
  alpha = sqrt(k * mu); % alpha value
  omega = sqrt(k / mu); % angular frequency
  % energy values
  E = e_scale * omega * (N + 1/2);
  % normalization constant
  Nnorm = (2.^N .* factorial(N)).^(-1/2) * (alpha / pi)^(1/4);
  % hermite polynomials
  xi = sqrt(alpha) * x;
  H1 = 1;
  H2 = 2*xi;
  H3 = 4*xi.^2 - 2;
  H4 = 8*xi.^3 - 12*xi;
  H5 = 16*xi.^4 - 48*xi.^2 + 12;
  H = [H1, H2, H3, H4, H5];
  % time dependence
  t_end = 4*pi/(e_scale * omega); % complete one full cycle
  n_frames = t_end/t_step + 1;
  t = 0; % start at time zero
  % calculate figure plot for each frame
  for frame_num = 1:n_frames
      RePsi = 0; % real part
      ImPsi = 0; % imaginary part
      for n = 1:length(C) % add amplitude from each eigenfunction
          psi_n = Nnorm(n)*exp(-alpha*x.^2 / 2);
          if n == 1
              psi_n = 1 * psi_n;
          elseif n == 2
              psi_n = 2 * xi .* psi_n;
          elseif n == 3
              psi_n = (4*xi.^2 - 2) .* psi_n;
          elseif n == 4
              psi_n = (8*xi.^3 - 12*xi) .* psi_n;
          elseif n == 5
              psi_n = (16*xi.^4 - 48*xi.^2 + 12) .* psi_n;
          end
          RePsi = RePsi + C(n)*cos(E(n)*t)*psi_n;
          ImPsi = ImPsi + C(n)*sin(E(n)*t)*psi_n;
      end
      Psi2 = RePsi.^2 + ImPsi.^2;
      t = t + t_step; % propogate time
      plot(x, RePsi, x, ImPsi, x, Psi2, x, V); % this frame's figure
      axis([-width/2 width/2 -(sum(C.^2)) (sum(C.^2))]);
      xlabel('Bond Displacement (x)');
      ylabel('Function Amplitude');
      text(-0.4*width, 0.9*(sum(C.^2)), ...
          'Time Dependent Harmonic Oscillator Wavefunction Animation');
      movie_frames(frame_num) = getframe; % add frame to movie
  end
  % play resulting movie
  n_cycles = 50;
  framerate = 60;
  movie(movie_frames, n_cycles, framerate);
end