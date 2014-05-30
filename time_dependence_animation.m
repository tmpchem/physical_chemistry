% prepares a movie of Psi(x,t) over time with specified conditions
function [] = time_dependence_animation(C, e_scale)
  L = 1; % relative units
  x_step = 0.01; % resolution of spatial grid
  t_step = 0.02; % time between movie frames
  x = linspace(0, L, L/x_step + 1);
  t_end = 2*L/e_scale; % complete one full cycle
  n_frames = t_end/t_step + 1;
  t = 0; % start at time zero
  % calculate figure plot for each frame
  for frame_num = 1:n_frames
      RePsi = 0; % real part
      ImPsi = 0; % imaginary part
      for n = 1:length(C) % add amplitude from each eigenfunction
          psi_n = sqrt(2/L)*sin(n*pi*x/L);
          RePsi = RePsi + C(n)*cos(n^2 * pi*e_scale*t/L)*psi_n;
          ImPsi = ImPsi + C(n)*sin(n^2 * pi*e_scale*t/L)*psi_n;
      end
      Psi2 = RePsi.^2 + ImPsi.^2;
      t = t + t_step; % propogate time
      plot(x, RePsi, x, ImPsi, x, Psi2); % prepare this frame's figure
      axis([0 L -(2/L)*(sum(C.^2)) (2/L)*(sum(C.^2))]);
      xlabel('Position (x / L)');
      ylabel('Function Amplitude');
      text(0.12*L, 1.7*(sum(C.^2)), ...
          'Time Dependent Particle in a Box Wavefunction Animation');
      movie_frames(frame_num) = getframe; % add frame to movie
  end
  % play resulting movie
  n_cycles = 50;
  framerate = 60;
  movie(movie_frames, n_cycles, framerate);
end