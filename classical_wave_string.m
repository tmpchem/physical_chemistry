% prepares a movie of U(x,t) over time with specified conditions
function [] = classical_wave_string(A, v, L)
  x_step = 0.01; % resolution of spatial grid
  t_step = 0.02; % time between movie frames
  x = linspace(0, L, L/x_step + 1);
  t_end = 2*L/v; % complete one full cycle
  n_frames = t_end/t_step + 1;
  t = 0; % start at time zero
  % calculate figure plot for each frame
  for frame_num = 1:n_frames
      U = 0; % initialize amplitude to zero
      % add amplitude from each normal mode
      for n = 1:length(A)
          U = U + A(n)*cos(n*pi*v*t/L)*sin(n*pi*x/L);
      end
      t = t + t_step; % propogate time
      plot(x, U); % prepare this frame's figure
      axis([0 L -sum(A) sum(A)]);
      xlabel('Position (x / L)');
      ylabel('Wave Amplitude U(x,t) / A');
      text(0.18 * L, 0.9 * sum(A), ...
          'Classical Wave Equation for a Vibrating String');
      movie_frames(frame_num) = getframe; % add frame to movie
  end
  % play resulting movie
  n_cycles = 50;
  framerate = 60;
  movie(movie_frames, n_cycles, framerate);
end