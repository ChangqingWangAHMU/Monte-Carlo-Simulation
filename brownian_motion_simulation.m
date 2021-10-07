function [x dx] = brownian_motion_simulation (n,d,dt)

%*****************************************************************************80
%% BROWNIAN_MOTION_SIMULATION simulates Brownian motion.
%  Parameters:
%    Input, integer n, the number of time steps to take, plus 1. 
%    Input, real d, the diffusion coefficient.
%    Input, real dt, the time step;
%    Output, x and dx, the position at every time and interval between two
%    times
%
  rng('shuffle');
%  Compute the individual steps.
  x = zeros (n,3);
%  Direction is random.
    a = randn ( 3, n - 1 );
    v = sqrt ( 6*d * dt ) * randn ( 1, n - 1 )./ sqrt ( sum ( a.^2 ) );
    b = spdiags ( v', 0, n-1, n-1 );
    dx(1:3,1:n-1) = a * b;
%  Each position is the sum of the previous steps.
  x(2:n,1:3) = cumsum ( dx(1:3,1:n-1), 2 )';
  return
end