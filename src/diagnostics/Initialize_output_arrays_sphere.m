function [vec_Ene, vec_L_z, ...
          vec_position_b_North, vec_position_b_South, ...
          vec_t, ...
          vec_Ec, vec_Ei, vec_Ehd, vec_Eqp] = ...
    Initialize_output_arrays_sphere(n_iterations, sample_frequency)

%INITIALIZE_OUTPUT_ARRAYS_SPHERE Preallocate sampled diagnostics.
%
%   Allocates arrays for total energy, angular momentum, component-b peak
%   positions, sample times, and the compressible/incompressible kinetic-
%   energy decomposition.
%
%   Inputs
%   ------
%   n_iterations
%       Total number of real-time integration steps.
%   sample_frequency
%       Number of integration steps between stored samples.

    n_elements = floor(n_iterations/sample_frequency) + 1;

    % --- Standard outputs ---
    vec_Ene = zeros(1,n_elements);
    vec_L_z = zeros(1,n_elements);

    vec_position_b_North = zeros(n_elements,2);   % [theta, phi]
    vec_position_b_South = zeros(n_elements,2);   % [theta, phi]

    vec_t = zeros(1,n_elements);

    % --- New: hydrodynamic decomposition ---
    vec_Ec = zeros(1,n_elements);
    vec_Ei = zeros(1,n_elements);

    vec_Ehd = zeros(1,n_elements);
    vec_Eqp = zeros(1,n_elements);

end
