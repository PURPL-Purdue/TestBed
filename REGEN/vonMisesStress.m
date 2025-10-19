function sigma_vm = vonMisesStress(varargin)
% VONMISESSTRESS Compute von Mises equivalent stress (2D plane-stress or 3D)
%
% Usage:
%   sigma_vm = vonMisesStress(sxx, syy, sxy)                % 2D plane stress
%   sigma_vm = vonMisesStress(sxx, syy, szz, sxy, syz, sxz) % full 3D
%
% Inputs may be scalars, vectors or arrays of identical size (vectorized).
% Units: Pascals (Pa)
%
% Examples:
%   vm = vonMisesStress(100e6, 0, 0); % uniaxial 100 MPa -> vm = 100e6
%   vm = vonMisesStress(sxx, syy, sxy); % vectorized inputs

n = nargin;
if n == 3
    sxx = varargin{1}; syy = varargin{2}; sxy = varargin{3};
    szz = 0; syz = 0; sxz = 0;
elseif n == 6
    sxx = varargin{1}; syy = varargin{2}; szz = varargin{3};
    sxy = varargin{4}; syz = varargin{5}; sxz = varargin{6};
else
    error('vonMisesStress:InvalidInput','Call with 3 (plane stress) or 6 (3D) components.');
end

% Basic size checks (vectorized)
if ~isequal(size(sxx), size(syy)) || ~isequal(size(sxx), size(szz)) || ...
   ~isequal(size(sxx), size(sxy)) || ~isequal(size(sxx), size(syz)) || ...
   ~isequal(size(sxx), size(sxz))
    error('vonMisesStress:SizeMismatch','All inputs must have the same size.');
end

% --- Computation: von Mises (engineering shear convention) ---
% 3D formula:
%   sigma_vm = sqrt( 0.5*((sxx-syy).^2 + (syy-szz).^2 + (szz-sxx).^2) + 3*(sxy.^2 + syz.^2 + sxz.^2) );
sigma_vm = sqrt( 0.5*((sxx - syy).^2 + (syy - szz).^2 + (szz - sxx).^2) ...
                 + 3*(sxy.^2 + syz.^2 + sxz.^2) );
end
