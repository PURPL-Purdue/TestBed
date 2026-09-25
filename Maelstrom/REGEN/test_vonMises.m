% test_vonMises.m - quick checks for vonMisesStress
% 1) Uniaxial stress: sxx = 100e6, others zero -> vm = 100e6
% 2) Pure shear: sxy = 50e6, plane stress -> vm = sqrt(3)*sxy

% Test 1 - uniaxial
sxx = 100e6; syy = 0; sxy = 0;
vm1 = vonMisesStress(sxx, syy, sxy);
fprintf('Uniaxial: computed VM = %g Pa (expected %g Pa)\n', vm1, sxx);

% Test 2 - pure shear (plane stress)
sxx = 0; syy = 0; sxy = 50e6;
vm2 = vonMisesStress(sxx, syy, sxy);
expected_vm2 = sqrt(3)*sxy;
fprintf('Pure shear: computed VM = %g Pa (expected %g Pa)\n', vm2, expected_vm2);

% Simple pass/fail
tol = 1e-6;
if abs(vm1 - 100e6) < tol && abs(vm2 - expected_vm2) < 1e-3
    fprintf('Basic tests passed.\n');
else
    fprintf('Basic tests failed.\n');
end
