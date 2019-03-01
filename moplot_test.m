%% test for moplot
clear; clc;
% H2O
% Load test data file
load('motest','testdata');
test = testdata(6);
atoms = test.atoms;
xyz_a0 = test.xyz_a0;
charge = test.charge;
out = mocalc(atoms, xyz_a0, charge, test);
%%
iMO = 5;
level = 0.45;
moplot(atoms, xyz_a0, out, iMO, level);