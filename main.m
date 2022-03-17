clc;
close all;
clear all;

currentPath = fileparts(mfilename('fullpath'));
% this path could be used to include some user-defined functions in future

L = 60;
NUM_f = 600;
% fracture size distribution
% uniform
min = 1; max = 10;

Center_f = zeros(NUM_f, 3);
R_f = zeros(NUM_f, 1);
Normal_f = zeros(NUM_f, 3);
Verts_f = zeros(NUM_f, 12);
alpha_f = zeros(NUM_f, 1);
beta_f = zeros(NUM_f, 1);

%----------generate DFN
Center_f = [-0.5 * L + L .* rand(NUM_f, 1), -0.5 * L + L .* rand(NUM_f, 1), -0.5 * L + L .* rand(NUM_f, 1)]; % uniform fracture center positions
R_f = min + (max - min) .* rand(NUM_f, 1); % uniform sizes
Normal_f = [-1 + 2 .* rand(NUM_f, 1), -1 + 2 .* rand(NUM_f, 1), -1 + 2 .* rand(NUM_f, 1)]; % uniform orientation
tem_Norm = vecnorm(Normal_f')';
Normal_f = Normal_f ./ tem_Norm; clear tem_Norm
% now vertices of square fractures
Verts_f(:, [1 2]) = [-1 + 2 .* rand(NUM_f, 1), -1 + 2 .* rand(NUM_f, 1)];
Verts_f(:, 3) = -1.0 .* (Normal_f(:, 1) .* Verts_f(:, 1) + Normal_f(:, 2) .* Verts_f(:, 2)) ./ Normal_f(:, 3);
Norm_k = vecnorm(Verts_f(:, [1 2 3])')';
Ratio_k = R_f ./ Norm_k;
Verts_f(:, [1 2 3]) = [Ratio_k, Ratio_k, Ratio_k] .* Verts_f(:, [1 2 3]);
Verts_f(:, [7 8 9]) = -Verts_f(:, [1 2 3]);
tem_Vec = cross(Verts_f(:, [1 2 3]) ./ vecnorm(Verts_f(:, [1 2 3])')', Normal_f, 2);
Norm_k = vecnorm(tem_Vec')';
Ratio_k = R_f ./ Norm_k;
Verts_f(:, [4 5 6]) = [Ratio_k, Ratio_k, Ratio_k] .* tem_Vec;
Verts_f(:, [10 11 12]) = -Verts_f(:, [4 5 6]);
Verts_f(:, [1 2 3]) = Verts_f(:, [1 2 3]) + Center_f;
Verts_f(:, [4 5 6]) = Verts_f(:, [4 5 6]) + Center_f;
Verts_f(:, [7 8 9]) = Verts_f(:, [7 8 9]) + Center_f;
Verts_f(:, [10 11 12]) = Verts_f(:, [10 11 12]) + Center_f; clear Ratio_k Norm_k tem_Vec
% now orientations of fractures
sign_tmp = Normal_f(:, 3) ./ abs(Normal_f(:, 3));
Normal_tmp = [sign_tmp .* Normal_f(:, 1), sign_tmp .* Normal_f(:, 2), sign_tmp .* Normal_f(:, 3)]; clear sign_tmp
alpha_f = atan2(Normal_tmp(:, 2), Normal_tmp(:, 1));
A = find(alpha_f < 0);
alpha_f(A) = alpha_f(A) + 2 * pi;
beta_f = acos(Normal_tmp(:, 3) ./ vecnorm(Normal_tmp')'); clear A Normal_tmp

figure(1)
cube_frame = [-0.5 * L, -0.5 * L, 0.5 * L; -0.5 * L, 0.5 * L, 0.5 * L; 0.5 * L, 0.5 * L, 0.5 * L; 0.5 * L -0.5 * L, 0.5 * L; -0.5 * L, -0.5 * L, -0.5 * L; -0.5 * L, 0.5 * L, -0.5 * L; 0.5 * L, 0.5 * L, -0.5 * L; 0.5 * L -0.5 * L, -0.5 * L; -0.5 * L, 0.5 * L, 0.5 * L; -0.5 * L, 0.5 * L, -0.5 * L; -0.5 * L, -0.5 * L, -0.5 * L; -0.5 * L, -0.5 * L, 0.5 * L; 0.5 * L, 0.5 * L, 0.5 * L; 0.5 * L, 0.5 * L, -0.5 * L; 0.5 * L, -0.5 * L, -0.5 * L; 0.5 * L, -0.5 * L, 0.5 * L; 0.5 * L, -0.5 * L, 0.5 * L; 0.5 * L, -0.5 * L, -0.5 * L; -0.5 * L, -0.5 * L, -0.5 * L; -0.5 * L, -0.5 * L, 0.5 * L; 0.5 * L, 0.5 * L, 0.5 * L; 0.5 * L, 0.5 * L, -0.5 * L; -0.5 * L, 0.5 * L, -0.5 * L; -0.5 * L, 0.5 * L, 0.5 * L];
figure(1); view(3); title('Discete fracture network'); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); hold on
patch('Vertices', cube_frame, 'Faces', [1, 2, 3, 4; 5 6 7 8; 9 10 11 12; 13 14 15 16], 'FaceVertexCData', zeros(size(cube_frame, 1), 1), 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 0); hold on
Verts_show = zeros(NUM_f * 4, 3);
Verts_show([1:4:NUM_f * 4], :) = Verts_f(:, [1 2 3]);
Verts_show([2:4:NUM_f * 4], :) = Verts_f(:, [4 5 6]);
Verts_show([3:4:NUM_f * 4], :) = Verts_f(:, [7 8 9]);
Verts_show([4:4:NUM_f * 4], :) = Verts_f(:, [10 11 12]);
patch('Vertices', Verts_show, 'Faces', reshape([1:1:NUM_f * 4], [4 NUM_f])', 'FaceVertexCData', [Verts_show(:, 3)], 'FaceColor', 'interp', 'EdgeAlpha', 1, 'facealpha', 1); hold on;
axis([-1.1 * 0.5 * L, 1.1 * 0.5 * L, -1.1 * 0.5 * L, 1.1 * 0.5 * L, -1.1 * 0.5 * L, 1.1 * 0.5 * L]);
view(3)
figure(2)
polarscatter(alpha_f, beta_f, 's', 'filled'); rlim([0 0.5*pi]);
rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 * pi / 12, 6 * pi / 12 ]);
title(['Fractures', '''',' orientations']); hold on
set(gca,'rticklabel',[]);
