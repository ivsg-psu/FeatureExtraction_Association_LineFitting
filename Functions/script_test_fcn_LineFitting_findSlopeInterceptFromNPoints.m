% script_test_fcn_LineFitting_findSlopeInterceptFromNPoints
% Exercises the function: fcn_LineFitting_findSlopeInterceptFromNPoints

% Revision history:
% 2020_06_25 
% -- wrote the code
% 2023_07_18 - Sean Brennan, sbrennan@psu.edu
% -- better comments, reorganization

close all;
clc;


%% Simplest example
points = [2 3; 4 5];
[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points);
fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);


%% Example 1: a basic test of the figure using simplest example
fig_num = 1;
points = [2 3; 4 5];
[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points,fig_num);
title(sprintf('Fig %.0d: simplest example',fig_num));
fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);


%% Example 2: many points randomly determined
fig_num = fig_num + 1;
slope = 3;
intercept = 4;

Npoints = 1000;
x_data = linspace(-2,5,Npoints)';
y_data = x_data*slope + intercept + intercept*0.2*randn(Npoints,1);
points = [x_data,y_data];

[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points,fig_num);
title(sprintf('Fig %.0d: random points showing regression fit',fig_num));

fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Example 3: a vertical line
fig_num = fig_num + 1;
points = [2 0; 2 2];
[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points,fig_num);
title(sprintf('Fig %.0d: vertical line example',fig_num));

fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 4: many vertical points
fig_num = fig_num + 1;

Npoints = 1000;
x_data = 2*ones(Npoints,1);
y_data = linspace(-1,10,Npoints)';
points = [x_data,y_data];

[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points,fig_num);
title(sprintf('Fig %.0d: many vertical points example',fig_num));

fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);

%% Test 5: Test a singular situation (determinant is zero - which gives b = 0)
% The A matrix is singular, but should still generate an answer.
fig_num = fig_num + 1;
points = [6 4; 3 2];
[slope,intercept] = fcn_LineFitting_findSlopeInterceptFromNPoints(points,fig_num);
title(sprintf('Fig %.0d: zero determinant example',fig_num));

fprintf(1,'Slope is: %.2f, Intercept is: %.2f\n',slope,intercept);


%% Fill in some test data - 45 degree line
slope = 1;
angle_in_radians = atan(slope);
intercept = 23;
x = (1:100)';
y = slope*x + intercept + randn(length(x(:,1)),1);
ENU_data_to_fit = [x y];

%% BASIC call - fit 45 degree line
fig_num = fig_num + 1;
[slope_fitted, intercept_fitted, angle_in_radians_fitted, ...
    unit_tangent_vector, unit_orthogonal_vector] = ...
    fcn_LineFitting_findSlopeInterceptFromNPoints(ENU_data_to_fit, fig_num); %#ok<*ASGLU> 
title(sprintf('Fig %.0d: 45 degree line example',fig_num));


percent_error = 0.02;
assert(abs(slope_fitted - slope)<slope*percent_error+0.01);
assert(abs(intercept_fitted - intercept)<intercept*percent_error+0.01);
assert(abs(angle_in_radians_fitted - angle_in_radians)<angle_in_radians*percent_error+0.01);
assert(abs(unit_tangent_vector(1,1) - cos(angle_in_radians))<cos(angle_in_radians)*percent_error+0.01);
assert(abs(unit_tangent_vector(1,2) - sin(angle_in_radians))<sin(angle_in_radians)*percent_error+0.01);

%% Fill in some test data - horizontal line
slope = 0;
angle_in_radians = atan(slope);
intercept = 23;
x = (1:100)';
y = slope*x + intercept + randn(length(x(:,1)),1);
ENU_data_to_fit = [x y];

%% BASIC call - fit horizontal line
fig_num = fig_num + 1;
[slope_fitted, intercept_fitted, angle_in_radians_fitted, ...
    unit_tangent_vector, unit_orthogonal_vector] = ...
    fcn_LineFitting_findSlopeInterceptFromNPoints(ENU_data_to_fit, fig_num);
title(sprintf('Fig %.0d: 0 degree line example',fig_num));

percent_error = 0.02;
assert(abs(slope_fitted - slope)<slope*percent_error+0.01);
assert(abs(intercept_fitted - intercept)<intercept*percent_error+0.01);
assert(abs(angle_in_radians_fitted - angle_in_radians)<angle_in_radians*percent_error+0.01);
assert(abs(unit_tangent_vector(1,1) - cos(angle_in_radians))<cos(angle_in_radians)*percent_error+0.01);
assert(abs(unit_tangent_vector(1,2) - sin(angle_in_radians))<sin(angle_in_radians)*percent_error+0.01);
