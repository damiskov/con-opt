%{

Generates a contour plot of the Himmelblau function.

f(x)=(x21 +x2 −11)2 +(x1 +x2 −7)2

with the following constraints

f(x)=(x1^2 +x2 -11)^2 +(x1 +x2^2 -7)^2

s.t:
(x1+2)^2 -x2 >= 0
-4x1 +10x2 >= 0

and 

-5 <= x1 <= 5
-5 <= x2 <= 5


Along with extrema.
%}


lb = [-5; -5];
ub = [5; 5];


x1 = linspace(-6, 6, 100);
x2 = linspace(-6, 6, 100);
[X1, X2] = meshgrid(x1, x2);


Z = zeros(size(X1));
for i = 1:size(X1, 1)
    for j = 1:size(X1, 2)
        Z(i, j) = himmelblau([X1(i, j); X2(i, j)]);
    end
end

% Plot the contour
figure;
v = [-10:0.5:0, 0:1:10 10:5:100 100:20:200, 200:50:400];
contour(X1, X2, Z, v);

% Shading the feasible region

yc1 = (x1+2).^2;
yc2 = 0.4.*x1;

hold on
    fill(x1,yc1,[0.7 0.7 0.7],'facealpha',0.5)
    fill([x1 x1(end) x1(1)],[yc2 -5 -5],[0.7 0.7 0.7],'facealpha',0.5)
    % Shading area below lower bound
    fill([x1(1), x1(1), x1(end), x1(end)], [-6, -5, -5, -6], [0.7 0.7 0.7], 'facealpha', 0.5)
    % Shading area above upper bound
    fill([x1(1), x1(1), x1(end), x1(end)], [5, 6, 6, 5], [0.7 0.7 0.7], 'facealpha', 0.5)
    % Shading area to the left of the left bound
    fill([-6, -5, -5, -6], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
    % Shading area to the right of the right bound
    fill([5, 6, 6, 5], [x2(1), x2(1), x2(end), x2(end)], [0.7 0.7 0.7], 'facealpha', 0.5)
    
hold off

xlim([-6 6])
ylim([-6 6])
colorbar
xlabel('x_1');
ylabel('x_2');

% Reading data from csv, and plotting the extrema
data = readmatrix('data/himmelblau_extrema.csv');
% plotting the extrema, if Classification (last column) is maximum=1 then red, minimum=2 then blue, saddle=0 then green
hold on

h_max = plot(NaN, NaN, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', "#A2142F", "DisplayName", "Maximum");
h_min = plot(NaN, NaN, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', "#0072BD", "DisplayName", "Minimum");
h_saddle = plot(NaN, NaN, 'go', 'MarkerSize', 10, 'MarkerFaceColor', "#77AC30", "DisplayName", "Saddle Point");


for i = 1:size(data, 1)
    if data(i, end) == 1
        plot(data(i, 1), data(i, 2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', "#A2142F", "DisplayName", "Maximum"); 
    elseif data(i, end) == 2
        plot(data(i, 1), data(i, 2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', "#0072BD", "DisplayName", "Minimum");
    else
        plot(data(i, 1), data(i, 2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', "#77AC30", "DisplayName", "Saddle Point");
    end
end

hold off

legend('Location', 'northwest');
legend('Himmelblau', '', '', '', '', '', '', 'Maximum', 'Minimum', 'Saddle') % A bit finicky, but it works


% Saving the plot
fp = '/Users/davidmiles-skov/Desktop/Academics/Optimisation/02612 - Constrained Optimisation/Exam Assignment/problem4/figures/himmelblauExtrema.png';
f = gcf;
save_figure(f, fp, 300, 8, 6);
