function ContourLineSpeed(LinePtsX, LinePtsY, Speed)
ax = gca;
HG = hgtransform(ax);
[M,c] = contourf(LinePtsX, LinePtsY, Speed, 100, ...
                 'LineColor', 'none', 'Parent', HG);
% HG.Matrix = makehgtform('xrotate', pi/2);