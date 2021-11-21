function ContourCrossSectionSpeed(CrossSectionX, CrossSectionY, Speed, ...
                                                                TranslateX)
ax = gca;
HG = hgtransform(ax);
[M,c] = contourf(CrossSectionX, CrossSectionY, Speed, 100, ...
                 'LineColor', 'none', 'Parent', HG);
HG.Matrix = makehgtform('translate',TranslateX,0,0)* ...
            makehgtform('yrotate', pi/2)* ...
            makehgtform('zrotate', pi/2);