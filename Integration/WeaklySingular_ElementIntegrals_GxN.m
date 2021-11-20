function GxN = WeaklySingular_ElementIntegrals_GxN(chi, xi_e, mu, ...
                                                   xnodenum, ...
                                                   grx, grw, gtx, gtw)
%%
% chi           coordinates of target point
% xi_e          matrix of element node coord (3x9)
% mu            viscosity
% xnodenum      node number of target point in element
% grx, grw      gauss node coordinates and weights for radial direction r
% gtx, gtw      gauss node coordinates and weights for polar angle theta 
% GxN           weakly-singular element integral
%%
switch xnodenum
    case 1
        zetachi1 = -1; 
        zetachi2 = -1;
        % Triangle 1
        thet_min = 0; 
        thet_max = pi/4; 
        rmax_fn = 1;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = pi/4; 
        thet_max = pi/2; 
        rmax_fn = 2;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2;
    case 2
        zetachi1 = 1; 
        zetachi2 = -1;
        % Triangle 1
        thet_min = pi/2; 
        thet_max = 3*pi/4; 
        rmax_fn = 3;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = 3*pi/4; 
        thet_max = pi; 
        rmax_fn = 4;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2;
    case 3
        zetachi1 = 1; 
        zetachi2 = 1;
        % Triangle 1
        thet_min = pi; 
        thet_max = 5*pi/4; 
        rmax_fn = 5;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = 5*pi/4; 
        thet_max = 3*pi/2; 
        rmax_fn = 6;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2;
    case 4
        zetachi1 = -1; 
        zetachi2 = 1;
        % Triangle 1
        thet_min = 3*pi/2; 
        thet_max = 7*pi/4; 
        rmax_fn = 7;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = 7*pi/4; 
        thet_max = 2*pi; 
        rmax_fn = 8;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2;
    case 5
        zetachi1 = 0; 
        zetachi2 = -1;
        % Triangle 1
        thet_min = 0; 
        thet_max = atan2(2,1); 
        rmax_fn = 9;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = atan2(2,1); 
        thet_max = atan2(2,-1); 
        rmax_fn = 10;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 3
        thet_min = atan2(2,-1); 
        thet_max = pi; 
        rmax_fn = 11;
        GxN_tr3 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2 + GxN_tr3;
    case 6
        zetachi1 = 1; 
        zetachi2 = 0;
        % Triangle 1
        thet_min = pi/2; 
        thet_max = atan2(1,-2); 
        rmax_fn = 12;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = atan2(1,-2); 
        thet_max = atan2(-1,-2) + 2*pi;
        rmax_fn = 13;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 3
        thet_min = atan2(-1,-2) + 2*pi; 
        thet_max = 3*pi/2;
        rmax_fn = 14;
        GxN_tr3 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2 + GxN_tr3;
    case 7
        zetachi1 = 0; 
        zetachi2 = 1;
        % Triangle 1
        thet_min = pi; 
        thet_max = atan2(-2,-1) + 2*pi; 
        rmax_fn = 15;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = atan2(-2,-1) + 2*pi; 
        thet_max = atan2(-2,1) + 2*pi;
        rmax_fn = 16;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 3
        thet_min = atan2(-2,1) + 2*pi; 
        thet_max = 2*pi;
        rmax_fn = 17;
        GxN_tr3 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2 + GxN_tr3;
    case 8
        zetachi1 = -1; 
        zetachi2 = 0;
        % Triangle 1
        thet_min = -pi/2; 
        thet_max = atan2(-1,2); 
        rmax_fn = 18;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = atan2(-1,2); 
        thet_max = atan2(1,2);
        rmax_fn = 19;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 3
        thet_min = atan2(1,2); 
        thet_max = pi/2;
        rmax_fn = 20;
        GxN_tr3 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2 + GxN_tr3;
    case 9
        zetachi1 = 0; 
        zetachi2 = 0;
        % Triangle 1
        thet_min = -pi/4; 
        thet_max = pi/4; 
        rmax_fn = 21;
        GxN_tr1 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 2
        thet_min = pi/4; 
        thet_max = 3*pi/4;
        rmax_fn = 22;
        GxN_tr2 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 3
        thet_min = 3*pi/4; 
        thet_max = 5*pi/4;
        rmax_fn = 23;
        GxN_tr3 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        % Triangle 4
        thet_min = 5*pi/4; 
        thet_max = 7*pi/4;
        rmax_fn = 24;
        GxN_tr4 = WeaklySingularIntegrals_GxN(xi_e, chi, mu, ...
                                              zetachi1, zetachi2, ...
                                              thet_min, thet_max, ...
                                              rmax_fn,...
                                              grx, grw, gtx, gtw);
        GxN = GxN_tr1 + GxN_tr2 + GxN_tr3 + GxN_tr4;
end
