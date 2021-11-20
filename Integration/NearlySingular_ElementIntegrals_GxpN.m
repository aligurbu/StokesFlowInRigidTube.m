function GxpN = NearlySingular_ElementIntegrals_GxpN(xi_e, chi, mu, ...
                                                     grx, grw, gtx, gtw, ...
                                                     zetap1, zetap2)

% chi           coordinates of target point
% xi_e          matrix of element node coord (3x9)
% mu            viscosity
% grx, grw      gauss node coordinates and weights for radial direction r
% gtx, gtw      gauss node coordinates and weights for polar angle theta 
% zetap1        origin of polar coordinate system
% zetap2                            in intrinsic coordinates(nearest point)
% GxpN          nearly-singular element integral

if (zetap1 == -1 && zetap2 == -1) % Node 1
    % Triangle 1
    thet_min = 0;
    thet_max = pi/4;
    rmax_fn = 1;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = pi/4;
    thet_max = pi/2;
    rmax_fn = 2;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2;
elseif (zetap1 == 1 && zetap2 == -1)  % Node 2
    % Triangle 1
    thet_min = pi/2; 
    thet_max = 3*pi/4; 
    rmax_fn = 3;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 3*pi/4; 
    thet_max = pi; 
    rmax_fn = 4;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2;
elseif (zetap1 == 1 && zetap2 == 1)   % Node 3
    % Triangle 1
    thet_min = pi; 
    thet_max = 5*pi/4; 
    rmax_fn = 5;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 5*pi/4; 
    thet_max = 3*pi/2; 
    rmax_fn = 6;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2; 
elseif (zetap1 == -1 && zetap2 == 1)  % Node 4
    % Triangle 1
    thet_min = 3*pi/2; 
    thet_max = 7*pi/4; 
    rmax_fn = 7;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 7*pi/4; 
    thet_max = 2*pi; 
    rmax_fn = 8;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2;
elseif (zetap2 == -1)  % Generalized Node 5
    % Triangle 1
    thet_min = 0; 
    thet_max = atan2(2,(1-zetap1)); 
    rmax_fn = 9;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(2,(1-zetap1)); 
    thet_max = atan2(2,-(zetap1+1)); 
    rmax_fn = 10;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(2,-(zetap1+1)); 
    thet_max = pi; 
    rmax_fn = 11;
    GxpN_tr3 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2 + GxpN_tr3;
elseif (zetap1 == 1)   % Generalized Node 6
    % Triangle 1
    thet_min = pi/2; 
    thet_max = atan2((1-zetap2),-2); 
    rmax_fn = 12;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2((1-zetap2),-2); 
    thet_max = atan2(-(zetap2+1),-2) + 2*pi;
    rmax_fn = 13;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(-(zetap2+1),-2) + 2*pi; 
    thet_max = 3*pi/2;
    rmax_fn = 14;
    GxpN_tr3 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2 + GxpN_tr3;
elseif (zetap2 == 1)   % Generalized Node 7
    % Triangle 1
    thet_min = pi; 
    thet_max = atan2(-2,-(zetap1+1)) + 2*pi; 
    rmax_fn = 15;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(-2,-(zetap1+1)) + 2*pi; 
    thet_max = atan2(-2,(1-zetap1)) + 2*pi;
    rmax_fn = 16;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(-2,(1-zetap1)) + 2*pi; 
    thet_max = 2*pi;
    rmax_fn = 17;
    GxpN_tr3 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2 + GxpN_tr3;
elseif (zetap1 == -1)  % Generalized Node 8
    % Triangle 1
    thet_min = -pi/2; 
    thet_max = atan2(-(zetap2+1),2); 
    rmax_fn = 18;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(-(zetap2+1),2); 
    thet_max = atan2((1-zetap2),2);
    rmax_fn = 19;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2((1-zetap2),2); 
    thet_max = pi/2;
    rmax_fn = 20;
    GxpN_tr3 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2 + GxpN_tr3;
else  % Generalized Node 9
    % Triangle 1
    thet_min = atan2(-(zetap2+1),(1-zetap1)); 
    thet_max = atan2((1-zetap2),(1-zetap1)); 
    rmax_fn = 21;
    GxpN_tr1 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2((1-zetap2),(1-zetap1)); 
    thet_max = atan2((1-zetap2),-(zetap1+1));
    rmax_fn = 22;
    GxpN_tr2 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2((1-zetap2),-(zetap1+1)); 
    thet_max = atan2(-(zetap2+1),-(zetap1+1))+2*pi;
    rmax_fn = 23;
    GxpN_tr3 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 4
    thet_min = atan2(-(zetap2+1),-(zetap1+1))+2*pi; 
    thet_max = atan2(-(zetap2+1),(1-zetap1))+2*pi;
    rmax_fn = 24;
    GxpN_tr4 = NearlySingularIntegrals_GxpN(xi_e, chi, mu, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    GxpN = GxpN_tr1 + GxpN_tr2 + GxpN_tr3 + GxpN_tr4;
end
