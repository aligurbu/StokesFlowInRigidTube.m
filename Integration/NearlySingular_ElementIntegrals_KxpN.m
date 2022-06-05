function KxpN = NearlySingular_ElementIntegrals_KxpN(xi_e, chi, ...
                                                     grx, grw, gtx, gtw, ...
                                                     zetap1, zetap2)

% chi           coordinates of target point
% xi_e          matrix of element node coord (3x9)
% grx, grw      gauss node coordinates and weights for radial direction r
% gtx, gtw      gauss node coordinates and weights for polar angle theta 
% zetap1        origin of polar coordinate system
% zetap2                            in intrinsic coordinates(nearest point)
% KxpN          nearly-singular element integral

if (zetap1 == -1 && zetap2 == -1) % Node 1
    % Triangle 1
    thet_min = 0;
    thet_max = pi/4;
    rmax_fn = 1;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = pi/4;
    thet_max = pi/2;
    rmax_fn = 2;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2;
elseif (zetap1 == 1 && zetap2 == -1)  % Node 2
    % Triangle 1
    thet_min = pi/2; 
    thet_max = 3*pi/4; 
    rmax_fn = 3;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 3*pi/4; 
    thet_max = pi; 
    rmax_fn = 4;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2;
elseif (zetap1 == 1 && zetap2 == 1)   % Node 3
    % Triangle 1
    thet_min = pi; 
    thet_max = 5*pi/4; 
    rmax_fn = 5;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 5*pi/4; 
    thet_max = 3*pi/2; 
    rmax_fn = 6;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2;
elseif (zetap1 == -1 && zetap2 == 1)  % Node 4
    % Triangle 1
    thet_min = 3*pi/2; 
    thet_max = 7*pi/4; 
    rmax_fn = 7;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = 7*pi/4; 
    thet_max = 2*pi; 
    rmax_fn = 8;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2;
elseif (zetap2 == -1)  % Generalized Node 5
    % Triangle 1
    thet_min = 0; 
    thet_max = atan2(2,(1-zetap1)); 
    rmax_fn = 9;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(2,(1-zetap1)); 
    thet_max = atan2(2,-(zetap1+1)); 
    rmax_fn = 10;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(2,-(zetap1+1)); 
    thet_max = pi; 
    rmax_fn = 11;
    KxpN_tr3 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2 + KxpN_tr3;
elseif (zetap1 == 1)   % Generalized Node 6
    % Triangle 1
    thet_min = pi/2; 
    thet_max = atan2((1-zetap2),-2); 
    rmax_fn = 12;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2((1-zetap2),-2); 
    thet_max = atan2(-(zetap2+1),-2) + 2*pi;
    rmax_fn = 13;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(-(zetap2+1),-2) + 2*pi; 
    thet_max = 3*pi/2;
    rmax_fn = 14;
    KxpN_tr3 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2 + KxpN_tr3;
elseif (zetap2 == 1)   % Generalized Node 7
    % Triangle 1
    thet_min = pi; 
    thet_max = atan2(-2,-(zetap1+1)) + 2*pi; 
    rmax_fn = 15;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(-2,-(zetap1+1)) + 2*pi; 
    thet_max = atan2(-2,(1-zetap1)) + 2*pi;
    rmax_fn = 16;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2(-2,(1-zetap1)) + 2*pi; 
    thet_max = 2*pi;
    rmax_fn = 17;
    KxpN_tr3 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2 + KxpN_tr3;
elseif (zetap1 == -1)  % Generalized Node 8
    % Triangle 1
    thet_min = -pi/2; 
    thet_max = atan2(-(zetap2+1),2); 
    rmax_fn = 18;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2(-(zetap2+1),2); 
    thet_max = atan2((1-zetap2),2);
    rmax_fn = 19;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2((1-zetap2),2); 
    thet_max = pi/2;
    rmax_fn = 20;
    KxpN_tr3 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2 + KxpN_tr3;
else  % Generalized Node 9
    % Triangle 1
    thet_min = atan2(-(zetap2+1),(1-zetap1)); 
    thet_max = atan2((1-zetap2),(1-zetap1)); 
    rmax_fn = 21;
    KxpN_tr1 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 2
    thet_min = atan2((1-zetap2),(1-zetap1)); 
    thet_max = atan2((1-zetap2),-(zetap1+1));
    rmax_fn = 22;
    KxpN_tr2 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 3
    thet_min = atan2((1-zetap2),-(zetap1+1)); 
    thet_max = atan2(-(zetap2+1),-(zetap1+1))+2*pi;
    rmax_fn = 23;
    KxpN_tr3 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    % Triangle 4
    thet_min = atan2(-(zetap2+1),-(zetap1+1))+2*pi; 
    thet_max = atan2(-(zetap2+1),(1-zetap1))+2*pi;
    rmax_fn = 24;
    KxpN_tr4 = NearlySingularIntegrals_KxpN(xi_e, chi, ...
                                            zetap1, zetap2, ...
                                            thet_min, thet_max, ...
                                            rmax_fn,...
                                            grx, grw, gtx, gtw);
    KxpN = KxpN_tr1 + KxpN_tr2 + KxpN_tr3 + KxpN_tr4;
end
