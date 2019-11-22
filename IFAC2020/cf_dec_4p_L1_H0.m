clear;
warning('off','Control:ltiobject:RepeatedChannelNames');
warning('off','MATLAB:nargchk:deprecated');
 
init_simulink
 
Ixx = I(1);
Iyy = I(2);
Izz = I(3);
Dii = 1e-5;
 
 
[Amat, Bmat] = getLinSys(0,0,0,0,0,0,m*g,0,0,0,m,Ixx,Iyy,Izz,Dii);
Cmat = eye(12);
Dmat = zeros(12,4);
G1 = c2d(ss(Amat,Bmat,Cmat,Dmat), controller_h);
G2 = c2d(ss(Amat,Bmat,Cmat,Dmat), controller_h);
G3 = c2d(ss(Amat,Bmat,Cmat,Dmat), controller_h);
G4 = c2d(ss(Amat,Bmat,Cmat,Dmat), controller_h);

n1 = 16; ny1 = 12; nu1 = 4; nz1 = 16; nw1 = 28; %n1 = 12+4 low-pass-filters, nz1 = 12+4 (errors+control weights)
n2 = 16; ny2 = 12; nu2 = 4; nz2 = 16; nw2 = 28;
n3 = 16; ny3 = 12; nu3 = 4; nz3 = 16; nw3 = 28;
n4 = 16; ny4 = 12; nu4 = 4; nz4 = 16; nw4 = 28;
n = n1+n2+n3+n4; ny = ny1+ny2+ny3+ny4; nu = nu1+nu2+nu3+nu4; nz = nz1+nz2+nz3+nz4; nw = nw1+nw2+nw3+nw4;
 
%%
 
Wxy1 = makeweight(5,1,0.1,controller_h);
Wpos = blkdiag(Wxy1,Wxy1,makeweight(5,10,0.01,controller_h));
We = blkdiag(Wpos, 0, 0, 5, 0.5, 0.5, makeweight(5,0.5,0.1,controller_h), 0, 0, 0);
Wu = blkdiag(10, 40, 40, 10);
Wr = blkdiag(1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0);
Wn = diag([2e-2, 2e-2, 7e-2, 5e-2, 5e-2, 0.1, 1e-4, 1e-4, 0.000256, 0.1, 0.1, 0.1]);
Wd = diag([1e-2, 0, 0, 1e-4]);
 
Weps = blkdiag(1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0);
Wxy = makeweight(7,0.5,0.1,controller_h);
Wpos2 = blkdiag(Wxy,Wxy,makeweight(5,10,0.01,controller_h));
We2 = blkdiag(Wpos2, 1, 1, 2, 1, 1, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu2 = blkdiag(10, 40, 40, 5);

Wxy = makeweight(7,0.3,0.1,controller_h);
Wpos3 = blkdiag(Wxy,Wxy,makeweight(10,20,0.01,controller_h));
We3 = blkdiag(Wpos3, 1, 1, 2, 1, 1, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu3 = blkdiag(10, 40, 40, 5);

Wpos4 = blkdiag(Wxy,Wxy,makeweight(10,20,0.01,controller_h));
We4 = blkdiag(Wpos4, 1.5, 1.5, 2, 1, 1, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu4 = blkdiag(10, 40, 40, 5);
 
systemnames = 'G1 We Wu Wr Wn Wd';
inputvar = '[r1{12};n1{12};d1{4};u1{4}]';
outputvar = '[We; Wu; Wr-G1-Wn; G1]'; %Ze, Zu, r-(x+n), x
input_to_G1 = '[u1+Wd]';
input_to_We = '[Wr-G1]';
input_to_Wu = '[u1]';
input_to_Wr = '[r1]';
input_to_Wn = '[n1]';
input_to_Wd = '[d1]';
sysoutname = 'P1'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G2 We2 Wu2 Wr Wn Wd Weps';
inputvar = '[r2{12};n2{12};d2{4};u2{4};eps{12}]';
outputvar = '[We2; Wu2; Wr-G2-Wn+Weps; G2]'; %Ze, Zu, r-(x+n)+eps
input_to_G2 = '[u2+Wd]';
input_to_We2 = '[Wr-G2]';
input_to_Wu2 = '[u2]';
input_to_Wr = '[r2]';
input_to_Wn = '[n2]';
input_to_Wd = '[d2]';
input_to_Weps = '[eps]';
sysoutname = 'P2'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G3 We3 Wu3 Wr Wn Wd Weps';
inputvar = '[r3{12};n3{12};d3{4};u3{4};eps{12}]';
outputvar = '[We3; Wu3; Wr-G3-Wn+Weps; G3]'; %Ze, Zu, r-(x+n)+eps
input_to_G3 = '[u3+Wd]';
input_to_We3 = '[Wr-G3]';
input_to_Wu3 = '[u3]';
input_to_Wr = '[r3]';
input_to_Wn = '[n3]';
input_to_Wd = '[d3]';
input_to_Weps = '[eps]';
sysoutname = 'P3'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G4 We4 Wu4 Wr Wn Wd Weps';
inputvar = '[r4{12};n4{12};d4{4};u4{4};eps{12}]';
outputvar = '[We4; Wu4; Wr-G4-Wn+Weps]'; %Ze, Zu, r-(x+n)+eps
input_to_G4 = '[u4+Wd]';
input_to_We4 = '[Wr-G4]';
input_to_Wu4 = '[u4]';
input_to_Wr = '[r4]';
input_to_Wn = '[n4]';
input_to_Wd = '[d4]';
input_to_Weps = '[eps]';
sysoutname = 'P4'; cleanupsysic = 'yes'; sysic

systemnames = 'P1 P2 P3 P4';
inputvar = '[r{48};n{48};d{16};u{16}]';
outputvar = '[P1(1:16); P2(1:16); P3(1:16); P4(1:16); P1(17:28); P2(17:28); P3(17:28); P4(17:28)]';
input_to_P1 = '[r(1:12); n(1:12); d(1:4); u(1:4)]';
input_to_P2 = '[r(13:24); n(13:24); d(5:8); u(5:8); P1(29:40)]';
input_to_P3 = '[r(25:36); n(25:36); d(9:12); u(9:12); P2(29:40)]';
input_to_P4 = '[r(37:48); n(37:48); d(13:16); u(13:16); P3(29:40)]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic
 
A1 = P.a;
Bw1 = P.b(:,1:nw);
Bu1 = P.b(:,nw+1:end);
Cz1 = P.c(1:nz,:);
Cy1 = P.c(nz+1:end,:);
Dzw1 = P.d(1:nz,1:nw);
Dzu1 = P.d(1:nz,nw+1:end);
Dyw1 = P.d(nz+1:end,1:nw);
 
%%
Wxy1 = makeweight(12,5,0.1,controller_h);
Wpos = blkdiag(Wxy1,Wxy1,makeweight(10,10,0.01,controller_h));
We = blkdiag(Wpos, 5, 5, 1, 1, 1, makeweight(5,0.5,0.1,controller_h), 0, 0, 0);
Wu = blkdiag(10, 40, 40, 10);
Wr = blkdiag(1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0);
Wd = diag([1e-2, 0, 0, 1e-4]);
 
Wxy = makeweight(15,5,0.1,controller_h);
Wpos2 = blkdiag(Wxy,Wxy,makeweight(10,10,0.01,controller_h));
We2 = blkdiag(Wpos2, 5, 5, 1, 3, 3, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu2 = blkdiag(10, 40, 40, 5);
Weps = blkdiag(1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0);
 
Wpos3 = blkdiag(Wxy,Wxy,makeweight(10,10,0.01,controller_h));
We3 = blkdiag(Wpos3, 5, 5, 1, 3, 3, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu3 = blkdiag(10, 40, 40, 5);

Wpos4 = blkdiag(Wxy,Wxy,makeweight(10,10,0.01,controller_h));
We4 = blkdiag(Wpos4, 5, 5, 1, 5, 5, makeweight(5,1.0,0.1,controller_h), 0, 0, 0);
Wu4 = blkdiag(10, 40, 40, 5);
 
systemnames = 'G1 We Wu Wr Wn Wd';
inputvar = '[r1{12};n1{12};d1{4};u1{4}]';
outputvar = '[We; Wu; Wr-G1-Wn; G1]'; %Ze, Zu, r-(x+n), x
input_to_G1 = '[u1+Wd]';
input_to_We = '[Wr-G1]';
input_to_Wu = '[u1]';
input_to_Wr = '[r1]';
input_to_Wn = '[n1]';
input_to_Wd = '[d1]';
sysoutname = 'P1'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G2 We2 Wu2 Wr Wn Wd Weps';
inputvar = '[r2{12};n2{12};d2{4};u2{4};eps{12}]';
outputvar = '[We2; Wu2; Wr-G2-Wn+Weps; G2]'; %Ze, Zu, r-(x+n)+eps
input_to_G2 = '[u2+Wd]';
input_to_We2 = '[Wr-G2]';
input_to_Wu2 = '[u2]';
input_to_Wr = '[r2]';
input_to_Wn = '[n2]';
input_to_Wd = '[d2]';
input_to_Weps = '[eps]';
sysoutname = 'P2'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G3 We3 Wu3 Wr Wn Wd Weps';
inputvar = '[r3{12};n3{12};d3{4};u3{4};eps{12}]';
outputvar = '[We3; Wu3; Wr-G3-Wn+Weps; G3]'; %Ze, Zu, r-(x+n)+eps
input_to_G3 = '[u3+Wd]';
input_to_We3 = '[Wr-G3]';
input_to_Wu3 = '[u3]';
input_to_Wr = '[r3]';
input_to_Wn = '[n3]';
input_to_Wd = '[d3]';
input_to_Weps = '[eps]';
sysoutname = 'P3'; cleanupsysic = 'yes'; sysic
 
systemnames = 'G4 We4 Wu4 Wr Wn Wd Weps';
inputvar = '[r4{12};n4{12};d4{4};u4{4};eps{12}]';
outputvar = '[We4; Wu4; Wr-G4-Wn+Weps]'; %Ze, Zu, r-(x+n)+eps
input_to_G4 = '[u4+Wd]';
input_to_We4 = '[Wr-G4]';
input_to_Wu4 = '[u4]';
input_to_Wr = '[r4]';
input_to_Wn = '[n4]';
input_to_Wd = '[d4]';
input_to_Weps = '[eps]';
sysoutname = 'P4'; cleanupsysic = 'yes'; sysic

systemnames = 'P1 P2 P3 P4';
inputvar = '[r{48};n{48};d{16};u{16}]';
outputvar = '[P1(1:16); P2(1:16); P3(1:16); P4(1:16); P1(17:28); P2(17:28); P3(17:28); P4(17:28)]';
input_to_P1 = '[r(1:12); n(1:12); d(1:4); u(1:4)]';
input_to_P2 = '[r(13:24); n(13:24); d(5:8); u(5:8); P1(29:40)]';
input_to_P3 = '[r(25:36); n(25:36); d(9:12); u(9:12); P2(29:40)]';
input_to_P4 = '[r(37:48); n(37:48); d(13:16); u(13:16); P3(29:40)]';
sysoutname = 'P'; cleanupsysic = 'yes'; sysic
 
A2 = P.a;
Bw2 = P.b(:,1:nw);
Bu2 = P.b(:,nw+1:end);
Cz2 = P.c(1:nz,:);
Cy2 = P.c(nz+1:end,:);
Dzw2 = P.d(1:nz,1:nw);
Dzu2 = P.d(1:nz,nw+1:end);
Dyw2 = P.d(nz+1:end,1:nw);
%% 
 
Ey_0 = zeros(ny,0);
Ey_1 = [eye(ny1);zeros(ny2+ny3+ny4,ny1)];
Ey_2 = [eye(ny1+ny2);zeros(ny3+ny4,ny1+ny2)];
Ey_3 = [eye(ny1+ny2+ny3);zeros(ny4,ny1+ny2+ny3)];
Ey_4 = eye(ny);

Eu_0 = zeros(nu,0);
Eu_1 = [eye(nu1);zeros(nu2+nu3+nu4,nu1)];
Eu_2 = [eye(nu1+nu2);zeros(nu3+nu4,nu1+nu2)];
Eu_3 = [eye(nu1+nu2+nu3);zeros(nu4,nu1+nu2+nu3)];
Eu_4 = eye(nu);

E_0 = zeros(n,0);
E_1 = [eye(n1);zeros(n2+n3+n4,n1)];
E_2 = [eye(n1+n2);zeros(n3+n4,n1+n2)];
E_3 = [eye(n1+n2+n3);zeros(n4,n1+n2+n3)];
E_4 = eye(n);

Eybar_0 = null(Ey_0'); Eybar_1 = null(Ey_1'); Eybar_2 = null(Ey_2'); Eybar_3 = null(Ey_3'); Eybar_4 = null(Ey_4'); 
Eubar_0 = null(Eu_0'); Eubar_1 = null(Eu_1'); Eubar_2 = null(Eu_2'); Eubar_3 = null(Eu_3'); Eubar_4 = null(Eu_4'); 
Ebar_0 = null(E_0'); Ebar_1 = null(E_1'); Ebar_2 = null(E_2'); Ebar_3 = null(E_3'); Ebar_4 = null(E_4'); 

Ny_0_1 = null([Ey_0'*Cy1 Ey_0'*Dyw1]);
Ny_0_2 = null([Ey_0'*Cy2 Ey_0'*Dyw2]);
Ny_1_1 = null([Ey_1'*Cy1 Ey_1'*Dyw1]);
Ny_1_2 = null([Ey_1'*Cy2 Ey_1'*Dyw2]);
Ny_2_1 = null([Ey_2'*Cy1 Ey_2'*Dyw1]);
Ny_2_2 = null([Ey_2'*Cy2 Ey_2'*Dyw2]);
Ny_3_1 = null([Ey_3'*Cy1 Ey_3'*Dyw1]);
Ny_3_2 = null([Ey_3'*Cy2 Ey_3'*Dyw2]);
Ny_4_1 = null([Ey_4'*Cy1 Ey_4'*Dyw1]);
Ny_4_2 = null([Ey_4'*Cy2 Ey_4'*Dyw2]);
A11_0_1 = E_0'*A1*E_0; A21_0_1 = Ebar_0'*A1*E_0; A22_0_1 = Ebar_0'*A1*Ebar_0;
A11_0_2 = E_0'*A2*E_0; A21_0_2 = Ebar_0'*A2*E_0; A22_0_2 = Ebar_0'*A2*Ebar_0;
A11_1_1 = E_1'*A1*E_1; A21_1_1 = Ebar_1'*A1*E_1; A22_1_1 = Ebar_1'*A1*Ebar_1;
A11_1_2 = E_1'*A2*E_1; A21_1_2 = Ebar_1'*A2*E_1; A22_1_2 = Ebar_1'*A2*Ebar_1;
A11_2_1 = E_2'*A1*E_2; A21_2_1 = Ebar_2'*A1*E_2; A22_2_1 = Ebar_2'*A1*Ebar_2;
A11_2_2 = E_2'*A2*E_2; A21_2_2 = Ebar_2'*A2*E_2; A22_2_2 = Ebar_2'*A2*Ebar_2;
A11_3_1 = E_3'*A1*E_3; A21_3_1 = Ebar_3'*A1*E_3; A22_3_1 = Ebar_3'*A1*Ebar_3;
A11_3_2 = E_3'*A2*E_3; A21_3_2 = Ebar_3'*A2*E_3; A22_3_2 = Ebar_3'*A2*Ebar_3;
A11_4_1 = E_4'*A1*E_4; A21_4_1 = Ebar_4'*A1*E_4; A22_4_1 = Ebar_4'*A1*Ebar_4;
A11_4_2 = E_4'*A2*E_4; A21_4_2 = Ebar_4'*A2*E_4; A22_4_2 = Ebar_4'*A2*Ebar_4;
Czf1 = Cz1; Cyf1 = Cy1; Dzwf1 = Dzw1; Dzuf1 = Dzu1; Dywf1 = Dyw1;
Czf2 = Cz2; Cyf2 = Cy2; Dzwf2 = Dzw2; Dzuf2 = Dzu2; Dywf2 = Dyw2;

%% 
 
lb = 20;
ub = 60;
tol = 1e-2;
cvx_flag = 1;
eta = 1e-6;

% while(cvx_flag)
% gamm = (ub+lb)/2
gamm = 45;
	cvx_begin sdp quiet 
		cvx_solver Mosek 
		variable Za_1_1(n1,n1) symmetric 
		variable Za_1_2(n1,n1) symmetric 
		variable Za_2_1(n1+n2,n1+n2) symmetric 
		variable Za_2_2(n1+n2,n1+n2) symmetric 
		variable Za_3_1(n1+n2+n3,n1+n2+n3) symmetric 
		variable Za_3_2(n1+n2+n3,n1+n2+n3) symmetric 
		variable Za_4_1(n1+n2+n3+n4,n1+n2+n3+n4) symmetric 
		variable Za_4_2(n1+n2+n3+n4,n1+n2+n3+n4) symmetric 
		variable Zb_1_1(n1,n2+n3+n4) 
		variable Zb_1_2(n1,n2+n3+n4) 
		variable Zb_2_1(n1+n2,n3+n4) 
		variable Zb_2_2(n1+n2,n3+n4) 
		variable Zb_3_1(n1+n2+n3,n4) 
		variable Zb_3_2(n1+n2+n3,n4) 
		variable Zc_0_1(n1+n2+n3+n4,n1+n2+n3+n4) symmetric 
		variable Zc_0_2(n1+n2+n3+n4,n1+n2+n3+n4) symmetric 
		variable Zc_1_1(n2+n3+n4,n2+n3+n4) symmetric 
		variable Zc_1_2(n2+n3+n4,n2+n3+n4) symmetric 
		variable Zc_2_1(n3+n4,n3+n4) symmetric 
		variable Zc_2_2(n3+n4,n3+n4) symmetric 
		variable Zc_3_1(n4,n4) symmetric 
		variable Zc_3_2(n4,n4) symmetric 

        Za_1_1 >= eta*eye(n1);
        Za_1_2 >= eta*eye(n1);
        Za_2_1 >= eta*eye(n1+n2);
        Za_2_2 >= eta*eye(n1+n2);
        Za_3_1 >= eta*eye(n1+n2+n3);
        Za_3_2 >= eta*eye(n1+n2+n3);
        Za_4_1 >= eta*eye(n1+n2+n3+n4);
        Za_4_2 >= eta*eye(n1+n2+n3+n4);
        Zc_0_1 >= eta*eye(n1+n2+n3+n4);
        Zc_0_2 >= eta*eye(n1+n2+n3+n4);
        Zc_1_1 >= eta*eye(n2+n3+n4);
        Zc_1_2 >= eta*eye(n2+n3+n4);
        Zc_2_1 >= eta*eye(n3+n4);
        Zc_2_2 >= eta*eye(n3+n4);
        Zc_3_1 >= eta*eye(n4);
        Zc_3_2 >= eta*eye(n4);

		Cz1 = Czf1/gamm; Cy1 = Cyf1/gamm; Dzw1 = Dzwf1/gamm; Dzu1 = Dzuf1/gamm; Dyw1 = Dywf1/gamm; 
		Cz2 = Czf2/gamm; Cy2 = Cyf2/gamm; Dzw2 = Dzwf2/gamm; Dzu2 = Dzuf2/gamm; Dyw2 = Dywf2/gamm; 

		Nu_0_1 = null([Eubar_0'*Bu1' Eubar_0'*Dzu1']);
		Nu_0_2 = null([Eubar_0'*Bu2' Eubar_0'*Dzu2']);
		Nu_1_1 = null([Eubar_1'*Bu1' Eubar_1'*Dzu1']);
		Nu_1_2 = null([Eubar_1'*Bu2' Eubar_1'*Dzu2']);
		Nu_2_1 = null([Eubar_2'*Bu1' Eubar_2'*Dzu1']);
		Nu_2_2 = null([Eubar_2'*Bu2' Eubar_2'*Dzu2']);
		Nu_3_1 = null([Eubar_3'*Bu1' Eubar_3'*Dzu1']);
		Nu_3_2 = null([Eubar_3'*Bu2' Eubar_3'*Dzu2']);
		Nu_4_1 = null([Eubar_4'*Bu1' Eubar_4'*Dzu1']);
		Nu_4_2 = null([Eubar_4'*Bu2' Eubar_4'*Dzu2']);

		G_0_1 = blkdiag(Nu_0_1, Ny_0_1);
		G_0_2 = blkdiag(Nu_0_2, Ny_0_2);
		G_1_1 = blkdiag(Nu_1_1, Ny_1_1);
		G_1_2 = blkdiag(Nu_1_2, Ny_1_2);
		G_2_1 = blkdiag(Nu_2_1, Ny_2_1);
		G_2_2 = blkdiag(Nu_2_2, Ny_2_2);
		G_3_1 = blkdiag(Nu_3_1, Ny_3_1);
		G_3_2 = blkdiag(Nu_3_2, Ny_3_2);
		G_4_1 = blkdiag(Nu_4_1, Ny_4_1);
		G_4_2 = blkdiag(Nu_4_2, Ny_4_2);

		Zl_0_1 = Zc_0_1; Zu_0_1 = eye(n); 
		Zl_0_2 = Zc_0_2; Zu_0_2 = eye(n); 
		Zl_1_1 = [eye(n1) zeros(n1,n2+n3+n4); -Zb_1_1' Zc_1_1]; Zu_1_1 = [Za_1_1 Zb_1_1; zeros(n2+n3+n4,n1) eye(n2+n3+n4)];
		Zl_1_2 = [eye(n1) zeros(n1,n2+n3+n4); -Zb_1_2' Zc_1_2]; Zu_1_2 = [Za_1_2 Zb_1_2; zeros(n2+n3+n4,n1) eye(n2+n3+n4)];
		Zl_2_1 = [eye(n1+n2) zeros(n1+n2,n3+n4); -Zb_2_1' Zc_2_1]; Zu_2_1 = [Za_2_1 Zb_2_1; zeros(n3+n4,n1+n2) eye(n3+n4)];
		Zl_2_2 = [eye(n1+n2) zeros(n1+n2,n3+n4); -Zb_2_2' Zc_2_2]; Zu_2_2 = [Za_2_2 Zb_2_2; zeros(n3+n4,n1+n2) eye(n3+n4)];
		Zl_3_1 = [eye(n1+n2+n3) zeros(n1+n2+n3,n4); -Zb_3_1' Zc_3_1]; Zu_3_1 = [Za_3_1 Zb_3_1; zeros(n4,n1+n2+n3) eye(n4)];
		Zl_3_2 = [eye(n1+n2+n3) zeros(n1+n2+n3,n4); -Zb_3_2' Zc_3_2]; Zu_3_2 = [Za_3_2 Zb_3_2; zeros(n4,n1+n2+n3) eye(n4)];
		Zl_4_1 = eye(n); Zu_4_1 = Za_4_1; 
		Zl_4_2 = eye(n); Zu_4_2 = Za_4_2; 

		S_0_1 = Zc_0_1;
		S_0_2 = Zc_0_2;
		S_1_1 = blkdiag(Za_1_1, Zc_1_1);
		S_1_2 = blkdiag(Za_1_2, Zc_1_2);
		S_2_1 = blkdiag(Za_2_1, Zc_2_1);
		S_2_2 = blkdiag(Za_2_2, Zc_2_2);
		S_3_1 = blkdiag(Za_3_1, Zc_3_1);
		S_3_2 = blkdiag(Za_3_2, Zc_3_2);
		S_4_1 = Za_4_1;
		S_4_2 = Za_4_2;

		ZAZ_0_11 = A22_0_1*Zc_0_1; ZAZ_0_11_t = Zc_0_1*A22_0_1';
		ZAZ_0_12 = A22_0_2*Zc_0_1; ZAZ_0_12_t = Zc_0_1*A22_0_2';
		ZAZ_0_21 = A22_0_1*Zc_0_2; ZAZ_0_21_t = Zc_0_2*A22_0_1';
		ZAZ_0_22 = A22_0_2*Zc_0_2; ZAZ_0_22_t = Zc_0_2*A22_0_2';
		ZAZ_1_11 = [Za_1_1*A11_1_1 zeros(n1,n2+n3+n4); Zb_1_1'*A11_1_1+A21_1_1-A22_1_1*Zb_1_1' A22_1_1*Zc_1_1]; ZAZ_1_11_t = [A11_1_1'*Za_1_1 A11_1_1'*Zb_1_1+A21_1_1'-Zb_1_1*A22_1_1'; zeros(n2+n3+n4,n1) Zc_1_1*A22_1_1'];
		ZAZ_1_12 = [Za_1_2*A11_1_2 zeros(n1,n2+n3+n4); Zb_1_2'*A11_1_2+A21_1_2-A22_1_2*Zb_1_1' A22_1_2*Zc_1_1]; ZAZ_1_12_t = [A11_1_2'*Za_1_2 A11_1_2'*Zb_1_2+A21_1_2'-Zb_1_1*A22_1_2'; zeros(n2+n3+n4,n1) Zc_1_1*A22_1_2'];
		ZAZ_1_21 = [Za_1_1*A11_1_1 zeros(n1,n2+n3+n4); Zb_1_1'*A11_1_1+A21_1_1-A22_1_1*Zb_1_2' A22_1_1*Zc_1_2]; ZAZ_1_21_t = [A11_1_1'*Za_1_1 A11_1_1'*Zb_1_1+A21_1_1'-Zb_1_2*A22_1_1'; zeros(n2+n3+n4,n1) Zc_1_2*A22_1_1'];
		ZAZ_1_22 = [Za_1_2*A11_1_2 zeros(n1,n2+n3+n4); Zb_1_2'*A11_1_2+A21_1_2-A22_1_2*Zb_1_2' A22_1_2*Zc_1_2]; ZAZ_1_22_t = [A11_1_2'*Za_1_2 A11_1_2'*Zb_1_2+A21_1_2'-Zb_1_2*A22_1_2'; zeros(n2+n3+n4,n1) Zc_1_2*A22_1_2'];
		ZAZ_2_11 = [Za_2_1*A11_2_1 zeros(n1+n2,n3+n4); Zb_2_1'*A11_2_1+A21_2_1-A22_2_1*Zb_2_1' A22_2_1*Zc_2_1]; ZAZ_2_11_t = [A11_2_1'*Za_2_1 A11_2_1'*Zb_2_1+A21_2_1'-Zb_2_1*A22_2_1'; zeros(n3+n4,n1+n2) Zc_2_1*A22_2_1'];
		ZAZ_2_12 = [Za_2_2*A11_2_2 zeros(n1+n2,n3+n4); Zb_2_2'*A11_2_2+A21_2_2-A22_2_2*Zb_2_1' A22_2_2*Zc_2_1]; ZAZ_2_12_t = [A11_2_2'*Za_2_2 A11_2_2'*Zb_2_2+A21_2_2'-Zb_2_1*A22_2_2'; zeros(n3+n4,n1+n2) Zc_2_1*A22_2_2'];
		ZAZ_2_21 = [Za_2_1*A11_2_1 zeros(n1+n2,n3+n4); Zb_2_1'*A11_2_1+A21_2_1-A22_2_1*Zb_2_2' A22_2_1*Zc_2_2]; ZAZ_2_21_t = [A11_2_1'*Za_2_1 A11_2_1'*Zb_2_1+A21_2_1'-Zb_2_2*A22_2_1'; zeros(n3+n4,n1+n2) Zc_2_2*A22_2_1'];
		ZAZ_2_22 = [Za_2_2*A11_2_2 zeros(n1+n2,n3+n4); Zb_2_2'*A11_2_2+A21_2_2-A22_2_2*Zb_2_2' A22_2_2*Zc_2_2]; ZAZ_2_22_t = [A11_2_2'*Za_2_2 A11_2_2'*Zb_2_2+A21_2_2'-Zb_2_2*A22_2_2'; zeros(n3+n4,n1+n2) Zc_2_2*A22_2_2'];
		ZAZ_3_11 = [Za_3_1*A11_3_1 zeros(n1+n2+n3,n4); Zb_3_1'*A11_3_1+A21_3_1-A22_3_1*Zb_3_1' A22_3_1*Zc_3_1]; ZAZ_3_11_t = [A11_3_1'*Za_3_1 A11_3_1'*Zb_3_1+A21_3_1'-Zb_3_1*A22_3_1'; zeros(n4,n1+n2+n3) Zc_3_1*A22_3_1'];
		ZAZ_3_12 = [Za_3_2*A11_3_2 zeros(n1+n2+n3,n4); Zb_3_2'*A11_3_2+A21_3_2-A22_3_2*Zb_3_1' A22_3_2*Zc_3_1]; ZAZ_3_12_t = [A11_3_2'*Za_3_2 A11_3_2'*Zb_3_2+A21_3_2'-Zb_3_1*A22_3_2'; zeros(n4,n1+n2+n3) Zc_3_1*A22_3_2'];
		ZAZ_3_21 = [Za_3_1*A11_3_1 zeros(n1+n2+n3,n4); Zb_3_1'*A11_3_1+A21_3_1-A22_3_1*Zb_3_2' A22_3_1*Zc_3_2]; ZAZ_3_21_t = [A11_3_1'*Za_3_1 A11_3_1'*Zb_3_1+A21_3_1'-Zb_3_2*A22_3_1'; zeros(n4,n1+n2+n3) Zc_3_2*A22_3_1'];
		ZAZ_3_22 = [Za_3_2*A11_3_2 zeros(n1+n2+n3,n4); Zb_3_2'*A11_3_2+A21_3_2-A22_3_2*Zb_3_2' A22_3_2*Zc_3_2]; ZAZ_3_22_t = [A11_3_2'*Za_3_2 A11_3_2'*Zb_3_2+A21_3_2'-Zb_3_2*A22_3_2'; zeros(n4,n1+n2+n3) Zc_3_2*A22_3_2'];
		ZAZ_4_11 = Za_4_1*A11_4_1; ZAZ_4_11_t = A11_4_1'*Za_4_1;
		ZAZ_4_12 = Za_4_2*A11_4_2; ZAZ_4_12_t = A11_4_2'*Za_4_2;
		ZAZ_4_21 = Za_4_1*A11_4_1; ZAZ_4_21_t = A11_4_1'*Za_4_1;
		ZAZ_4_22 = Za_4_2*A11_4_2; ZAZ_4_22_t = A11_4_2'*Za_4_2;

		M0_11 = [S_0_1 zeros(n,nz) ZAZ_0_11 Zu_0_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_0_1 Dzw1; ZAZ_0_11_t Zl_0_1'*Cz1' S_0_1 zeros(n,nw);Bw1'*Zu_0_1 Dzw1' zeros(nw,n) eye(nw)];
		M0_12 = [S_0_2 zeros(n,nz) ZAZ_0_12 Zu_0_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_0_1 Dzw2; ZAZ_0_12_t Zl_0_1'*Cz2' S_0_1 zeros(n,nw);Bw2'*Zu_0_2 Dzw2' zeros(nw,n) eye(nw)];
		M0_21 = [S_0_1 zeros(n,nz) ZAZ_0_21 Zu_0_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_0_2 Dzw1; ZAZ_0_21_t Zl_0_2'*Cz1' S_0_2 zeros(n,nw);Bw1'*Zu_0_1 Dzw1' zeros(nw,n) eye(nw)];
		M0_22 = [S_0_2 zeros(n,nz) ZAZ_0_22 Zu_0_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_0_2 Dzw2; ZAZ_0_22_t Zl_0_2'*Cz2' S_0_2 zeros(n,nw);Bw2'*Zu_0_2 Dzw2' zeros(nw,n) eye(nw)];
		M1_11 = [S_1_1 zeros(n,nz) ZAZ_1_11 Zu_1_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_1_1 Dzw1; ZAZ_1_11_t Zl_1_1'*Cz1' S_1_1 zeros(n,nw);Bw1'*Zu_1_1 Dzw1' zeros(nw,n) eye(nw)];
		M1_12 = [S_1_2 zeros(n,nz) ZAZ_1_12 Zu_1_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_1_1 Dzw2; ZAZ_1_12_t Zl_1_1'*Cz2' S_1_1 zeros(n,nw);Bw2'*Zu_1_2 Dzw2' zeros(nw,n) eye(nw)];
		M1_21 = [S_1_1 zeros(n,nz) ZAZ_1_21 Zu_1_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_1_2 Dzw1; ZAZ_1_21_t Zl_1_2'*Cz1' S_1_2 zeros(n,nw);Bw1'*Zu_1_1 Dzw1' zeros(nw,n) eye(nw)];
		M1_22 = [S_1_2 zeros(n,nz) ZAZ_1_22 Zu_1_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_1_2 Dzw2; ZAZ_1_22_t Zl_1_2'*Cz2' S_1_2 zeros(n,nw);Bw2'*Zu_1_2 Dzw2' zeros(nw,n) eye(nw)];
		M2_11 = [S_2_1 zeros(n,nz) ZAZ_2_11 Zu_2_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_2_1 Dzw1; ZAZ_2_11_t Zl_2_1'*Cz1' S_2_1 zeros(n,nw);Bw1'*Zu_2_1 Dzw1' zeros(nw,n) eye(nw)];
		M2_12 = [S_2_2 zeros(n,nz) ZAZ_2_12 Zu_2_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_2_1 Dzw2; ZAZ_2_12_t Zl_2_1'*Cz2' S_2_1 zeros(n,nw);Bw2'*Zu_2_2 Dzw2' zeros(nw,n) eye(nw)];
		M2_21 = [S_2_1 zeros(n,nz) ZAZ_2_21 Zu_2_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_2_2 Dzw1; ZAZ_2_21_t Zl_2_2'*Cz1' S_2_2 zeros(n,nw);Bw1'*Zu_2_1 Dzw1' zeros(nw,n) eye(nw)];
		M2_22 = [S_2_2 zeros(n,nz) ZAZ_2_22 Zu_2_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_2_2 Dzw2; ZAZ_2_22_t Zl_2_2'*Cz2' S_2_2 zeros(n,nw);Bw2'*Zu_2_2 Dzw2' zeros(nw,n) eye(nw)];
		M3_11 = [S_3_1 zeros(n,nz) ZAZ_3_11 Zu_3_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_3_1 Dzw1; ZAZ_3_11_t Zl_3_1'*Cz1' S_3_1 zeros(n,nw);Bw1'*Zu_3_1 Dzw1' zeros(nw,n) eye(nw)];
		M3_12 = [S_3_2 zeros(n,nz) ZAZ_3_12 Zu_3_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_3_1 Dzw2; ZAZ_3_12_t Zl_3_1'*Cz2' S_3_1 zeros(n,nw);Bw2'*Zu_3_2 Dzw2' zeros(nw,n) eye(nw)];
		M3_21 = [S_3_1 zeros(n,nz) ZAZ_3_21 Zu_3_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_3_2 Dzw1; ZAZ_3_21_t Zl_3_2'*Cz1' S_3_2 zeros(n,nw);Bw1'*Zu_3_1 Dzw1' zeros(nw,n) eye(nw)];
		M3_22 = [S_3_2 zeros(n,nz) ZAZ_3_22 Zu_3_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_3_2 Dzw2; ZAZ_3_22_t Zl_3_2'*Cz2' S_3_2 zeros(n,nw);Bw2'*Zu_3_2 Dzw2' zeros(nw,n) eye(nw)];
		M4_11 = [S_4_1 zeros(n,nz) ZAZ_4_11 Zu_4_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_4_1 Dzw1; ZAZ_4_11_t Zl_4_1'*Cz1' S_4_1 zeros(n,nw);Bw1'*Zu_4_1 Dzw1' zeros(nw,n) eye(nw)];
		M4_12 = [S_4_2 zeros(n,nz) ZAZ_4_12 Zu_4_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_4_1 Dzw2; ZAZ_4_12_t Zl_4_1'*Cz2' S_4_1 zeros(n,nw);Bw2'*Zu_4_2 Dzw2' zeros(nw,n) eye(nw)];
		M4_21 = [S_4_1 zeros(n,nz) ZAZ_4_21 Zu_4_1'*Bw1;zeros(nz,n) eye(nz) Cz1*Zl_4_2 Dzw1; ZAZ_4_21_t Zl_4_2'*Cz1' S_4_2 zeros(n,nw);Bw1'*Zu_4_1 Dzw1' zeros(nw,n) eye(nw)];
		M4_22 = [S_4_2 zeros(n,nz) ZAZ_4_22 Zu_4_2'*Bw2;zeros(nz,n) eye(nz) Cz2*Zl_4_2 Dzw2; ZAZ_4_22_t Zl_4_2'*Cz2' S_4_2 zeros(n,nw);Bw2'*Zu_4_2 Dzw2' zeros(nw,n) eye(nw)];

		G_0_1'*M0_11*G_0_1 >= eta*eye(288);
		G_0_2'*M0_12*G_0_2 >= eta*eye(288);
		G_0_1'*M0_21*G_0_1 >= eta*eye(288);
		G_0_2'*M0_22*G_0_2 >= eta*eye(288);
		G_1_1'*M1_11*G_1_1 >= eta*eye(280);
		G_1_2'*M1_12*G_1_2 >= eta*eye(280);
		G_1_1'*M1_21*G_1_1 >= eta*eye(280);
		G_1_2'*M1_22*G_1_2 >= eta*eye(280);
		G_2_1'*M2_11*G_2_1 >= eta*eye(272);
		G_2_2'*M2_12*G_2_2 >= eta*eye(272);
		G_2_1'*M2_21*G_2_1 >= eta*eye(272);
		G_2_2'*M2_22*G_2_2 >= eta*eye(272);
		G_3_1'*M3_11*G_3_1 >= eta*eye(264);
		G_3_2'*M3_12*G_3_2 >= eta*eye(264);
		G_3_1'*M3_21*G_3_1 >= eta*eye(264);
		G_3_2'*M3_22*G_3_2 >= eta*eye(264);
		G_4_1'*M4_11*G_4_1 >= eta*eye(256);
		G_4_2'*M4_12*G_4_2 >= eta*eye(256);
		G_4_1'*M4_21*G_4_1 >= eta*eye(256);
		G_4_2'*M4_22*G_4_2 >= eta*eye(256);

		S_2_1_1 = blkdiag(Za_1_1,eye(n2),Zc_2_1)+[zeros(n1) Zb_1_1;zeros(n2+n3+n4,n1) zeros(n2+n3+n4)]-[zeros(n1+n2) Zb_2_1;zeros(n3+n4,n1+n2) zeros(n3+n4)];
		S_2_1_2 = blkdiag(Za_1_2,eye(n2),Zc_2_2)+[zeros(n1) Zb_1_2;zeros(n2+n3+n4,n1) zeros(n2+n3+n4)]-[zeros(n1+n2) Zb_2_2;zeros(n3+n4,n1+n2) zeros(n3+n4)];
		S_3_2_1 = blkdiag(Za_2_1,eye(n3),Zc_3_1)+[zeros(n1+n2) Zb_2_1;zeros(n3+n4,n1+n2) zeros(n3+n4)]-[zeros(n1+n2+n3) Zb_3_1;zeros(n4,n1+n2+n3) zeros(n4)];
		S_3_2_2 = blkdiag(Za_2_2,eye(n3),Zc_3_2)+[zeros(n1+n2) Zb_2_2;zeros(n3+n4,n1+n2) zeros(n3+n4)]-[zeros(n1+n2+n3) Zb_3_2;zeros(n4,n1+n2+n3) zeros(n4)];

		[S_1_1 Zl_1_1';Zl_1_1 S_0_1] >= 0;
		[S_1_2 Zl_1_2';Zl_1_2 S_0_2] >= 0;
		[S_2_1 S_2_1_1;S_2_1_1' S_1_1] >= 0;
		[S_2_2 S_2_1_2;S_2_1_2' S_1_2] >= 0;
		[S_3_1 S_3_2_1;S_3_2_1' S_2_1] >= 0;
		[S_3_2 S_3_2_2;S_3_2_2' S_2_2] >= 0;
		[S_4_1 Zu_3_1;Zu_3_1' S_3_1] >= 0;
		[S_4_2 Zu_3_2;Zu_3_2' S_3_2] >= 0;
	cvx_end

% 	if (strcmp(cvx_status, 'Solved'))
% 		ub = gamm;
% 		if (abs(gamm-lb)<=tol)
% 			cvx_flag = 0;
% 		end
% 	else
% 		lb = gamm;
% 	end
% end

cvx_status

%% 
 
Z0_1 = full(Zl_0_1*(Zu_0_1)^(-1));
Z0_2 = full(Zl_0_2*(Zu_0_2)^(-1));
Z1_1 = full(Zl_1_1*(Zu_1_1)^(-1));
Z1_2 = full(Zl_1_2*(Zu_1_2)^(-1));
Z2_1 = full(Zl_2_1*(Zu_2_1)^(-1));
Z2_2 = full(Zl_2_2*(Zu_2_2)^(-1));
Z3_1 = full(Zl_3_1*(Zu_3_1)^(-1));
Z3_2 = full(Zl_3_2*(Zu_3_2)^(-1));
Z4_1 = full(Zl_4_1*(Zu_4_1)^(-1));
Z4_2 = full(Zl_4_2*(Zu_4_2)^(-1));

nK_1 = n; nK_2 = n; nK_3 = n; nK_4 = n; nK = 4*n;

YC_0_1 = Z0_1;

% i = 1
R11 = Z1_1^-1; S1 = YC_0_1;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_1); R23 = zeros(mn,nK_1);
S11_bar = Z0_1 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_1_1 = (S+S')/2;

% i = 2
R11 = Z2_1^-1; S1 = YC_1_1;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_2); R23 = zeros(mn,nK_2);
S11_bar = Z1_1 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_2_1 = (S+S')/2;

% i = 3
R11 = Z3_1^-1; S1 = YC_2_1;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_3); R23 = zeros(mn,nK_3);
S11_bar = Z2_1 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_3_1 = (S+S')/2;

% i = 4
R11 = Z4_1^-1; S1 = YC_3_1;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_4); R23 = zeros(mn,nK_4);
S11_bar = Z3_1 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_4_1 = (S+S')/2;

YC_1 = YC_4_1; XC_1 = R;

YC_0_2 = Z0_2;

% i = 1
R11 = Z1_2^-1; S1 = YC_0_2;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_1); R23 = zeros(mn,nK_1);
S11_bar = Z0_2 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_1_2 = (S+S')/2;

% i = 2
R11 = Z2_2^-1; S1 = YC_1_2;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_2); R23 = zeros(mn,nK_2);
S11_bar = Z1_2 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_2_2 = (S+S')/2;

% i = 3
R11 = Z3_2^-1; S1 = YC_2_2;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_3); R23 = zeros(mn,nK_3);
S11_bar = Z2_2 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_3_2 = (S+S')/2;

% i = 4
R11 = Z4_2^-1; S1 = YC_3_2;
S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;
S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);
R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_4); R23 = zeros(mn,nK_4);
S11_bar = Z3_2 ^-1;
R13 = chol(R11-S11_bar, 'lower');
R = [R11, R12, R13; R12', R22, R23; R13', R23', R33]; S = R^-1;
YC_4_2 = (S+S')/2;

YC_2 = YC_4_2; XC_2 = R;


%% 
 
R_1 = [A1 zeros(n,nK) Bw1; zeros(nK,n) zeros(nK,nK) zeros(nK,nw); Cz1 zeros(nz,nK) Dzw1];
R_2 = [A2 zeros(n,nK) Bw2; zeros(nK,n) zeros(nK,nK) zeros(nK,nw); Cz2 zeros(nz,nK) Dzw2];
EK_0 = zeros(nK,0);
EK_1 = [eye(nK_1);zeros(nK_2+nK_3+nK_4,nK_1)];
EK_2 = [eye(nK_1+nK_2);zeros(nK_3+nK_4,nK_1+nK_2)];
EK_3 = [eye(nK_1+nK_2+nK_3);zeros(nK_4,nK_1+nK_2+nK_3)];
EK_4 = eye(nK);
EKbar_0 = null(EK_0'); EKbar_1 = null(EK_1'); EKbar_2 = null(EK_2'); EKbar_3 = null(EK_3'); EKbar_4 = null(EK_4'); 

UC_1_1 = [zeros(nK_1+nK_2+nK_3+nK_4,n) EKbar_0' zeros(nK_1+nK_2+nK_3+nK_4,nz); Eubar_0'*Bu1' zeros(nu1+nu2+nu3+nu4,nK) Eubar_0'*Dzu1'];
UC_1_2 = [zeros(nK_1+nK_2+nK_3+nK_4,n) EKbar_0' zeros(nK_1+nK_2+nK_3+nK_4,nz); Eubar_0'*Bu2' zeros(nu1+nu2+nu3+nu4,nK) Eubar_0'*Dzu2'];
UC_2_1 = [zeros(nK_2+nK_3+nK_4,n) EKbar_1' zeros(nK_2+nK_3+nK_4,nz); Eubar_1'*Bu1' zeros(nu2+nu3+nu4,nK) Eubar_1'*Dzu1'];
UC_2_2 = [zeros(nK_2+nK_3+nK_4,n) EKbar_1' zeros(nK_2+nK_3+nK_4,nz); Eubar_1'*Bu2' zeros(nu2+nu3+nu4,nK) Eubar_1'*Dzu2'];
UC_3_1 = [zeros(nK_3+nK_4,n) EKbar_2' zeros(nK_3+nK_4,nz); Eubar_2'*Bu1' zeros(nu3+nu4,nK) Eubar_2'*Dzu1'];
UC_3_2 = [zeros(nK_3+nK_4,n) EKbar_2' zeros(nK_3+nK_4,nz); Eubar_2'*Bu2' zeros(nu3+nu4,nK) Eubar_2'*Dzu2'];
UC_4_1 = [zeros(nK_4,n) EKbar_3' zeros(nK_4,nz); Eubar_3'*Bu1' zeros(nu4,nK) Eubar_3'*Dzu1'];
UC_4_2 = [zeros(nK_4,n) EKbar_3' zeros(nK_4,nz); Eubar_3'*Bu2' zeros(nu4,nK) Eubar_3'*Dzu2'];

VC_1_1 = [zeros(nK_1,n) EK_1' zeros(nK_1,nw); Ey_1'*Cy1 zeros(ny1,nK) Ey_1'*Dyw1];
VC_1_2 = [zeros(nK_1,n) EK_1' zeros(nK_1,nw); Ey_1'*Cy2 zeros(ny1,nK) Ey_1'*Dyw2];
VC_2_1 = [zeros(nK_1+nK_2,n) EK_2' zeros(nK_1+nK_2,nw); Ey_2'*Cy1 zeros(ny1+ny2,nK) Ey_2'*Dyw1];
VC_2_2 = [zeros(nK_1+nK_2,n) EK_2' zeros(nK_1+nK_2,nw); Ey_2'*Cy2 zeros(ny1+ny2,nK) Ey_2'*Dyw2];
VC_3_1 = [zeros(nK_1+nK_2+nK_3,n) EK_3' zeros(nK_1+nK_2+nK_3,nw); Ey_3'*Cy1 zeros(ny1+ny2+ny3,nK) Ey_3'*Dyw1];
VC_3_2 = [zeros(nK_1+nK_2+nK_3,n) EK_3' zeros(nK_1+nK_2+nK_3,nw); Ey_3'*Cy2 zeros(ny1+ny2+ny3,nK) Ey_3'*Dyw2];
VC_4_1 = [zeros(nK_1+nK_2+nK_3+nK_4,n) EK_4' zeros(nK_1+nK_2+nK_3+nK_4,nw); Ey_4'*Cy1 zeros(ny1+ny2+ny3+ny4,nK) Ey_4'*Dyw1];
VC_4_2 = [zeros(nK_1+nK_2+nK_3+nK_4,n) EK_4' zeros(nK_1+nK_2+nK_3+nK_4,nw); Ey_4'*Cy2 zeros(ny1+ny2+ny3+ny4,nK) Ey_4'*Dyw2];
VC_1_1_n = null(VC_1_1); VC_1_2_n = null(VC_1_2); 
VC_2_1_n = null(VC_2_1); VC_2_2_n = null(VC_2_2); 
VC_3_1_n = null(VC_3_1); VC_3_2_n = null(VC_3_2); 
VC_4_1_n = null(VC_4_1); VC_4_2_n = null(VC_4_2); 
VC_0_1_n = eye(n+nK+nw); VC_0_2_n = eye(n+nK+nw); 

eK_1 = [eye(nK_1); zeros(nK_2+nK_3+nK_4,nK_1)];
eK_2 = [zeros(nK_1,nK_2); eye(nK_2); zeros(nK_3+nK_4,nK_2)];
eK_3 = [zeros(nK_1+nK_2,nK_3); eye(nK_3); zeros(nK_4,nK_3)];
eK_4 = [zeros(nK_1+nK_2+nK_3,nK_4); eye(nK_4)];
ey_1 = [eye(ny1); zeros(ny2+ny3+ny4,ny1)];
ey_2 = [zeros(ny1,ny2); eye(ny2); zeros(ny3+ny4,ny2)];
ey_3 = [zeros(ny1+ny2,ny3); eye(ny3); zeros(ny4,ny3)];
ey_4 = [zeros(ny1+ny2+ny3,ny4); eye(ny4)];
VC_1_1_t = [zeros(nK_1,n) eK_1' zeros(nK_1,nw); ey_1'*Cy1 zeros(ny1,nK) ey_1'*Dyw1];
VC_1_2_t = [zeros(nK_1,n) eK_1' zeros(nK_1,nw); ey_1'*Cy2 zeros(ny1,nK) ey_1'*Dyw2];
VC_2_1_t = [zeros(nK_2,n) eK_2' zeros(nK_2,nw); ey_2'*Cy1 zeros(ny2,nK) ey_2'*Dyw1];
VC_2_2_t = [zeros(nK_2,n) eK_2' zeros(nK_2,nw); ey_2'*Cy2 zeros(ny2,nK) ey_2'*Dyw2];
VC_3_1_t = [zeros(nK_3,n) eK_3' zeros(nK_3,nw); ey_3'*Cy1 zeros(ny3,nK) ey_3'*Dyw1];
VC_3_2_t = [zeros(nK_3,n) eK_3' zeros(nK_3,nw); ey_3'*Cy2 zeros(ny3,nK) ey_3'*Dyw2];
VC_4_1_t = [zeros(nK_4,n) eK_4' zeros(nK_4,nw); ey_4'*Cy1 zeros(ny4,nK) ey_4'*Dyw1];
VC_4_2_t = [zeros(nK_4,n) eK_4' zeros(nK_4,nw); ey_4'*Cy2 zeros(ny4,nK) ey_4'*Dyw2];

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_4_11_t(nK_4+nu4,nK_4+ny4)
	B1 = [blkdiag(YC_1,eye(nz)) (R_1+UC_4_1'*QK_4_11_t*VC_4_1_t)*VC_3_1_n; VC_3_1_n'*(VC_4_1_t'*QK_4_11_t'*UC_4_1+ R_1') VC_3_1_n'*blkdiag(XC_1, eye(nw))*VC_3_1_n];
	(B1+B1')/2 >= 0;
	%minimize(norm(QK_4_11_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_4_12_t(nK_4+nu4,nK_4+ny4)
	B2 = [blkdiag(YC_2,eye(nz)) (R_2+UC_4_2'*QK_4_12_t*VC_4_2_t)*VC_3_2_n; VC_3_2_n'*(VC_4_2_t'*QK_4_12_t'*UC_4_2+ R_2') VC_3_2_n'*blkdiag(XC_1, eye(nw))*VC_3_2_n];
	(B2+B2')/2 >= 0;
	%minimize(norm(QK_4_12_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_4_21_t(nK_4+nu4,nK_4+ny4)
	B3 = [blkdiag(YC_1,eye(nz)) (R_1+UC_4_1'*QK_4_21_t*VC_4_1_t)*VC_3_1_n; VC_3_1_n'*(VC_4_1_t'*QK_4_21_t'*UC_4_1+ R_1') VC_3_1_n'*blkdiag(XC_2, eye(nw))*VC_3_1_n];
	(B3+B3')/2 >= 0;
	%minimize(norm(QK_4_21_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_4_22_t(nK_4+nu4,nK_4+ny4)
	B4 = [blkdiag(YC_2,eye(nz)) (R_2+UC_4_2'*QK_4_22_t*VC_4_2_t)*VC_3_2_n; VC_3_2_n'*(VC_4_2_t'*QK_4_22_t'*UC_4_2+ R_2') VC_3_2_n'*blkdiag(XC_2, eye(nw))*VC_3_2_n];
	(B4+B4')/2 >= 0;
	%minimize(norm(QK_4_22_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_3_11_t(nK_3+nK_4+nu3+nu4,nK_3+ny3)
	B5 = [blkdiag(YC_1,eye(nz)) (R_1+UC_3_1'*QK_3_11_t*VC_3_1_t+UC_4_1'*QK_4_11_t*VC_4_1_t)*VC_2_1_n; VC_2_1_n'*(VC_3_1_t'*QK_3_11_t'*UC_3_1+VC_4_1_t'*QK_4_11_t'*UC_4_1+ R_1') VC_2_1_n'*blkdiag(XC_1, eye(nw))*VC_2_1_n];
	(B5+B5')/2 >= 0;
	%minimize(norm(QK_3_11_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_3_12_t(nK_3+nK_4+nu3+nu4,nK_3+ny3)
	B6 = [blkdiag(YC_2,eye(nz)) (R_2+UC_3_2'*QK_3_12_t*VC_3_2_t+UC_4_2'*QK_4_12_t*VC_4_2_t)*VC_2_2_n; VC_2_2_n'*(VC_3_2_t'*QK_3_12_t'*UC_3_2+VC_4_2_t'*QK_4_12_t'*UC_4_2+ R_2') VC_2_2_n'*blkdiag(XC_1, eye(nw))*VC_2_2_n];
	(B6+B6')/2 >= 0;
	%minimize(norm(QK_3_12_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_3_21_t(nK_3+nK_4+nu3+nu4,nK_3+ny3)
	B7 = [blkdiag(YC_1,eye(nz)) (R_1+UC_3_1'*QK_3_21_t*VC_3_1_t+UC_4_1'*QK_4_21_t*VC_4_1_t)*VC_2_1_n; VC_2_1_n'*(VC_3_1_t'*QK_3_21_t'*UC_3_1+VC_4_1_t'*QK_4_21_t'*UC_4_1+ R_1') VC_2_1_n'*blkdiag(XC_2, eye(nw))*VC_2_1_n];
	(B7+B7')/2 >= 0;
	%minimize(norm(QK_3_21_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_3_22_t(nK_3+nK_4+nu3+nu4,nK_3+ny3)
	B8 = [blkdiag(YC_2,eye(nz)) (R_2+UC_3_2'*QK_3_22_t*VC_3_2_t+UC_4_2'*QK_4_22_t*VC_4_2_t)*VC_2_2_n; VC_2_2_n'*(VC_3_2_t'*QK_3_22_t'*UC_3_2+VC_4_2_t'*QK_4_22_t'*UC_4_2+ R_2') VC_2_2_n'*blkdiag(XC_2, eye(nw))*VC_2_2_n];
	(B8+B8')/2 >= 0;
	%minimize(norm(QK_3_22_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_2_11_t(nK_2+nK_3+nK_4+nu2+nu3+nu4,nK_2+ny2)
	B9 = [blkdiag(YC_1,eye(nz)) (R_1+UC_2_1'*QK_2_11_t*VC_2_1_t+UC_3_1'*QK_3_11_t*VC_3_1_t+UC_4_1'*QK_4_11_t*VC_4_1_t)*VC_1_1_n; VC_1_1_n'*(VC_2_1_t'*QK_2_11_t'*UC_2_1+VC_3_1_t'*QK_3_11_t'*UC_3_1+VC_4_1_t'*QK_4_11_t'*UC_4_1+ R_1') VC_1_1_n'*blkdiag(XC_1, eye(nw))*VC_1_1_n];
	(B9+B9')/2 >= 0;
	%minimize(norm(QK_2_11_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_2_12_t(nK_2+nK_3+nK_4+nu2+nu3+nu4,nK_2+ny2)
	B10 = [blkdiag(YC_2,eye(nz)) (R_2+UC_2_2'*QK_2_12_t*VC_2_2_t+UC_3_2'*QK_3_12_t*VC_3_2_t+UC_4_2'*QK_4_12_t*VC_4_2_t)*VC_1_2_n; VC_1_2_n'*(VC_2_2_t'*QK_2_12_t'*UC_2_2+VC_3_2_t'*QK_3_12_t'*UC_3_2+VC_4_2_t'*QK_4_12_t'*UC_4_2+ R_2') VC_1_2_n'*blkdiag(XC_1, eye(nw))*VC_1_2_n];
	(B10+B10')/2 >= 0;
	%minimize(norm(QK_2_12_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_2_21_t(nK_2+nK_3+nK_4+nu2+nu3+nu4,nK_2+ny2)
	B11 = [blkdiag(YC_1,eye(nz)) (R_1+UC_2_1'*QK_2_21_t*VC_2_1_t+UC_3_1'*QK_3_21_t*VC_3_1_t+UC_4_1'*QK_4_21_t*VC_4_1_t)*VC_1_1_n; VC_1_1_n'*(VC_2_1_t'*QK_2_21_t'*UC_2_1+VC_3_1_t'*QK_3_21_t'*UC_3_1+VC_4_1_t'*QK_4_21_t'*UC_4_1+ R_1') VC_1_1_n'*blkdiag(XC_2, eye(nw))*VC_1_1_n];
	(B11+B11')/2 >= 0;
	%minimize(norm(QK_2_21_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_2_22_t(nK_2+nK_3+nK_4+nu2+nu3+nu4,nK_2+ny2)
	B12 = [blkdiag(YC_2,eye(nz)) (R_2+UC_2_2'*QK_2_22_t*VC_2_2_t+UC_3_2'*QK_3_22_t*VC_3_2_t+UC_4_2'*QK_4_22_t*VC_4_2_t)*VC_1_2_n; VC_1_2_n'*(VC_2_2_t'*QK_2_22_t'*UC_2_2+VC_3_2_t'*QK_3_22_t'*UC_3_2+VC_4_2_t'*QK_4_22_t'*UC_4_2+ R_2') VC_1_2_n'*blkdiag(XC_2, eye(nw))*VC_1_2_n];
	(B12+B12')/2 >= 0;
	%minimize(norm(QK_2_22_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_1_11_t(nK_1+nK_2+nK_3+nK_4+nu1+nu2+nu3+nu4,nK_1+ny1)
	B13 = [blkdiag(YC_1,eye(nz)) (R_1+UC_1_1'*QK_1_11_t*VC_1_1_t+UC_2_1'*QK_2_11_t*VC_2_1_t+UC_3_1'*QK_3_11_t*VC_3_1_t+UC_4_1'*QK_4_11_t*VC_4_1_t)*VC_0_1_n; VC_0_1_n'*(VC_1_1_t'*QK_1_11_t'*UC_1_1+VC_2_1_t'*QK_2_11_t'*UC_2_1+VC_3_1_t'*QK_3_11_t'*UC_3_1+VC_4_1_t'*QK_4_11_t'*UC_4_1+ R_1') VC_0_1_n'*blkdiag(XC_1, eye(nw))*VC_0_1_n];
	(B13+B13')/2 >= 0;
	%minimize(norm(QK_1_11_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_1_12_t(nK_1+nK_2+nK_3+nK_4+nu1+nu2+nu3+nu4,nK_1+ny1)
	B14 = [blkdiag(YC_2,eye(nz)) (R_2+UC_1_2'*QK_1_12_t*VC_1_2_t+UC_2_2'*QK_2_12_t*VC_2_2_t+UC_3_2'*QK_3_12_t*VC_3_2_t+UC_4_2'*QK_4_12_t*VC_4_2_t)*VC_0_2_n; VC_0_2_n'*(VC_1_2_t'*QK_1_12_t'*UC_1_2+VC_2_2_t'*QK_2_12_t'*UC_2_2+VC_3_2_t'*QK_3_12_t'*UC_3_2+VC_4_2_t'*QK_4_12_t'*UC_4_2+ R_2') VC_0_2_n'*blkdiag(XC_1, eye(nw))*VC_0_2_n];
	(B14+B14')/2 >= 0;
	%minimize(norm(QK_1_12_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_1_21_t(nK_1+nK_2+nK_3+nK_4+nu1+nu2+nu3+nu4,nK_1+ny1)
	B15 = [blkdiag(YC_1,eye(nz)) (R_1+UC_1_1'*QK_1_21_t*VC_1_1_t+UC_2_1'*QK_2_21_t*VC_2_1_t+UC_3_1'*QK_3_21_t*VC_3_1_t+UC_4_1'*QK_4_21_t*VC_4_1_t)*VC_0_1_n; VC_0_1_n'*(VC_1_1_t'*QK_1_21_t'*UC_1_1+VC_2_1_t'*QK_2_21_t'*UC_2_1+VC_3_1_t'*QK_3_21_t'*UC_3_1+VC_4_1_t'*QK_4_21_t'*UC_4_1+ R_1') VC_0_1_n'*blkdiag(XC_2, eye(nw))*VC_0_1_n];
	(B15+B15')/2 >= 0;
	%minimize(norm(QK_1_21_t));
cvx_end
cvx_status

cvx_begin sdp quiet
	cvx_solver Mosek
	variable QK_1_22_t(nK_1+nK_2+nK_3+nK_4+nu1+nu2+nu3+nu4,nK_1+ny1)
	B16 = [blkdiag(YC_2,eye(nz)) (R_2+UC_1_2'*QK_1_22_t*VC_1_2_t+UC_2_2'*QK_2_22_t*VC_2_2_t+UC_3_2'*QK_3_22_t*VC_3_2_t+UC_4_2'*QK_4_22_t*VC_4_2_t)*VC_0_2_n; VC_0_2_n'*(VC_1_2_t'*QK_1_22_t'*UC_1_2+VC_2_2_t'*QK_2_22_t'*UC_2_2+VC_3_2_t'*QK_3_22_t'*UC_3_2+VC_4_2_t'*QK_4_22_t'*UC_4_2+ R_2') VC_0_2_n'*blkdiag(XC_2, eye(nw))*VC_0_2_n];
	(B16+B16')/2 >= 0;
	%minimize(norm(QK_1_22_t));
cvx_end
cvx_status


%% 
 
QK_11 = blkdiag(EKbar_0,Eubar_0)*QK_1_11_t*blkdiag(eK_1,ey_1)'+blkdiag(EKbar_1,Eubar_1)*QK_2_11_t*blkdiag(eK_2,ey_2)'+blkdiag(EKbar_2,Eubar_2)*QK_3_11_t*blkdiag(eK_3,ey_3)'+blkdiag(EKbar_3,Eubar_3)*QK_4_11_t*blkdiag(eK_4,ey_4)';
QK_12 = blkdiag(EKbar_0,Eubar_0)*QK_1_12_t*blkdiag(eK_1,ey_1)'+blkdiag(EKbar_1,Eubar_1)*QK_2_12_t*blkdiag(eK_2,ey_2)'+blkdiag(EKbar_2,Eubar_2)*QK_3_12_t*blkdiag(eK_3,ey_3)'+blkdiag(EKbar_3,Eubar_3)*QK_4_12_t*blkdiag(eK_4,ey_4)';
QK_21 = blkdiag(EKbar_0,Eubar_0)*QK_1_21_t*blkdiag(eK_1,ey_1)'+blkdiag(EKbar_1,Eubar_1)*QK_2_21_t*blkdiag(eK_2,ey_2)'+blkdiag(EKbar_2,Eubar_2)*QK_3_21_t*blkdiag(eK_3,ey_3)'+blkdiag(EKbar_3,Eubar_3)*QK_4_21_t*blkdiag(eK_4,ey_4)';
QK_22 = blkdiag(EKbar_0,Eubar_0)*QK_1_22_t*blkdiag(eK_1,ey_1)'+blkdiag(EKbar_1,Eubar_1)*QK_2_22_t*blkdiag(eK_2,ey_2)'+blkdiag(EKbar_2,Eubar_2)*QK_3_22_t*blkdiag(eK_3,ey_3)'+blkdiag(EKbar_3,Eubar_3)*QK_4_22_t*blkdiag(eK_4,ey_4)';
AK(:,:,1) = QK_11(1:nK,1:nK); BK(:,:,1) = QK_11(1:nK,nK+1:end)/gamm; CK(:,:,1) = QK_11(nK+1:end,1:nK); DK(:,:,1) = QK_11(nK+1:end,nK+1:end)/gamm;
AK(:,:,2) = QK_12(1:nK,1:nK); BK(:,:,2) = QK_12(1:nK,nK+1:end)/gamm; CK(:,:,2) = QK_12(nK+1:end,1:nK); DK(:,:,2) = QK_12(nK+1:end,nK+1:end)/gamm;
AK(:,:,3) = QK_21(1:nK,1:nK); BK(:,:,3) = QK_21(1:nK,nK+1:end)/gamm; CK(:,:,3) = QK_21(nK+1:end,1:nK); DK(:,:,3) = QK_21(nK+1:end,nK+1:end)/gamm;
AK(:,:,4) = QK_22(1:nK,1:nK); BK(:,:,4) = QK_22(1:nK,nK+1:end)/gamm; CK(:,:,4) = QK_22(nK+1:end,1:nK); DK(:,:,4) = QK_22(nK+1:end,nK+1:end)/gamm;

AK1 = AK(1:nK_1,1:nK_1,:); BK1 = BK(1:nK_1,1:ny1,:); CK1 = CK(+1:nu1,1:nK_1,:); DK1 = DK(+1:nu1,1:ny1,:);
AK2 = AK(1:nK_1+nK_2,1:nK_1+nK_2,:); BK2 = BK(1:nK_1+nK_2,1:ny1+ny2,:); CK2 = CK(nu1+1:nu1+nu2,1:nK_1+nK_2,:); DK2 = DK(nu1+1:nu1+nu2,1:ny1+ny2,:);
AK3 = AK(1:nK_1+nK_2+nK_3,1:nK_1+nK_2+nK_3,:); BK3 = BK(1:nK_1+nK_2+nK_3,1:ny1+ny2+ny3,:); CK3 = CK(nu1+nu2+1:nu1+nu2+nu3,1:nK_1+nK_2+nK_3,:); DK3 = DK(nu1+nu2+1:nu1+nu2+nu3,1:ny1+ny2+ny3,:);
AK4 = AK(1:nK_1+nK_2+nK_3+nK_4,1:nK_1+nK_2+nK_3+nK_4,:); BK4 = BK(1:nK_1+nK_2+nK_3+nK_4,1:ny1+ny2+ny3+ny4,:); CK4 = CK(nu1+nu2+nu3+1:nu1+nu2+nu3+nu4,1:nK_1+nK_2+nK_3+nK_4,:); DK4 = DK(nu1+nu2+nu3+1:nu1+nu2+nu3+nu4,1:ny1+ny2+ny3+ny4,:);

% AK1(abs(AK1)<1e-6) = 0; BK1(abs(BK1)<1e-6) = 0; CK1(abs(CK1)<1e-6) = 0; DK1(abs(DK1)<1e-6) = 0;
% AK2(abs(AK2)<1e-6) = 0; BK2(abs(BK2)<1e-6) = 0; CK2(abs(CK2)<1e-6) = 0; DK2(abs(DK2)<1e-6) = 0;
% AK3(abs(AK3)<1e-6) = 0; BK3(abs(BK3)<1e-6) = 0; CK3(abs(CK3)<1e-6) = 0; DK3(abs(DK3)<1e-6) = 0;
% AK4(abs(AK4)<1e-6) = 0; BK4(abs(BK4)<1e-6) = 0; CK4(abs(CK4)<1e-6) = 0; DK4(abs(DK4)<1e-6) = 0;

all(abs(eig(AK4(:,:,1)))<1)
all(abs(eig(AK4(:,:,2)))<1)
all(abs(eig(AK4(:,:,3)))<1)
all(abs(eig(AK4(:,:,4)))<1)

%% Model parameters of the crazyflie 2.0
init_simulink
Psi = [11; 12; 21; 22];