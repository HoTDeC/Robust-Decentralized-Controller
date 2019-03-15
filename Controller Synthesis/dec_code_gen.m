

possible_switch = num2str([11 12 21 22]'); % add possible switching sequences
start = ['1', '2']'; % starting modes (normally all modes)
h = '12'; %all modes


L = 1;  % Memory length
H = 1;  % horizon length
Psi = switching_seq_gen(start,L+H+1,possible_switch);
Phi = switching_seq_gen(start,L+H,possible_switch);

Phi_lbar = Psi(:,1:end-1);
Phi_ubar = Psi(:,2:end);
Phi_star = Psi(:,L+1);

p = size(Psi,1);
q = size(Phi,1);


% 
M = 4; %number of robots
J = '1234'; % string from 1 to M
J_bar = '01234'; % string from 0 to M

%%
% Add file name here
fid=fopen(['cf_dec_' num2str(M) 'p_L' num2str(L) '_H' num2str(H) '_' num2str(length(h)) 'modes.m'],'w');

%%

for i = J
   if i == J(1)
       nb = ['n' i]; nby = ['ny' i]; nbu = ['nu' i]; nbK = ['nK_' i]; nbz = ['nz' i]; nbw = ['nw' i];
   else
       nb = [nb '+n' i]; nby = [nby '+ny' i]; nbu = [nbu '+nu' i]; nbK = [nbK '+nK_' i]; nbz = [nbz '+nz' i]; nbw = [nbw '+nw' i];
   end
end


fprintf(fid,'%%%% Add own model here \r\n \r\n');

%%
fprintf(fid,'\r\n%%%% \r\n \r\n');

fprintf(fid, 'Ey_0 = zeros(ny,0);\r\n');
for i = J(1:end-1)
    fprintf(fid, ['Ey_' i ' = [eye(' nby(1:4*str2double(i)-1) ');zeros(' nby(4*str2double(i)+1:end) ',' nby(1:4*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['Ey_' J(end) ' = eye(ny);\r\n\r\n']);

fprintf(fid, 'Eu_0 = zeros(nu,0);\r\n');
for i = J(1:end-1)
    fprintf(fid, ['Eu_' i ' = [eye(' nbu(1:4*str2double(i)-1) ');zeros(' nbu(4*str2double(i)+1:end) ',' nbu(1:4*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['Eu_' J(end) ' = eye(nu);\r\n\r\n']);

fprintf(fid, 'E_0 = zeros(n,0);\r\n');
for i = J(1:end-1)
    fprintf(fid, ['E_' i ' = [eye(' nb(1:3*str2double(i)-1) ');zeros(' nb(3*str2double(i)+1:end) ',' nb(1:3*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['E_' J(end) ' = eye(n);\r\n\r\n']);

for i = J_bar
    fprintf(fid, ['Eybar_' i ' = null(Ey_' i '''); ']);
end
fprintf(fid,'\r\n');
for i = J_bar
    fprintf(fid, ['Eubar_' i ' = null(Eu_' i '''); ']);
end
fprintf(fid,'\r\n');
for i = J_bar
    fprintf(fid, ['Ebar_' i ' = null(E_' i '''); ']);
end
fprintf(fid,'\r\n\r\n');

for i = J_bar %number of players
    for j = h %number of modes
        fprintf(fid, ['Ny_' i '_' j ' = null([Ey_' i '''*Cy' j ' Ey_' i '''*Dyw' j ']);\r\n']);
    end
end
for i = J_bar %number of players
    for j = h %number of modes
        fprintf(fid, ['A11_' i '_' j ' = E_' i '''*A' j '*E_' i '; A21_' i '_' j ' = Ebar_' i '''*A' j '*E_' i '; A22_' i '_' j ' = Ebar_' i '''*A' j '*Ebar_' i ';\r\n']);
    end
end

for j = h
    fprintf(fid,['Czf' j ' = Cz' j '; Cyf' j ' = Cy' j '; Dzwf' j ' = Dzw' j '; Dzuf' j ' = Dzu' j '; Dywf' j ' = Dyw' j ';\r\n']);
end
%%
fprintf(fid,'\r\n%%%% \r\n \r\n');

fprintf(fid,'lb = 0;\r\nub = 10;\r\ntol = 1e-3;\r\ncvx_flag = 1;\r\n\r\nwhile(cvx_flag)\r\ngamm = (ub+lb)/2\r\n');
fprintf(fid,'\tcvx_begin sdp quiet \r\n\t\tcvx_solver Mosek \r\n');

for i = J_bar(2:end)
    for j = 1:q
        fprintf(fid, ['\t\tvariable Za_' i '_' Phi(j,:) '(' nb(1:3*str2double(i)-1) ',' nb(1:3*str2double(i)-1) ') semidefinite \r\n']);
    end
end

for i = J_bar(2:end-1)
    for j = 1:q
        fprintf(fid, ['\t\tvariable Zb_' i '_' Phi(j,:) '(' nb(1:3*str2double(i)-1) ',' nb(3*str2double(i)+1:end) ') \r\n']);
    end
end

for i = J_bar(1:end-1)
    for j = 1:q
        fprintf(fid, ['\t\tvariable Zc_' num2str(i) '_' Phi(j,:) '(' nb(3*str2double(i)+1:end) ',' nb(3*str2double(i)+1:end) ') semidefinite \r\n']);
    end
end

fprintf(fid, '\r\n');
for j = h
    fprintf(fid, ['\t\tCz' j ' = Czf' j '/gamm; Cy' j ' = Cyf' j '/gamm; Dzw' j ' = Dzwf' j '/gamm; Dzu' j ' = Dzuf' j '/gamm; Dyw' j ' = Dywf' j '/gamm; \r\n']);
end
fprintf(fid, '\r\n');

for i = J_bar %number of players
    for j = h %number of modes
        fprintf(fid, ['\t\tNu_' i '_' j ' = null([Eubar_' i '''*Bu' j ''' Eubar_' i '''*Dzu' j ''']);\r\n']);
    end
end

fprintf(fid, '\r\n');

% Fix extreme cases (i=0, i=end)
% for j = h
%     fprintf(fid, ['\t\tG_0_' j ' = [Nu_0_' j ' zeros(n+nz,n+nw); zeros(n+nw,n+nz-(' nbu(1:4*str2double(J_bar(end))-1) ')) Ny_0_' j '];\r\n']);
% end
% for i = J_bar(2:end-1)
%     u = num2str(M-str2double(i));
%     for j = h
%         fprintf(fid, ['\t\tG_' i '_' j ' = [Nu_' i '_' j ' zeros(n+nz,n+nw-(' nby(1:4*str2double(i)-1) ')); zeros(n+nw,n+nz-(' nbu(1:4*(str2double(J_bar(end))-str2double(i))-1) ')) Ny_' i '_' j '];\r\n']);
%     end
% end
% for j = h
%     i = J_bar(end);
%     fprintf(fid, ['\t\tG_' i '_' j ' = [Nu_' i '_' j ' zeros(n+nz,n+nw-(' nby(1:4*str2double(i)-1) ')); zeros(n+nw,n+nz) Ny_' i '_' j '];\r\n']);
% end
%original:
% for i = J_bar
%     u = num2str(M-str2double(i));
%     for j = h
%         fprintf(fid, ['\t\tG_' i '_' j ' = [Nu_' i '_' j ' zeros(n+nz,n+nw-(' nby(1:4*str2double(i)-1) ')); zeros(n+nw,n+nz-(' nbu(1:4*(str2double(J_bar(end))-str2double(i))-1) ')) Ny_' i '_' j '];\r\n']);
%     end
% end
for i = J_bar
    u = num2str(M-str2double(i));
    for j = h
        fprintf(fid, ['\t\tG_' i '_' j ' = blkdiag(Nu_' i '_' j ', Ny_' i '_' j ');\r\n']);
    end
end


fprintf(fid, '\r\n');

for j = 1:q
    fprintf(fid, ['\t\tZl_' J_bar(1) '_' Phi(j,:) ' = Zc_' J_bar(1) '_' Phi(j,:) '; Zu_' J_bar(1) '_' Phi(j,:) ' = eye(n); \r\n']);
end
for i = J_bar(2:end-1)
    for j = 1:q
        fprintf(fid, ['\t\tZl_' i '_' Phi(j,:) ' = [eye(' nb(1:3*str2double(i)-1) ') zeros(' nb(1:3*str2double(i)-1) ',' nb(3*str2double(i)+1:end) '); -Zb_' i '_' Phi(j,:) ''' '...
            'Zc_' i '_' Phi(j,:) ']; Zu_' i '_' Phi(j,:) ' = [Za_' i '_' Phi(j,:) ' Zb_' i '_' Phi(j,:) '; zeros(' nb(3*str2double(i)+1:end) ',' nb(1:3*str2double(i)-1) ') eye(' nb(3*str2double(i)+1:end) ')];\r\n']);
    end
end
for j = 1:q
    fprintf(fid, ['\t\tZl_' J_bar(end) '_' Phi(j,:) ' = eye(n); Zu_' J_bar(end) '_' Phi(j,:) ' = Za_' J_bar(end) '_' Phi(j,:) '; \r\n']);
end

fprintf(fid, '\r\n');

for j = 1:q
    fprintf(fid, ['\t\tS_' J_bar(1) '_' Phi(j,:) ' = Zc_' J_bar(1) '_' Phi(j,:) ';\r\n']);
end
% for i = J_bar(2:end-1)
%     for j = 1:q
%         fprintf(fid, ['\t\tS_' i '_' Phi(j,:) ' = [Za_' i '_' Phi(j,:) ' zeros(' nb(1:3*str2double(i)-1) ',' nb(3*str2double(i)+1:end) '); '...
%             'zeros(' nb(3*str2double(i)+1:end) ',' nb(1:3*str2double(i)-1) ') Zc_' i '_' Phi(j,:) '];\r\n']);
%     end
% end
for i = J_bar(2:end-1)
    for j = 1:q
        fprintf(fid, ['\t\tS_' i '_' Phi(j,:) ' = blkdiag(Za_' i '_' Phi(j,:) ', Zc_' i '_' Phi(j,:) ');\r\n']);
    end
end
for j = 1:q
    fprintf(fid, ['\t\tS_' J_bar(end) '_' Phi(j,:) ' = Za_' J_bar(end) '_' Phi(j,:) ';\r\n']);
end

fprintf(fid, '\r\n');

for j = 1:p
    fprintf(fid, ['\t\tZAZ_' J_bar(1) '_' Psi(j,:) ' = A22_' J_bar(1) '_' Phi_star(j,:) '*Zc_' J_bar(1) '_' Phi_lbar(j,:) '; ZAZ_' J_bar(1) '_' Psi(j,:) '_t = Zc_' J_bar(1) '_' Phi_lbar(j,:) '*A22_' J_bar(1) '_' Phi_star(j,:) ''';\r\n']);
end
for i = J_bar(2:end-1)
    for j = 1:p
        fprintf(fid, ['\t\tZAZ_' i '_' Psi(j,:) ' = [Za_' i '_' Phi_ubar(j,:) '*A11_' i '_' Phi_star(j,:) ' zeros(' nb(1:3*str2double(i)-1) ',' ...
            nb(3*str2double(i)+1:end) '); Zb_' i '_' Phi_ubar(j,:) '''*A11_' i '_' Phi_star(j,:) '+A21_' i '_' Phi_star(j,:) '-A22_' i '_' Phi_star(j,:) '*Zb_' i '_' Phi_lbar(j,:) ''' A22_' i '_' Phi_star(j,:) '*Zc_' i '_' Phi_lbar(j,:) ']; ' ...
            'ZAZ_' i '_' Psi(j,:) '_t = [A11_' i '_' Phi_star(j,:) '''*Za_' i '_' Phi_ubar(j,:) ' A11_' i '_' Phi_star(j,:) '''*Zb_' i '_' Phi_ubar(j,:) '+A21_' i '_' Phi_star(j,:) ...
            '''-Zb_' i '_' Phi_lbar(j,:) '*A22_' i '_' Phi_star(j,:) '''; zeros(' nb(3*str2double(i)+1:end) ',' nb(1:3*str2double(i)-1) ') Zc_' i '_' Phi_lbar(j,:) '*A22_' i '_' Phi_star(j,:) '''];\r\n']);
    end
end
for j = 1:p
    fprintf(fid, ['\t\tZAZ_' J_bar(end) '_' Psi(j,:) ' = Za_' J_bar(end) '_' Phi_ubar(j,:) '*A11_' J_bar(end) '_' Phi_star(j,:) '; ZAZ_' J_bar(end) '_' Psi(j,:) '_t = A11_' J_bar(end) '_' Phi_star(j,:) '''*Za_' J_bar(end) '_' Phi_ubar(j,:) ';\r\n']);
end

fprintf(fid, '\r\n');

for i = J_bar %number of players
    for j = 1:p %number of possible sequences
        fprintf(fid, ['\t\tM' i '_' Psi(j,:) ' = [S_' i '_' Phi_ubar(j,:) ' zeros(n,nz) ZAZ_' i '_' Psi(j,:) ' Zu_' i '_' Phi_ubar(j,:) '''*Bw' Phi_star(j,:) ';'...
            'zeros(nz,n) eye(nz) Cz' Phi_star(j,:) '*Zl_' i '_' Phi_lbar(j,:) ' Dzw' Phi_star(j,:) '; ZAZ_' i '_' Psi(j,:) '_t Zl_' i '_' Phi_lbar(j,:) '''*Cz' Phi_star(j,:) ''' S_' i '_' Phi_lbar(j,:) ' zeros(n,nw);'...
            'Bw' Phi_star(j,:) '''*Zu_' i '_' Phi_ubar(j,:) ' Dzw' Phi_star(j,:) ''' zeros(nw,n) eye(nw)];\r\n']);
    end
end

fprintf(fid, '\r\n');

for i = J_bar %number of players
    for j = 1:p %number of possible sequences
        k = num2str(j);
        fprintf(fid, ['\t\tG_' i '_' Phi_star(j,:) '''*M' i '_' Psi(j,:) '*G_' i '_' Phi_star(j,:) ' >= 0;\r\n']);
    end
end

fprintf(fid, '\r\n');

for i = J(2:end-1)
    for j = 1:q
        i_1 = num2str(str2double(i)-1);
        fprintf(fid, ['\t\tS_' i '_' i_1 '_' Phi(j,:) ' = blkdiag(Za_' i_1 '_' Phi(j,:) ',eye(n' i '),Zc_' i '_' Phi(j,:) ')+[zeros('...
            nb(1:3*str2double(i_1)-1) ') Zb_' i_1 '_' Phi(j,:) ';zeros(' nb(3*str2double(i_1)+1:end) ',' nb(1:3*str2double(i_1)-1) ') zeros(' nb(3*str2double(i_1)+1:end) ...
            ')]-[zeros(' nb(1:3*str2double(i)-1) ') Zb_' i '_' Phi(j,:) ';zeros(' nb(3*str2double(i)+1:end) ',' nb(1:3*str2double(i)-1) ...
            ') zeros(' nb(3*str2double(i)+1:end) ')];\r\n']);
    end
end

fprintf(fid, '\r\n');

for j = 1:q
    fprintf(fid, ['\t\t[S_' J_bar(2) '_' Phi(j,:) ' Zl_' J_bar(2) '_' Phi(j,:) ''';Zl_' J_bar(2) '_' Phi(j,:) ' S_' J_bar(1) '_' Phi(j,:) '] >= 0;\r\n']);
end
for i = J(2:end-1) %number of players
    for j = 1:q
        i_1 = num2str(str2double(i)-1);
        fprintf(fid, ['\t\t[S_' i '_' Phi(j,:) ' S_' i '_' i_1 '_' Phi(j,:) ';S_' i '_' i_1 '_' Phi(j,:) ''' S_' i_1 '_' Phi(j,:) '] >= 0;\r\n']);
    end
end
for j = 1:q
    fprintf(fid, ['\t\t[S_' J_bar(end) '_' Phi(j,:) ' Zu_' J_bar(end-1) '_' Phi(j,:) ';Zu_' J_bar(end-1) '_' Phi(j,:) ''' S_' J_bar(end-1) '_' Phi(j,:) '] >= 0;\r\n']);
end

fprintf(fid,'\tcvx_end\r\n\r\n');
fprintf(fid,['\tif (strcmp(cvx_status, ''Solved''))\r\n\t\tub = gamm;\r\n\t\tif (abs(gamm-lb)<=tol)' ...
    '\r\n\t\t\tcvx_flag = 0;\r\n\t\tend\r\n\telse\r\n\t\tlb = gamm;\r\n\tend\r\n']);
fprintf(fid,'end\r\n\r\ncvx_status\r\n');
%%
fprintf(fid,'\r\n%%%% \r\n \r\n');

for i = J_bar
    for j = 1:q
        fprintf(fid, ['Z' i '_' Phi(j,:) ' = full(Zl_' i '_' Phi(j,:) '*(Zu_' i '_' Phi(j,:) ')^(-1));\r\n']);
    end
end

fprintf(fid,['\r\nnK_1 = n; nK_2 = n; nK_3 = n; nK_4 = n; nK = ' num2str(M) '*n;\r\n\r\n']);

for j = 1:q
    fprintf(fid, ['YC_0_' Phi(j,:) ' = Z0_' Phi(j,:) ';\r\n']);
    for i = J
        i_1 = num2str(str2double(i)-1);
        fprintf(fid, ['\r\n%% i = ' i '\r\n']);
        fprintf(fid, ['R11 = Z' i '_' Phi(j,:) '^-1; S1 = YC_' i_1 '_' Phi(j,:) ';\r\n']);
        fprintf(fid, 'S11 = S1(1:n,1:n); S12 = S1(1:n,n+1:end); S22 = S1(n+1:end,n+1:end); mn = size(S1,1)-n;\r\n');
        fprintf(fid, 'S1_inv = S1^-1; S11_bar = S1_inv(1:n,1:n); S12_bar = S1_inv(1:n,n+1:end); S22_bar = S1_inv(n+1:end,n+1:end);\r\n');
        fprintf(fid, ['R12 = S12_bar; R22 = S22_bar; R33 = eye(nK_' i '); R23 = zeros(mn,nK_' i ');\r\n']);
        fprintf(fid, ['S11_bar = Z' i_1 '_' Phi(j,:) ' ^-1;\r\n']);
        fprintf(fid, 'R13 = chol(R11-S11_bar, ''lower'');\r\n');
        fprintf(fid, 'R = [R11, R12, R13; R12'', R22, R23; R13'', R23'', R33]; S = R^-1;\r\n');
        fprintf(fid, ['YC_' i '_' Phi(j,:) ' = (S+S'')/2;\r\n']);
    end
    fprintf(fid, ['\r\nYC_' Phi(j,:) ' = YC_' J(end) '_' Phi(j,:) '; XC_' Phi(j,:) ' = R;\r\n\r\n']);
end

%%
fprintf(fid,'\r\n%%%% \r\n \r\n');

for j = h
    fprintf(fid, ['R_' j ' = [A' j ' zeros(n,nK) Bw' j '; zeros(nK,n) zeros(nK,nK) zeros(nK,nw); Cz' j ' zeros(nz,nK) Dzw' j '];\r\n']);
end

fprintf(fid, 'EK_0 = zeros(nK,0);\r\n');
for i = J(1:end-1)
    fprintf(fid, ['EK_' i ' = [eye(' nbK(1:5*str2double(i)-1) ');zeros(' nbK(5*str2double(i)+1:end) ',' nbK(1:5*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['EK_' J(end) ' = eye(nK);\r\n']);
for i = J_bar
    fprintf(fid, ['EKbar_' i ' = null(EK_' i '''); ']);
end
fprintf(fid,'\r\n\r\n');

for i = J %number of players
    i_1 = num2str(str2double(i)-1);
    for j = h %number of modes
        fprintf(fid, ['UC_' i '_' j ' = [zeros(' nbK(5*str2double(i_1)+1:end) ',n) EKbar_' i_1 ''' zeros(' nbK(5*str2double(i_1)+1:end) ',nz); '...
            'Eubar_' i_1 '''*Bu' j ''' zeros(' nbu(4*str2double(i_1)+1:end) ',nK) Eubar_' i_1 '''*Dzu' j '''];\r\n']);
    end
end
fprintf(fid,'\r\n');

for i = J %number of players
    for j = h %number of modes
        fprintf(fid, ['VC_' i '_' j ' = [zeros(' nbK(1:5*str2double(i)-1) ',n) EK_' i ''' zeros(' nbK(1:5*str2double(i)-1) ',nw); '...
            'Ey_' i '''*Cy' j ' zeros(' nby(1:4*str2double(i)-1) ',nK) Ey_' i '''*Dyw' j '];\r\n']);
    end
end

for i = J
    for j = h
        fprintf(fid, ['VC_' i '_' j '_n = null(VC_' i '_' j '); ']);
    end
    fprintf(fid, '\r\n');
end

for j = h
    fprintf(fid, ['VC_0_' j '_n = eye(n+nK+nw); ']);
end
fprintf(fid, '\r\n');

%%
fprintf(fid,'\r\n');

fprintf(fid, ['eK_' J(1) ' = [eye(nK_1); zeros(' nbK(5+1:end) ',' nbK(1:5*str2double(J(1))-1) ')];\r\n']);
for i = J(2:end-1)
    i_1 = num2str(str2double(i)-1);
    fprintf(fid, ['eK_' i ' = [zeros(' nbK(1:5*str2double(i_1)-1) ',' nbK(5*str2double(i_1)+1:5*str2double(i)-1) ...
        '); eye(' nbK(5*str2double(i_1)+1:5*str2double(i)-1) '); zeros(' nbK(5*str2double(i)+1:end) ',' nbK(5*str2double(i_1)+1:5*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['eK_' J(end) ' = [zeros(' nbK(1:5*str2double(J(end-1))-1) ',' nbK(5*str2double(J(end-1))+1:end) '); eye(' nbK(5*str2double(J(end-1))+1:end) ')];\r\n']);

fprintf(fid, ['ey_' J(1) ' = [eye(ny1); zeros(' nby(4*str2double(J(1))+1:end) ',' nby(1:4*str2double(J(1))-1) ')];\r\n']);
for i = J(2:end-1)
    i_1 = num2str(str2double(i)-1);
    fprintf(fid, ['ey_' i ' = [zeros(' nby(1:4*str2double(i_1)-1) ',' nby(4*str2double(i_1)+1:4*str2double(i)-1) ...
        '); eye(' nby(4*str2double(i_1)+1:4*str2double(i)-1) '); zeros(' nby(4*str2double(i)+1:end) ',' nby(4*str2double(i_1)+1:4*str2double(i)-1) ')];\r\n']);
end
fprintf(fid, ['ey_' J(end) ' = [zeros(' nby(1:4*str2double(J(end-1))-1) ',' nby(4*str2double(J(end-1))+1:end) '); eye(' nby(4*str2double(J(end-1))+1:end) ')];\r\n']);

for i = J %number of players
    i_1 = str2double(i)-1;
    for j = h %number of modes
        fprintf(fid, ['VC_' i '_' j '_t = [zeros(' nbK(5*i_1+1:5*str2double(i)-1) ',n) eK_' i ''' zeros(' nbK(5*i_1+1:5*str2double(i)-1) ',nw); '...
            'ey_' i '''*Cy' j ' zeros(' nby(4*i_1+1:4*str2double(i)-1) ',nK) ey_' i '''*Dyw' j '];\r\n']);
    end
end

fprintf(fid, '\r\n');

UQV = cell(1,p);
VQU = cell(1,p);
ii = 0;
for i = flip(J)
    i_1 = str2double(i)-1;
    for j = 1:p
        if i == J(end)
            UQV{j} = ['UC_' i '_' Phi_star(j,:) '''*QK_' i '_' Psi(j,:) '_t*VC_' i '_' Phi_star(j,:) '_t'];
            VQU{j} = ['VC_' i '_' Phi_star(j,:) '_t''*QK_' i '_' Psi(j,:) '_t''*UC_' i '_' Phi_star(j,:)];
        else
            UQV{j} = ['UC_' i '_' Phi_star(j,:) '''*QK_' i '_' Psi(j,:) '_t*VC_' i '_' Phi_star(j,:) '_t+' UQV{j}];
            VQU{j} = ['VC_' i '_' Phi_star(j,:) '_t''*QK_' i '_' Psi(j,:) '_t''*UC_' i '_' Phi_star(j,:) '+' VQU{j}];
        end
        fprintf(fid, 'cvx_begin sdp quiet\r\n\tcvx_solver Mosek\r\n');
        fprintf(fid, ['\tvariable QK_' i '_' Psi(j,:) '_t(' nbK(5*i_1+1:end) '+' nbu(4*i_1+1:end)...
            ',' nbK(5*i_1+1:5*str2double(i)-1) '+' nby(4*i_1+1:4*str2double(i)-1)  ')\r\n']);
        fprintf(fid, ['\tB' num2str(j+ii*p) ' = [blkdiag(YC_' Phi_ubar(j,:) ',eye(nz)) (R_' Phi_star(j,:) '+' UQV{j} ')*VC_' num2str(i_1) '_' Phi_star(j,:) '_n; ' ...
            'VC_' num2str(i_1) '_' Phi_star(j,:) '_n''*(' VQU{j} '+ R_' Phi_star(j,:) ''') VC_' num2str(i_1) '_' Phi_star(j,:) '_n'''...
            '*blkdiag(XC_' Phi_lbar(j,:) ', eye(nw))*VC_' num2str(i_1) '_' Phi_star(j,:) '_n];\r\n']);
        fprintf(fid, ['\t(B' num2str(j+ii*p) '+B' num2str(j+ii*p) ''')/2 >= 0;\r\n']);
        fprintf(fid, ['\t%%minimize(norm(QK_' i '_' Psi(j,:) '_t));\r\n']);
        fprintf(fid, 'cvx_end\r\n');
        fprintf(fid, 'cvx_status\r\n\r\n');
    end
    ii = ii+1;
end
%% Alternative way to calc
% for i = 1:p
%     fprintf(fid, ['R_' Psi(i,:) ' = [A' Phi_star(i,:) ' zeros(n,nK) Bw' Phi_star(i,:) '; zeros(nK,n) zeros(nK,nK) zeros(nK,nw); Cz' Phi_star(i,:) ' zeros(nz,nK) Dzw' Phi_star(i,:) '];\r\n']);
% end
% 
% for i=1:p
%     fprintf(fid, ['UC_' Psi(i,:) ' = [zeros(nK,n) eye(nK,nK) zeros(nK,nz); Bu' Phi_star(i,:) '''/sqrt(gamm(' num2str(i) ')) zeros(nu,nK) Dzu' Phi_star(i,:) '''/gamm(' num2str(i) ')];\r\n']);
% end
% for i=1:p
%     fprintf(fid, ['VC_' Psi(i,:) ' = [zeros(nK,n) eye(nK,nK) zeros(nK,nw); Cy' Phi_star(i,:) '/sqrt(gamm(' num2str(i) ')) zeros(ny,nK) Dyw' Phi_star(i,:) '/gamm(' num2str(i) ')];\r\n']);
% end
% 
% fprintf(fid, '\r\n');
% 
% for i=1:p
%     fprintf(fid, 'cvx_begin sdp quiet \r\n\tcvx_solver Mosek\r\n');
%     for ii = 1:M
%         fprintf(fid, ['\tvariable AK' num2str(ii) '(nK_' num2str(ii) ',' nbK(1:5*ii-1) ')\r\n']);
%         fprintf(fid, ['\tvariable BK' num2str(ii) '(nK_' num2str(ii) ',' nby(1:4*ii-1) ')\r\n']);
%         fprintf(fid, ['\tvariable CK' num2str(ii) '(nu' num2str(ii) ',' nbK(1:5*ii-1) ')\r\n']);
%         fprintf(fid, ['\tvariable DK' num2str(ii) '(nu' num2str(ii) ',' nby(1:4*ii-1) ')\r\n']);
%     end
% %     fprintf(fid, ['\tQK_' Psi(i,:) ' =  [AK11 zeros(nK_1,' nbK(5*str2double(2)+1:end) ') BK11 zeros(nK_1,ny2); AK21 AK22 BK21 BK22; CK11 zeros(nu1,nK_2) DK11 zeros(nu1,ny2); CK21 CK22 DK21 DK22];\r\n']);
%     fprintf(fid, ['\tQK_' Psi(i,:) ' =  [']);
%     for ii = J(1:end-1)
%         fprintf(fid, ['AK' ii ' zeros(nK_' ii ',' nbK(5*str2double(ii)+1:end) ') BK' ii ' zeros(nK_' ii ',' nby(4*str2double(ii)+1:end) '); ']);
%     end
%     fprintf(fid, ['AK' J(end) ' BK' J(end) ';\r\n\t\t']);
%     for ii = J(1:end-1)
%         fprintf(fid, ['CK' ii ' zeros(nu' ii ',' nbK(5*str2double(ii)+1:end) ') DK' ii ' zeros(nu' ii ',' nby(4*str2double(ii)+1:end) '); ']);
%     end
%     fprintf(fid, ['CK' J(end) ' DK' J(end)]);
%     fprintf(fid, '];\r\n');
%     fprintf(fid, ['\t[blkdiag(YC_' Phi_ubar(i,:) ',eye(nz)) R_' Psi(i,:) '+UC_' Psi(i,:) '''*QK_' Psi(i,:) '*VC_' Psi(i,:) '; VC_' Psi(i,:) '''*QK_' Psi(i,:) '''*UC_' Psi(i,:) '+R_' Psi(i,:) ''' blkdiag(XC_' Phi_lbar(i,:) ',eye(nw))] >= 0;\r\n']);
%     fprintf(fid, 'cvx_end\r\n');
%     fprintf(fid, 'cvx_status\r\n');
%     fprintf(fid, ['QK_' Psi(i,:) ' = full(QK_' Psi(i,:) ');\r\n\r\n']);
% end
% 
% fprintf(fid, '\r\n');


%%
fprintf(fid,'\r\n%%%% \r\n \r\n');

for j = 1:p
    fprintf(fid, ['QK_' Psi(j,:) ' = ']);
    for i = J
        i_1 = num2str(str2double(i)-1);
        fprintf(fid, ['blkdiag(EKbar_' i_1 ',Eubar_' i_1 ')*QK_' i '_' Psi(j,:) '_t*blkdiag(eK_' i ',ey_' i ')''']);
        if i ~= J(end)
            fprintf(fid,'+');
        end
    end
    fprintf(fid, ';\r\n');
end
% for j = 1:p
%     fprintf(fid, ['QK_' Psi(j,:) '(abs(QK_' Psi(j,:) ')<1e-2) = 0;\r\n']);
% end
% fprintf(fid, 'AK = QK(1:nK,1:nK); BK = QK(1:nK,nK+1:end)/gamm;\r\n');
% fprintf(fid, 'CK = QK(nK+1:end,1:nK); DK = QK(nK+1:end,nK+1:end)/gamm;\r\n');
% fprintf(fid, '\r\n');

for j = 1:p
    fprintf(fid, ['AK(:,:,' num2str(j) ') = QK_' Psi(j,:) '(1:nK,1:nK); BK(:,:,' num2str(j) ') = QK_' Psi(j,:) '(1:nK,nK+1:end)/gamm; CK(:,:,' num2str(j) ') = QK_' Psi(j,:) '(nK+1:end,1:nK); DK(:,:,' num2str(j) ') = QK_' Psi(j,:) '(nK+1:end,nK+1:end)/gamm;\r\n']);
end

% fprintf(fid, 'AK(abs(AK)<1e-8) = 0; BK(abs(BK)<1e-8) = 0; CK(abs(CK)<1e-8) = 0; DK(abs(DK)<1e-8) = 0;\r\n');
% fprintf(fid, ['save(''dec_ctrl.mat'', ''AK'',''BK'',''CK'',''DK'');']);

fprintf(fid, '\r\n');
for j = J
    i = str2double(j);
    i_1 = str2double(j)-1;
    fprintf(fid, ['AK' j ' = AK(1:' nbK(1:5*i-1) ',1:' nbK(1:5*i-1) ',:); BK' j ' = BK(1:' nbK(1:5*i-1) ',1:' nby(1:4*i-1) ',:); CK' j ' = CK(' nbu(1:4*i_1-1) '+1:' nbu(1:4*i-1) ',1:' nbK(1:5*i-1) ',:); DK' j ' = DK(' nbu(1:4*i_1-1) '+1:' nbu(1:4*i-1) ',1:' nby(1:4*i-1) ',:);\r\n']);
end

fprintf(fid, '\r\n');
for i = 1:p
    fprintf(fid, 'all(abs(eig(AK(:,:,%d)))<1)\r\n', i);
end

fprintf(fid, '\r\nPsi = [');
for i = 1:p
    fprintf(fid, [Psi(i,:) ' ']);
end
fprintf(fid,'];\r\n');

fclose(fid);
disp('Done');