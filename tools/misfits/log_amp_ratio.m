function [misfit,adstf] = log_amp_ratio(u,u_0,win,t)

% pi= if more precise pi value needed
dt = abs(t(2)-t(1));

% find the opposite window
win2 = fliplr(win);
if win2 == win
    win2 = flipud(win);
end

% compute energy for both windows of both observed and synthetic
u1 = win.*u;
u2 = win2.*u;
e_caus = trapz(u1.^2)*dt;
e_acaus = trapz(u2.^2)*dt;

u01 = win.*u_0;
u02 = win2.*u_0;
e0_caus = trapz(u01.^2)*dt;
e0_acaus = trapz(u02.^2)*dt;

% How to correct when e0_acaus should be 0?


% get dE
de_caus = win.^2 .* u;
% is time integration necessary?
de_acaus = win2.^2 .* u;

adstf = de_caus/e_caus - de_acaus/e_acaus;


if e0_acaus==0
    misfit = 0.5*(log(e_caus/e_acaus))^2;
else
    misfit = 0.5*(log(e_caus/e_acaus)-log(e0_caus/e0_acaus))^2;
end