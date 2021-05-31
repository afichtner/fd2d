function [X,Z,K_rho_0,K_mu_0,K_rho,K_beta]=run_adjoint

%==========================================================================
% run forward simulation
%
% output:
%--------
% K: sensitivity kernel with respect to density
%==========================================================================

%==========================================================================
% read input and set paths
%==========================================================================

path(path,'propagation/');
path(path,'../input/');

load cm_velocity;
load('../output/v_forward.mat');
load('../output/dx_u_forward.mat');
load('../output/dz_u_forward.mat');

input_parameters;
nt=5*round(nt/5);

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

[mu,rho]=define_material_parameters(nx,nz,model_type);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
    
%- dynamic fields ---------------------------------------------------------

v=zeros(nx,nz);
dx_u=zeros(nx,nz);
dz_u=zeros(nx,nz);

K_rho_0=zeros(nx,nz);
K_mu_0=zeros(nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

%- read adjoint source locations ------------------------------------------

fid=fopen([adjoint_source_path '/source_locations'],'r');
src_x=zeros(1);
src_z=zeros(1);

k=1;
while (feof(fid)==0)
    src_x(k)=fscanf(fid,'%g',1);
    src_z(k)=fscanf(fid,'%g',1);
    fgetl(fid);
    k=k+1;
end

fclose(fid);

%- read adjoint source time functions -------------------------------------

ns=length(src_x);
stf=zeros(ns,nt);

for n=1:ns
    fid=fopen([adjoint_source_path '/src_' num2str(n)],'r');
    stf(n,1:nt)=fscanf(fid,'%g',nt);
    stf(n,nt:-1:1)=stf(n,1:nt);
end

%- compute indices for adjoint source locations ---------------------------

src_x_id=zeros(1,ns);
src_z_id=zeros(1,ns);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:ns

    src_x_id(i)=min(find(min(abs(x-src_x(i)))==abs(x-src_x(i))));
    src_z_id(i)=min(find(min(abs(z-src_z(i)))==abs(z-src_z(i))));

end

%==========================================================================
% initialise absorbing boundaries.
%==========================================================================

absbound=ones(nx,nz);
init_absbound;
   
%==========================================================================
% iterate
%==========================================================================

for n=1:nt
    
    %- compute divergence of current stress tensor and add external forces
    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    for i=1:ns
        DS(src_x_id(i),src_z_id(i))=DS(src_x_id(i),src_z_id(i))+stf(i,n);
    end
    
    %- update velocity field ----------------------------------------------
    
    v=v+dt*DS./rho;
    
    %- apply absorbing boundary taper -------------------------------------
    
    v=v.*absbound;
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    sxy=sxy+dt*mu(1:nx-1,:).*dx_v(v,dx,dz,nx,nz,order);
    szy=szy+dt*mu(:,1:nz-1).*dz_v(v,dx,dz,nx,nz,order);
    
    %- compute strain field from the stress field -------------------------
    
    dx_u(1:nx-1,:)=sxy./mu(1:nx-1,:);
    dz_u(:,1:nz-1)=szy./mu(:,1:nz-1);
    
    %- plot and compute kernel every 10th iteration -----------------------
    
    if (mod(n+4,5)==0)
        
        %- adjoint field --------------------------------------------------
        subplot(2,2,1);
        pcolor(X,Z,v');
        hold on
        for k=1:ns
            plot(src_x(k),src_z(k),'kx');
        end
        hold off
        
        caxis([-0.8*max(max(abs(v))) 0.8*max(max(abs(v)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('adjoint velocity field [m/s]');
        
        %- forward field --------------------------------------------------
        v_forward_s=squeeze(v_forward((n+4)/5,:,:));
        dx_u_forward_s=squeeze(dx_u_forward((n+4)/5,:,:));
        dz_u_forward_s=squeeze(dz_u_forward((n+4)/5,:,:));
        
        subplot(2,2,2);
        pcolor(X,Z,v_forward_s');
        
        caxis([-0.8*max(max(abs(v_forward_s))) 0.8*max(max(abs(v_forward_s)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('forward velocity field [m/s]');
        
        %- interaction ----------------------------------------------------
        interaction_rho=v.*v_forward_s;
        interaction_mu=2*(dx_u.*dx_u_forward_s+dz_u.*dz_u_forward_s);
         
        subplot(2,2,3);
        pcolor(X,Z,interaction_rho');
        
        caxis([-0.8*max(max(abs(interaction_rho))) 0.8*max(max(abs(interaction_rho)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('interaction (forward \cdot adjoint)');
        
        %- kernel ---------------------------------------------------------
        K_rho_0=K_rho_0+interaction_rho*dt;
        K_mu_0=K_mu_0+interaction_mu*dt;
         
        subplot(2,2,4);
        pcolor(X,Z,K_rho_0');
        
        caxis([-0.5*max(max(abs(K_rho_0))) 0.5*max(max(abs(K_rho_0)))]);
        colormap(cm);
        axis image
        shading interp
        xlabel('x [m]');
        ylabel('z [m]');
        title('sensitivity kernel');
        
        pause(0.01)
        
    end
    
end

%==========================================================================
% compute relative and derived relative kernels kernels
%==========================================================================

K_rho=K_rho_0+mu.*K_mu_0./rho;
K_beta=2*sqrt(mu.*rho).*K_mu_0;

K_rho=rho.*K_rho;
K_beta=sqrt(mu./rho).*K_beta;

K_rho_0=rho.*K_rho_0;
K_mu_0=mu.*K_mu_0;

