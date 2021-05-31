function [u,t,rec_x,rec_z]=run_forward

%==========================================================================
% run forward simulation
%
% output:
%--------
% u: displacement seismograms
% t: time axis
% rec_x: x-coordinate of receiver positions
% rec_z: z-coordinate of receiver positions
%
%==========================================================================

%==========================================================================
% set paths and read input
%==========================================================================

path(path,'propagation/');
path(path,'../input/');
path(path,'../input/interferometry');

input_parameters;
nt=5*round(nt/5);

load cm_velocity;

%==========================================================================
% initialise simulation
%==========================================================================

%- material and domain ----------------------------------------------------

[mu,rho]=define_material_parameters(nx,nz,model_type);

[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);
 
plot_model;

%- forward simulations ('forward', 'forward_green') -----------------------

if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_store') || strcmp(simulation_mode,'forward_green'))

    %- time axis ----------------------------------------------------------
    
    t=0:dt:dt*(nt-1);
    
    %- compute indices for source locations -------------------------------

    ns=length(src_x);
    
    src_x_id=zeros(1,ns);
    src_z_id=zeros(1,ns);

    x=0:dx:Lx;
    z=0:dz:Lz;

    for i=1:ns

        src_x_id(i)=min(find(min(abs(x-src_x(i)))==abs(x-src_x(i))));
        src_z_id(i)=min(find(min(abs(z-src_z(i)))==abs(z-src_z(i))));

    end
    
    %- make source time function ------------------------------------------

    make_source_time_function;
    
    %- initialise interferometry ------------------------------------------
    
    if strcmp(simulation_mode,'forward_green')
    
        input_interferometry;
        w_sample=2*pi*f_sample;

        %- Fourier transform of the forward Greens function
        G_2=zeros(nx,nz,length(f_sample));
        
    end
    
%- forward simulation to compute correlation function ---------------------     
    
elseif strcmp(simulation_mode,'correlation')

    %- time axis ----------------------------------------------------------
    t=-(nt-1)*dt:dt:(nt-1)*dt;
    
    %- initialise interferometry ------------------------------------------
    input_interferometry;
    w_sample=2*pi*f_sample;
    
    %- Fourier transform of the correlation velocity field
    C_2=zeros(nx,nz,length(f_sample));
        
    %- load frequency-domain Greens function
    load('../output/interferometry/G_2.mat');
    
    %- initialise noise source locations and spectra
    make_noise_source;
    
end
 
%- dynamic fields and absorbing boundary field ----------------------------

v=zeros(nx,nz);
absbound=ones(nx,nz);

v_forward=zeros(nt/5,nx,nz);
dx_u_forward=zeros(nt/5,nx,nz);
dz_u_forward=zeros(nt/5,nx,nz);

sxy=zeros(nx-1,nz);
szy=zeros(nx,nz-1);

%- compute indices for receiver locations ---------------------------------

n_receivers=length(rec_x);

rec_x_id=zeros(1,n_receivers);
rec_z_id=zeros(1,n_receivers);

x=0:dx:Lx;
z=0:dz:Lz;

for i=1:n_receivers

    rec_x_id(i)=min(find(min(abs(x-rec_x(i)))==abs(x-rec_x(i))));
    rec_z_id(i)=min(find(min(abs(z-rec_z(i)))==abs(z-rec_z(i))));

end

%- initialise seismograms -------------------------------------------------

displacement_seismograms=zeros(n_receivers,nt);
velocity_seismograms=zeros(n_receivers,nt);

%- initialise absorbing boundary taper a la Cerjan ------------------------

init_absbound;

%==========================================================================
% iterate
%==========================================================================

figure;

for n=1:length(t)
    
    %- compute divergence of current stress tensor ------------------------
    
    DS=div_s(sxy,szy,dx,dz,nx,nz,order);
    
    %- add point sources --------------------------------------------------
    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_store') || strcmp(simulation_mode,'forward_green'))
    
        for i=1:ns
            DS(src_x_id(i),src_z_id(i))=DS(src_x_id(i),src_z_id(i))+stf(n);
        end
        
    end
    
    %- add sources of the correlation field -------------------------------
    
    if (strcmp(simulation_mode,'correlation') && (t(n)<0.0))
        
        %- transform on the fly to the time domain
        
        dw=w_sample(2)-w_sample(1);
        S=zeros(nx,nz,n_noise_sources);
        
        i=sqrt(-1);
        for ns=1:n_noise_sources
            
            %- inverse Fourier transform for each noise source region
            for k=2:length(f_sample)
                S(:,:,ns)=S(:,:,ns)+noise_spectrum(k,ns)*conj(G_2(:,:,k))*exp(i*w_sample(k)*t(n));
            end
            S(:,:,ns)=real(dw*S(:,:,ns)/pi);
            
            %- add sources
            DS=DS+noise_source_distribution(:,:,ns).*real(S(:,:,ns));
        end
           
    end
    
    %- update velocity field ----------------------------------------------
    
    v=v+dt*DS./rho;
    
    %- apply absorbing boundary taper -------------------------------------
    
    v=v.*absbound;
    
    %- compute derivatives of current velocity and update stress tensor ---
    
    sxy=sxy+dt*mu(1:nx-1,1:nz).*dx_v(v,dx,dz,nx,nz,order);
    szy=szy+dt*mu(:,1:nz-1).*dz_v(v,dx,dz,nx,nz,order);
    
    %- record velocity seismograms ----------------------------------------
    
    for k=1:n_receivers
        velocity_seismograms(k,n)=v(rec_x_id(k),rec_z_id(k));
    end
    
    %- store time-reversed history ----------------------------------------
    
    if (strcmp(simulation_mode,'forward_store'))
    
        if (mod(n,5)==0)
            v_forward(nt/5+1-n/5,:,:)=v(:,:);
            dx_u_forward(nt/5+1-n/5,1:nx-1,:)=sxy./mu(1:nx-1,:);
            dz_u_forward(nt/5+1-n/5,:,1:nz-1)=szy./mu(:,1:nz-1);
        end
        
    end
    
    %- accumulate Fourier transform of the displacement Greens function ---
    
    if strcmp(simulation_mode,'forward_green')
    
        i=sqrt(-1);
        for k=1:length(w_sample)
            G_2(:,:,k)=G_2(:,:,k)+v(:,:)*exp(-i*w_sample(k)*t(n))*dt;
        end
        
    end
    
    %- accumulate Fourier transform of the correlation velocity field -----
    
    if strcmp(simulation_mode,'correlation')
    
        i=sqrt(-1);
        for k=1:length(w_sample)
            C_2(:,:,k)=C_2(:,:,k)+v(:,:)*exp(-i*w_sample(k)*t(n))*dt;
        end
        
    end
    
    %- plot velocity field every 4th time step ----------------------------
    
    plot_velocity_field;
    
end

%==========================================================================
% output 
%==========================================================================

%- store time-reversed forward field --------------------------------------

if strcmp(simulation_mode,'forward_store')
    save('../output/v_forward','v_forward');
    save('../output/dz_u_forward','dz_u_forward');
    save('../output/dx_u_forward','dx_u_forward');
end

%- store Fourier transformed velocity Greens function ---------------------

if strcmp(simulation_mode,'forward_green')
    save('../output/interferometry/G_2','G_2');
end

%- store Fourier transformed correlation velocity field -------------------

if strcmp(simulation_mode,'correlation')
    save('../output/interferometry/C_2','C_2');
end

%- displacement seismograms -----------------------------------------------

u=cumsum(velocity_seismograms,2)*dt;

%- store the movie if wanted ----------------------------------------------

if strcmp(make_movie,'yes')
    writerObj=VideoWriter(movie_file,'MPEG-4');
    open(writerObj);
    writeVideo(writerObj,M);
    close(writerObj);
end