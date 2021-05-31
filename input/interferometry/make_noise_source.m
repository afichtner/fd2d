%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- number of noise sources
n_noise_sources=1;
sources_everywhere = true; %false

%- characteristics of the noise spectrum ----------------------------------
%- only needed in this routine --------------------------------------------
f_peak=[0.007,0.0085];       % peak frequency in Hz
bandwidth=[0.001,0.005];    % bandwidth in Hz

%- Geographic distribution of sources -------------------------------------
%- Location and width of a Gaussian 'blob' --------------------------------
x_sourcem=[2.6e6,0.4e6];
z_sourcem=[7.5e6,1.1e6];
sourcearea_width=[0.4e6,0.3e6];

%- Or how to set sources everywhere? --------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get model setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path(path,'..');
path(path,'../../code/propagation/');

input_parameters;
input_interferometry;
[X,Z,dx,dz]=define_computational_domain(Lx,Lz,nx,nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% source spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_spectrum=zeros(length(f_sample),n_noise_sources);


%- spectrum for source regions --------------------------------------------
for i=1:n_noise_sources
    noise_spectrum(:,i)=exp(-(abs(f_sample)-f_peak(i)).^2/bandwidth(i)^2);

    figure
    set(gca,'FontSize',20)
    hold on

    plot(f_sample,noise_spectrum(:,i),'k');
    xlabel('frequency [Hz]','FontSize',20);
    title('noise power-spectral density for source region ' ,'FontSize',20)


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geographic power-spectral density distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise_source_distribution=zeros(nx,nz,n_noise_sources);

%- if distribution homogeneous
if sources_everywhere == true
    noise_source_distribution = noise_source_distribution +1.;
    figure;
    set(gca,'FontSize',20);
    load cm_psd

    pcolor(X,Z,noise_source_distribution(:,:,i)');
    shading interp
    colormap(cm_psd)
    xlabel('x [m]','FontSize',20);
    ylabel('z [m]','FontSize',20);
    title('power-spectral density distribution of noise sources for region ','FontSize',20);    
    
else

    %- noise source geography for region 1 ------------------------------------

    for i=1:n_noise_sources
    noise_source_distribution(:,:,i)=(exp(-((X-x_sourcem(i)).^2+(Z-z_sourcem(i)).^2)/(sourcearea_width(i))^2))';

    figure;
    set(gca,'FontSize',20);
    load cm_psd

    pcolor(X,Z,noise_source_distribution(:,:,i)');
    shading interp
    colormap(cm_psd)
    xlabel('x [m]','FontSize',20);
    ylabel('z [m]','FontSize',20);
    title('power-spectral density distribution of noise sources for region ','FontSize',20);    
    end
end


