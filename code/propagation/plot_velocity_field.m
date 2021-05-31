if (mod(n,4)==0)
    
    hold on
    
    %- plot velocity field ------------------------------------------------
   
    pcolor(X,Z,v');
    set(gca,'FontSize',20);
    axis image
    
    %- scale, label, etc ... ----------------------------------------------
    
    if (n<0.8*length(t))
        scale=max(max(abs(v)));
    end
   
    caxis([-scale scale]);
    colormap(cm);
    shading interp
    xlabel('x [m]','FontSize',20);
    ylabel('z [m]','FontSize',20);
    
    if (strcmp(simulation_mode,'forward') || strcmp(simulation_mode,'forward_green'))
        title('velocity field [m/s]','FontSize',20);
    elseif (strcmp(simulation_mode,'correlation') && t(n)<0)
        title('acausal correlation field','FontSize',20);
    elseif (strcmp(simulation_mode,'correlation') && t(n)>=0)
        title('causal correlation field','FontSize',20);
    end
    
    %- plot source and receiver positions ---------------------------------
    
    for p=1:5
    
    for k=1:length(src_x)
        plot(src_x(k),src_z(k),'kx')
    end
    
    for k=1:length(rec_x)
        plot(rec_x(k),rec_z(k),'ko')
    end
    
    end
       
    %- record movie -------------------------------------------------------
    
    if strcmp(make_movie,'yes')
    
        if exist('movie_index','var')
            movie_index=movie_index+1;
        else
            movie_index=1;
        end
        
        M(movie_index)=getframe(gcf);
        
    end
    
    %- finish -------------------------------------------------------------
    
    pause(0.01)
    
    hold off
    clf
            
end