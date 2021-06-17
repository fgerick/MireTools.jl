
#create a latitude longitude grid for plotting
function latlongrid(nlon,nlat,delta=0.0)
    repeat(range(0,2pi,length=nlon),1,nlat),
    repeat(range(-pi/2+delta,pi/2-delta,length=nlat)',nlon,1)
 end
 
 #plot some data on the a projection map
 function plotsurface(fig,lons,lats,data,ncontours=20; mapprojection=cartopy.crs.Mollweide(), extr = maximum(abs.(data)), kwargs...)
    ax = fig.add_subplot(1, 1, 1, projection=mapprojection)
 
    c=ax.contourf(lons*180/pi, lats*180/pi, data,ncontours,
                transform=cartopy.crs.PlateCarree(),
                norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
 
    ax.coastlines()
    ax.gridlines()
    ax.set_global()
    return ax,c
 end
 
 # make surface velocity plot from polynomial vector field u, using streamlines and
 # uᵩ as colormap in background.
 function vel2surfaceplot(u,lons,lats,fig=figure(); densityv=1.4, velfac=2.0, normvelocity = 0, nlayers=30, kwargs...)
     uθ = map((phi,theta)->real(dot(u,[x*z,y*z,-x^2-y^2])(x=>cos(phi)*sin(theta),
         y=>sin(phi)*sin(theta),z=>cos(theta))/sin(theta)), lons,lats.+2eps().+pi/2)
     uϕ = map((phi,theta)->real(dot(u,[-y,x,0])(x=>cos(phi)*sin(theta),
         y=>sin(phi)*sin(theta),z=>cos(theta))/sin(theta)), lons,lats.+2eps().+pi/2)
     ax = fig.add_subplot(1, 1, 1; projection=cartopy.crs.Mollweide())
     extr = maximum(abs.(uϕ))
     data = @. sqrt(uϕ^2+uθ^2)
     if normvelocity>0
        if normvelocity!=1
            data/=normvelocity
            uϕ/=normvelocity
            uθ/=normvelocity
            extr/=normvelocity
        else
            umax = maximum(abs.(data))
            data/=umax
            uϕ/=umax
            uθ/=umax
            extr/=umax
        end
     end

     ax.streamplot(lons*180/pi, lats*180/pi, uϕ.+eps(), uθ.+eps(),transform=cartopy.crs.PlateCarree(), density=densityv,linewidth=velfac*data/maximum(data), color="k")
 
     c=ax.contourf(lons*180/pi, lats*180/pi, uϕ,nlayers,
                transform=cartopy.crs.PlateCarree(),
                norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
 
     ax.gridlines()
     ax.set_global()
     return fig,ax,c
 end
 
 # make surface velocity plot from polynomial vector field u without streamlines
 function vel2surfaceplotnostream(u,lons,lats,fig=figure();densityv=1.4,velfac=2.0, kwargs...)
     uθ = map((phi,theta)->real(dot(u,[x*z,y*z,-x^2-y^2])(x=>cos(phi)*sin(theta),
         y=>sin(phi)*sin(theta),z=>cos(theta))/sin(theta)), lons,lats.+2eps().+pi/2)
     uϕ = map((phi,theta)->real(dot(u,[-y,x,0])(x=>cos(phi)*sin(theta),
         y=>sin(phi)*sin(theta),z=>cos(theta))/sin(theta)), lons,lats.+2eps().+pi/2)
     ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.Mollweide())
     extr = maximum(abs.(uϕ))
     data = @. sqrt(uϕ^2+uθ^2)
 
     c=ax.contourf(lons*180/pi, lats*180/pi, uϕ,30,
                transform=cartopy.crs.PlateCarree(),
                norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
 
 
     ax.gridlines()
     ax.set_global()
     return fig,ax,c
 end
 
 #plot radial magnetic field component with coastlines in the plot
 function br2surfaceplot(b,lons,lats,fig=figure();nlayers=40, normbr=false, kwargs...)
 
     br=dot(b,[x,y,z])
     bsurf = map((phi,theta)->real(br(x=>cos(phi)*sin(theta),
             y=>sin(phi)*sin(theta),z=>cos(theta))), lons,lats.+pi/2)
     extr=maximum(abs.(real.(bsurf)))
     ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.Mollweide())
    
    if normbr
        bmax=maximum(abs.(bsurf))
        println(bmax)
        bsurf/=bmax
        extr/=bmax
    end

     c=ax.contourf(lons*180/pi, lats*180/pi, bsurf,nlayers,
                transform=cartopy.crs.PlateCarree(),
                norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
 
     ax.coastlines()
     ax.gridlines()
 
     ax.set_global()
     return ax,c
 end
 function plot_br_meridional(v1, ngrid=40;nlayers=40, kwargs...)
     phi_ = range(0,2pi,length=ngrid)
     r_ = range(0.0001,0.9999,length=ngrid)
     X=r_.*cos.(phi_)'
     Y=r_.*sin.(phi_)'
 
     uphi = dot(v1,[x,y,z])
     uphi = broadcast((X,Y,r)->uphi(x=>X,y=>0,z=>Y)/r,X,Y,r_)
     extr=maximum(abs.(uphi))
 
     f = figure()
     ax=f.add_subplot(1,1,1,projection="polar")
 
     maxval = maximum(abs.(real.(uphi)))
     ax.contourf(phi_,r_,uphi, nlayers, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...)
     # plot(range(0,2pi,length=100),ones(100),"k",lw=1)
 
     axis("off")
 end
 function plot_uphi_equator(v1, ngrid=40;nlayers=40, kwargs...)
     phi_ = range(0,2pi,length=ngrid)
     r_ = range(1e-4,1,length=ngrid)
     X=r_.*cos.(phi_)'
     Y=r_.*sin.(phi_)'
 
     uphi = dot(v1,[-y,x,0])
     uphi = broadcast((X,Y)->uphi(x=>X,y=>Y,z=>0)/sqrt(X^2+Y^2),X,Y)
     extr=maximum(abs.(uphi))
 
     f = figure()
     ax=f.add_subplot(1,1,1,projection="polar")
 
     maxval = maximum(abs.(real.(uphi)))
     ax.contourf(phi_,r_,uphi, nlayers, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...) #, norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
     axis("off")
    # PyPlot.plot(range(0,2pi,length=100),ones(100),"k",lw=1)
 end
 
 function plot_u_equator(v1, ngrid=40;nlayers=40, kwargs...)
     phi_ = range(0,2pi,length=ngrid)
     r_ = range(1e-4,1,length=ngrid)
     X=r_.*cos.(phi_)'
     Y=r_.*sin.(phi_)'
 
     uphi = dot(v1,v1)
     uphi = broadcast((X,Y)->sqrt(abs(real(uphi(x=>X,y=>Y,z=>0)))),X,Y)
     extr=maximum(abs.(uphi))
 
     f = figure()
     ax=f.add_subplot(1,1,1,projection="polar")
 
     maxval = maximum(abs.(real.(uphi)))
     ax.contourf(phi_,r_,uphi, nlayers) #, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...) #, norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
     axis("off")
    # PyPlot.plot(range(0,2pi,length=100),ones(100),"k",lw=1)
 end
 
 function plot_uphi_meridional(v1, ngrid=40;nlayers=40, kwargs...)
     phi_ = range(0,2pi,length=ngrid)
     r_ = range(0.0001,0.9999,length=ngrid)
     X=r_.*cos.(phi_)'
     Y=r_.*sin.(phi_)'
 
     uphi = dot(v1,[-y,x,0])
     uphi = broadcast((X,Y)->uphi(x=>X,y=>0,z=>Y)/abs.(X),X,Y)
     extr=maximum(abs.(uphi))
 
     f = figure()
     ax=f.add_subplot(1,1,1,projection="polar")
 
     maxval = maximum(abs.(real.(uphi)))
     ax.contourf(phi_,r_,uphi, nlayers, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...)
 
     # PyPlot.plot(range(0,2pi,length=100),ones(100),"k",lw=1)
     axis("off")
 end
 
 function plot_uphi_1d(v1, ngrid=400;nlayers=40, kwargs...)
     # phi_ = range(0,2pi,length=ngrid)
     r_ = range(0,0.9999,length=ngrid)
     # X=r_.*cos.(phi_)'
     # Y=r_.*sin.(phi_)'
 
     uphi = dot(v1,[-y,x,0])
     uphi = broadcast((X)->real(uphi(x=>X,y=>0,z=>0))/abs.(X),r_)
     extr=findmax(abs.(uphi[.!isnan.(uphi)]))[2]
     uphi/=uphi[extr]
     # @show uphi
     # f = figure()
     plot(r_,uphi; kwargs...)
     xlabel(L"s")
     ylabel(L"u_\phi")
     # ax=f.add_subplot(1,1,1,projection="polar")
     #
     # maxval = maximum(abs.(real.(uphi)))
     # ax.contourf(phi_,r_,uphi, nlayers, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...)
 
     # PyPlot.plot(range(0,2pi,length=100),ones(100),"k",lw=1)
     # axis("off")
 end
 
 function plot_psi_equator(psi1; ngrid=40,nlayers=40, kwargs...)
 #       X, Y = grid(x_,x_)
     phi_ = range(0,2pi,length=ngrid)
     r_ = range(1e-4,1,length=ngrid)
     X = r_.*cos.(phi_)'
     Y = r_.*sin.(phi_)'
 
     psi2 = broadcast((X,Y,R)->(1-R^2)^(3/2)*psi1(x=>X,y=>Y),X,Y,r_)
     # outsideellipse=(X.^2+Y.^2).>=1.0
     # psi2[outsideellipse].=0.
     extr=maximum(abs.(psi2))
 
     f = figure()
     ax=f.add_subplot(1,1,1,projection="polar")
     # maxval = maximum(abs.(real.(psi2)))
     ax.contourf(phi_,r_,real.(psi2), nlayers; kwargs...)  #, norm = PyPlot.matplotlib.colors.Normalize(vmin=-maxval,vmax=maxval); kwargs...) #, norm = PyPlot.matplotlib.colors.Normalize(vmin=-extr,vmax=extr); kwargs...)
     axis("off")
 
     # ax=f.add_subplot(1,1,1)
     # ax.contourf(X,Y,imag.(psi2))
    # PyPlot.plot(range(0,2pi,length=100),ones(100),"k",lw=1)
 #     return x_,uphi,maxval
 end
 
 
 function plot_velocity_equator(a,b,v1; ngrid=50, thickfac=2.5,kwargs...)
     function _grid(x,y)
         X = [i for j in y, i in x]
         Y = [j for j in y, i in x]
         return X,Y
     end
 
 
     X, Y = _grid(range(-a,stop=a,length=ngrid),range(-b,stop=b,length=ngrid))
     ux =  real.([v1[1](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
     uy =  real.([v1[2](Mire.x=>xt,Mire.y=>yt,Mire.z=>0) for (xt,yt) in zip(X,Y)])
     radius = .√(X.^2/a^2+Y.^2/b^2)
     u = .√(ux.^2+uy.^2)
     phi = atan.(Y/b,X/a);
 
     outsideellipse=(X.^2/a^2+Y.^2/b^2).>=1.0
     ux[outsideellipse].=0.
     uy[outsideellipse].=0.
     u[outsideellipse].=0.
     #
     # streamplot(Float64.(X), Float64.(Y), Float64.(ux),
     #         Float64.(uy) ; color=Float64.(u), kwargs...)
 
                 streamplot(Float64.(X), Float64.(Y), Float64.(ux),
                         Float64.(uy) ; kwargs...) # linewidth=thickfac*Float64.(abs.(u)/maximum(abs.(u))), kwargs...)
 
 
     ellipsex=a*cos.(range(0,stop=2π,length=100))
     ellipsey=b*sin.(range(0,stop=2π,length=100))
     PyPlot.plot(ellipsex,ellipsey,"k",linewidth=2)
     PyPlot.axis("equal")
     PyPlot.axis("off");
 end
 