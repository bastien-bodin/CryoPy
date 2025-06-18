import numpy as np

def Cuve2D(x0:float, y0:float, Lx:float, Ly:float, dp:float, nl:int) -> tuple:
    if nl <= 0:
        raise ValueError
    dp2 = dp*0.5
    eps = dp*0.1
    x = []
    y = []
    for i in range(nl):
        i_mod = i%2
        yl = np.arange(y0-dp2*i, y0+Ly+dp2*i_mod+eps, dp)
        xl = np.ones_like(yl)*(x0-dp2*i)
        xc = np.arange(x0+dp-i*dp2, x0+Lx-dp+dp2*i+eps, dp)
        yc = np.ones_like(xc)*(y0-i*dp2)
        yr = yl
        xr = np.ones_like(yr)*(x0+Lx+i*dp2)
        xi = np.concatenate([xl, xc, xr])
        yi = np.concatenate([yl, yc, yr])
        x.append(xi)
        y.append(yi)
    X = np.concatenate([x[i] for i in range(nl)])
    Y = np.concatenate([y[i] for i in range(nl)])
    return X, Y

def MyCuve2D(x0:float, y0:float, Lx:float, Ly:float, dp:float, nl:int) -> tuple:
    dp2 = dp/2.0
    x, y = Cuve2D(x0-dp2,y0-dp2, Lx+dp, Ly+dp2, dp, nl)
    return x,y

def Block2D(x0:float, y0:float, Lx:float, Ly:float, dp:float) -> tuple:
    eps = dp*0.1
    x, y = np.mgrid[x0:x0+Lx+eps:dp,
                    y0:y0+Ly+eps:dp]
    x = x.ravel(); y = y.ravel()
    return x, y

def MyBlock2D(x0:float, y0:float, Lx:float, Ly:float, dp:float) -> tuple:
    dp2 = 0.5*dp
    return Block2D(x0+dp2, y0+dp2, Lx-dp, Ly-dp, dp)

def MyPit2D(x0:float, y0:float, Lx:float, Ly:float,
            lx_pit:float, ly_pit:float, dp:float, nl:int) -> tuple:
    eps = dp*0.1
    x1, y1 = MyCuve2D(x0, y0, Lx, Ly, dp, nl)
    
    x0_pit = x0 + 0.5*(Lx-lx_pit)
    y0_pit = y0 - ly_pit
    x2, y2 = MyCuve2D(x0_pit, y0_pit, lx_pit, ly_pit, dp, nl)
    
    mask = (x1<x0_pit-dp-eps) + (x1>x0_pit+lx_pit+dp+eps)
    
    x1 = x1[mask]; y1 = y1[mask]
    
    x = np.concatenate([x1, x2]); y = np.concatenate([y1, y2])
    return x, y

def DoomMons(x0:float, y0:float,
             Lx:float, Ly:float,
             lx_vent:float, ly_vent:float,
             dp:float, nl:int=2) -> tuple:
    """Generates a geometry inspired from Doom Mons.

    Args:
        x0 (float): Beginning of the flow front /x.
        y0 (float): Beginning of the flow front /y.
        Lx (float): Length of the domain.
        Ly (float): Height of the domain.
        lx_vent (float): Length of the vent.
        ly_vent (float): Depth of the vent.
        dp (float): Spacing between particles.
        nl (int, optional): Number of layers. Defaults to 2.

    Returns:
        tuple: Coordinates of the substrate particles (Xs, Ys)
               and of the inflow particles (Xf, Yf)
    """
    dp2 = dp*0.5
    X1, Y1 = MyCuve2D(x0-lx_vent-dp*2.0, y0,
                      Lx+dp*2.0, Ly,
                      dp, nl)
    mask = (X1 < x0 - lx_vent - dp*1.1) + (X1 > x0 +dp * 1.1)# *(Y1>y0-dp*1.1)
    X1 = X1[mask]; Y1 = Y1[mask]
    X2, Y2 = MyCuve2D(x0-lx_vent, y0-ly_vent,
                      lx_vent, ly_vent-dp*0.6,
                      dp, nl)
    Xs = np.concatenate([X1, X2])
    Ys = np.concatenate([Y1, Y2])
    Xf, Yf = MyBlock2D(x0-lx_vent, y0-ly_vent,
                       lx_vent, ly_vent, dp)
    return Xs, Ys, Xf, Yf

def DoomMonsReservoir(x0:float, y0:float,
                      Lx:float, Ly:float,
                      lx_vent:float, ly_vent:float,
                      lx_cuve:float, ly_cuve:float,
                      dp:float, nl:int=2) -> tuple:
    eps = 0.1*dp; dp2 = dp*0.5
    Xc, Yc = MyCuve2D(x0 - (lx_cuve+2*lx_vent), y0 - 2*ly_vent,
                      lx_cuve, ly_cuve+2*ly_vent,
                      dp, nl)
    mask = (Yc > y0 - ly_vent+ dp2*nl + eps)*(Xc>x0 - 2.5*lx_vent) + (Xc<x0 - 2.0*lx_vent+eps)
    Xc = Xc[mask]; Yc = Yc[mask]
    
    X1 = np.arange(x0-2.0*lx_vent + dp2, x0+lx_vent+eps, dp)
    Y1 = np.ones_like(X1)*(y0-2*ly_vent-dp/2)
    
    Xp = np.arange(x0-2.0*lx_vent + dp2, x0-dp2+eps, dp)
    Yp = np.ones_like(Xp)*(y0-ly_vent+dp/2)
    
    Yh = np.arange(y0 - ly_vent+ dp2*(nl+(nl%2+1)), y0-dp2*(nl-1)+eps, dp)
    Xh = np.ones_like(Yh)*(x0-dp2)
    
    Yhb = np.arange(y0 - 2*ly_vent- dp2*(nl-(nl%2+1)), y0-dp2*(nl-1)+eps, dp)
    Xhb = np.ones_like(Yhb)*(x0+lx_vent+dp2)
    for i in range(1,nl):
        k = i%2
        X2 = np.arange(x0-2.0*lx_vent + dp2 + dp2*k, x0+lx_vent+eps + dp2*k, dp)
        Y2 = np.ones_like(X2)*(y0-2*ly_vent-dp/2 - dp2*i)
        
        Xp2 = np.arange(x0-2.0*lx_vent + dp2 + dp2*k, x0-dp2+eps- dp2*k, dp)
        Yp2 = np.ones_like(Xp2)*(y0-ly_vent+dp/2 + dp2*i)
        
        X1 = np.concatenate([X1, X2])
        Y1 = np.concatenate([Y1, Y2])
        
        Yh2 = np.arange(y0 - ly_vent+ dp2*(nl+(nl%2+1)-k),
                        y0-dp2*(nl-1-k)+eps,
                        dp)
        Xh2 = np.ones_like(Yh2)*(x0-dp2*(i+1))
        
        Yhb2 = np.arange(y0 - 2*ly_vent - dp2*(nl-(nl%2+1)+k), y0-dp2*(nl-1-k)+eps, dp)
        Xhb2 = np.ones_like(Yhb2)*(x0+lx_vent+dp2*(1+i))
        
        Xp = np.concatenate([Xp, Xp2])
        Yp = np.concatenate([Yp, Yp2])
        
        Xh = np.concatenate([Xh, Xh2])
        Yh = np.concatenate([Yh, Yh2])
        
        Xhb = np.concatenate([Xhb, Xhb2])
        Yhb = np.concatenate([Yhb, Yhb2])
    Xh = np.concatenate([Xh, Xhb])
    Yh = np.concatenate([Yh, Yhb])
    X1 = np.concatenate([X1, Xp])
    Y1 = np.concatenate([Y1, Yp])
    
    Xr = np.concatenate([Xc, X1, Xh])
    Yr = np.concatenate([Yc, Y1, Yh])
    
    # Domain
    Xs, Ys = MyCuve2D(x0-dp*2.0, y0,
                      Lx+dp*2.0, Ly,
                      dp, nl)
    mask = (Xs < x0 - dp*1.1) + (Xs > x0 +lx_vent+dp * 1.1)# *(Y1>y0-dp*1.1)
    Xs = Xs[mask]; Ys = Ys[mask]
    
    # Fluid
    Xf1, Yf1 = MyBlock2D(x0-(lx_cuve+2*lx_vent), y0-2*ly_vent,
                         lx_cuve, 2*ly_vent, dp)
    
    Xf2, Yf2 = MyBlock2D(x0-(lx_cuve+2*lx_vent), y0,
                         lx_cuve, ly_cuve, dp)
    Xf3, Yf3 = MyBlock2D(x0-2*lx_vent, y0-2*ly_vent,
                         2*lx_vent, ly_vent, dp)
    Xf4, Yf4 = MyBlock2D(x0, y0-2*ly_vent,
                         lx_vent, 2*ly_vent, dp)
    Xf = np.concatenate([Xf1, Xf2, Xf3, Xf4])
    Yf = np.concatenate([Yf1, Yf2, Yf3, Yf4])
    return Xr, Yr, Xs, Ys, Xf, Yf
def MyFissure2D(x0:float, y0:float, Lx:float, Ly:float,
                lx_pit:float, ly_pit:float, dp:float, nl:int) -> tuple:
    eps = dp*0.1; dp2 = dp*0.5
    x1, y1 = MyCuve2D(x0-lx_pit-nl*dp2, y0, Lx, Ly, dp, nl)
    
    x0_pit = x0 - lx_pit
    y0_pit = y0 - ly_pit
    x2, y2 = MyCuve2D(x0_pit, y0_pit, lx_pit, ly_pit, dp, nl)
    
    mask = (x1<x0_pit-eps-nl*dp2) + (x1>x0_pit+lx_pit+nl*dp2+eps)
    
    x1 = x1[mask]; y1 = y1[mask]
    
    x = np.concatenate([x1, x2]); y = np.concatenate([y1, y2])
    return x, y

if __name__=='__main__':
    Xb, Yb, Xf, Yf = DoomMons(0,0,200e3,1150,20e3,50,5)
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.plot(Xb, Yb, '.',
             Xf, Yf, '.')
    plt.show()