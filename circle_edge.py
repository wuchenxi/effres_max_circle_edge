#Usage: python3 circle_edge.py

from random import *

n=20
#0-1: on coordinate on circle, 2: on tip
points_conf=[2.]*n

#Calculate pairwise sum of effective resistence
def effres(pc):
    b=0
    y=[]
    for i in range(n):
        if pc[i]>1:
            b+=1
        else:
            y+=[pc[i]]            
    r=float(b*(n-b))
    for i in range(len(y)):
        r+=y[i]*(1-y[i])
        for j in range(i, len(y)):
            r+=abs(y[i]-y[j])*(1-abs(y[i]-y[j]))*b
    return r

#Move the i-th point around to maximize effective resistence
def optimize(pc, i):
    pc1=[a for a in pc]
    pc1[i]=2
    r0=effres(pc1)
    grid=[0]
    for j in range(n):
        if pc[j]<1 and j!=i:
            grid+=[pc[j]]
    grid+=[1]
    m=len(grid)-1
    rv=[]
    pos=[]
    for j in range(m):
        pc1[i]=grid[j]
        rl=effres(pc1)
        pc1[i]=0.5*(grid[j]+grid[j+1])
        rm=effres(pc1)
        pc1[i]=grid[j+1]
        rr=effres(pc1)
        posj=quadratic_opt(rl, rm, rr, grid[j], grid[j+1])
        pc1[i]=posj
        pos+=[posj]
        rv+=[effres(pc1)]
    rv+=[r0]
    pos+=[2]
    tuples=[[-rv[j], pos[j]] for j in range(m+1)]
    tuples.sort()
    pc[i]=tuples[0][1]
    pc.sort()
    return effres(pc)

#quadratic interpolation and optimization
def quadratic_opt(l, m, r, lx, rx):
    slope1=m-l
    slope2=r-m
    slopel=(3*slope1-slope2)/2
    sloper=(3*slope2-slope1)/2
    if slopel*sloper>=0:
        if slopel+sloper>0:
            return rx
        else:
            return lx
    else:
        wr=-slopel/(sloper-slopel)
        wl=sloper/(sloper-slopel)
        return lx*wl+rx*wr

#print(format(effres(points_conf), ".6f"))
for i in range(10000):
    id=randrange(n)
    r=optimize(points_conf, id)
    #print(format(r, ".6f"))
print(','.join([format(a, ".6f") for a in points_conf]))
