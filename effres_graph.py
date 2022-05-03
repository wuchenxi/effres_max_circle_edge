#usage: python3 effres_graph.py
#format for specifying graphs:
#   [number of vertices, number of edges, list of end points of edges,
#        list of edge lengths]
import scipy.linalg
import numpy
from math import *
from random import *

n=10
points_conf=[0.]*n


g=[2, 2, [(0,1), (1,1)], [1, 1]]
#g=[2, 3, [(0, 1), (0, 1), (0, 1)], [1, 2, 3]]

M=sum(g[3])+1

def coord(x):
    j=(int)(floor(x/M))
    t=x-j*M
    return j, t

#effective resistance calculation
def effres(x, y):
    jx, tx=coord(x)
    jy, ty=coord(y)
    ng=[g[0], g[1], [a for a in g[2]], [a for a in g[3]]]
    splitx=False
    splity=False
    if tx<=0:
        vx=ng[2][jx][0]
    elif tx>=ng[3][jx]:
        vx=ng[2][jx][1]
    else:
        splitx=True
    if ty<=0:
        vy=ng[2][jy][0]
    elif ty>=ng[3][jy]:
        vy=ng[2][jy][1]
    else:
        splity=True
    if splitx and splity:
        ng[0]+=2
        ng[1]+=2
        vx=ng[0]-2
        vy=ng[0]-1
        if jx==jy:
            if tx==ty:
                return 0
            elif tx>ty:
                tmp=ty
                ty=tx
                tx=tmp
            ng[2]+=[(ng[0]-2, ng[0]-1), (ng[0]-1, ng[2][jx][1])]
            ng[2][jx]=(ng[2][jx][0], ng[0]-2)
            ng[3]+=[ty-tx, ng[3][jx]-ty]
            ng[3][jx]=tx
        else:
            ng[2]+=[(ng[0]-2, ng[2][jx][1])]
            ng[2][jx]=(ng[2][jx][0], ng[0]-2)
            ng[3]+=[ng[3][jx]-tx]
            ng[3][jx]=tx
            ng[2]+=[(ng[0]-1, ng[2][jy][1])]
            ng[2][jy]=(ng[2][jy][0], ng[0]-1)
            ng[3]+=[ng[3][jy]-ty]
            ng[3][jy]=ty
    elif splitx:
        ng[0]+=1
        ng[1]+=1
        vx=ng[0]-1
        ng[2]+=[(ng[0]-1, ng[2][jx][1])]
        ng[2][jx]=(ng[2][jx][0], ng[0]-1)
        ng[3]+=[ng[3][jx]-tx]
        ng[3][jx]=tx
    elif splity:
        ng[0]+=1
        ng[1]+=1
        vy=ng[0]-1
        ng[2]+=[(ng[0]-1, ng[2][jy][1])]
        ng[2][jy]=(ng[2][jy][0], ng[0]-1)
        ng[3]+=[ng[3][jy]-ty]
        ng[3][jy]=ty
    bdry_mat=[[0. for i in range(ng[0])] for j in range(ng[1])]
    for i in range(ng[1]):
        bdry_mat[i][ng[2][i][0]]-=1.
        bdry_mat[i][ng[2][i][1]]+=1.
    b=numpy.array(bdry_mat)
    d=numpy.diag([1/r for r in ng[3]])
    lap=numpy.dot(b.transpose(), numpy.dot(d, b))
    b_l=[0 for i in range(ng[0])]
    b_l[vx]=1
    b_l[vy]=-1
    b_l[0]=0
    lap[0, 0]=1.
    for i in range(1, ng[0]):
        lap[0, i]=0
    res=scipy.linalg.solve(lap, numpy.array(b_l))
    return abs(res[vx]-res[vy])

def sum_effres(pc, i, loc):
    r=0
    for j in range(len(pc)):
        if j !=i:
            r+=effres(loc, pc[j])
    return r

#Move the i-th point around to maximize effective resistence
def optimize(pc, i):
    grid=[]
    break_pts=[[] for j in range(g[1])]
    for j in range(len(pc)):
        if j != i:
            k, t=coord(pc[j])
            if t>0 and t<g[3][k]:
                break_pts[k]+=[pc[j]]
    for j in range(g[1]):
        if len(break_pts[j])==0:
            grid+=[(j*M, j*M+g[3][j])]
        else:
            grid+=[(j*M, break_pts[j][0])]
            for k in range(len(break_pts[j])-1):
                grid+=[(break_pts[j][k], break_pts[j][k+1])]
            grid+=[(break_pts[j][-1], j*M+g[3][j])]
    m=len(grid)
    rv=[]
    pos=[]
    for j in range(m):
        rl=sum_effres(pc, i, grid[j][0])
        rm=sum_effres(pc, i, (grid[j][1]+grid[j][0])/2)
        rr=sum_effres(pc, i, grid[j][1])
        posj=quadratic_opt(rl, rm, rr, grid[j][0], grid[j][1])
        pos+=[posj]
        rv+=[sum_effres(pc, i, posj)]
    tuples=[[-rv[j], pos[j]] for j in range(m)]
    tuples.sort()
    pc[i]=tuples[0][1]
    pc.sort()
    return 0

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
for i in range(100):
    id=randrange(n)
    r=optimize(points_conf, id)
    #print(format(r, ".6f"))
for pt in points_conf:
    j, t=coord(pt)
    print("edge "+str(j)+": "+format(t, ".6f"))
    

    

