from sympy import *
import math

def algebras():
    ax,bx,cx,ay,by,cy,t,dx,dy = symbols('ax bx cx ay by cy t dx dy')
    init_printing()

    #https://stackoverflow.com/questions/14264431/expanding-algebraic-powers-in-python-sympy
    #smichr
    def sack(expr):            
        return expr.replace(
            lambda x: x.is_Pow and x.exp > 0,
            lambda x: Symbol('*'.join([x.base.name]*x.exp))
        )

    x = (1-t)*(1-t)*ax+2*(1-t)*t*bx+t*t*cx
    y = (1-t)*(1-t)*ay+2*(1-t)*t*by+t*t*cy

    expr = expand((bx-x)*diff(x,t)+(by-y)*diff(y,t))
    expr = collect(expr,t)
    print("Cubic for B")
    print(sack(expr))

    expr = expand((dx-x)*diff(x,t)+(dy-y)*diff(y,t))
    expr = collect(expr,t)
    print("Cubic for D")
    print(sack(expr))

    expr = expand((bx-x)**2 + (by-y)**2)
    expr = collect(expr,t)
    print("Quartic for B")
    print(sack(expr))

    expr = expand((dx-x)**2 + (dy-y)**2)
    expr = collect(expr,t)
    print("Quartic for D")
    print(sack(expr))
    
    expr = expand((1-t)*(1-t)*ax+2*(1-t)*t*bx+t*t*cx)
    expr = collect(expr,t)
    print("Bezier Polynomial")
    print(sack(expr))

#algebras()


#https://en.wikipedia.org/wiki/Golden-section_search    
"""Python program for golden section search.  This implementation
   reuses function evaluations, saving 1/2 of the evaluations per
   iteration, and returns a bounding interval."""

invphi = (math.sqrt(5) - 1) / 2  # 1 / phi
invphi2 = (3 - math.sqrt(5)) / 2  # 1 / phi^2

def gss(f, a, b, tol=1e-5):
    """Golden section search.

    Given a function f with a single local minimum in
    the interval [a,b], gss returns a subset interval
    [c,d] that contains the minimum with d-c <= tol.

    Example:
    >>> f = lambda x: (x-2)**2
    >>> a = 1
    >>> b = 5
    >>> tol = 1e-5
    >>> (c,d) = gss(f, a, b, tol)
    >>> print(c, d)
    1.9999959837979107 2.0000050911830893
    """

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return (a, b)

    # Required steps to achieve tolerance
    n = int(math.ceil(math.log(tol / h) / math.log(invphi)))
    print(n)
    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    for k in range(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            yc = f(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            yd = f(d)

    if yc < yd:
        t = (a+d)/2
        #return (a, d)
    else:
        t = (c+b)/2
        #return (c, b)
    return t

def add(X,Y):
    return (X[0]+Y[0],X[1]+Y[1])

def sub(X,Y):
    return (X[0]-Y[0],X[1]-Y[1])

def dot(X,Y):
    return X[0]*Y[0]+X[1]*Y[1]

#https://www.wyzant.com/resources/answers/607711/determine-which-side-of-a-line-a-point-lies
def det(P,X,Y):
    #d=(x-x1)(y2-y1) - (y-y1)(x2-x1)
    return (P[0]-X[0])*(Y[1]-X[1])-(P[1]-X[1])*(Y[0]-X[0])
    
def mul(k,X):
    return (X[0]*k,X[1]*k)

def dist(t,Alpha,Beta,Gamma):
    D = add(Gamma,mul(t,add(Beta,mul(t,Alpha))))
    return dot(D,D)

def bez(t,A,B,C):
    return add(add(mul((1-t)*(1-t),A),mul(2*(1-t)*t,B)),mul(t*t,C))

def minima(A,B,C,D):
    #first find critical point
    #this is garaunteed to be a unimodal function
    Alpha = add(sub(A,mul(2,B)),C)
    Beta  = mul(2,sub(B,A))
    f = lambda t: dist(t,Alpha,Beta,sub(A,B))
    tol = 0.01
    tc = gss(f,0,1,tol)
    P = bez(tc,A,B,C)
    print(P)
    #now find which side of the critial line D is on
    #then minimize distance to D on the sub-interval
    f = lambda t: dist(t,Alpha,Beta,sub(A,D))
    if det(A,B,P)*det(D,B,P) >= 0:
        print("left")
        tm = gss(f,0,tc,tol)
    else:
        print("right")
        tm = gss(f,tc,1,tol)
    
    return bez(tm,A,B,C)
    

    

A = (1.2,0.2)
B = (1.2,1.0)
C = (1.6,0.2)
D = (1.3,0.4)
P = minima(A,B,C,D)

print(P)

