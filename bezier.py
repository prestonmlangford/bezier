import sympy
from sympy import expand
from sympy import collect
from sympy import collect
from sympy import symbols
from sympy import init_printing
from sympy import Symbol
from sympy import diff
import math
import matplotlib.pyplot as plt
import numpy as np

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

    x = Symbol('x')
    

    # Factorial function
    def factorial(n):
        if n <= 0:
            return 1
        else:
            return n*factorial(n-1)

    # Taylor approximation at x0 of the function 'function'
    def taylor(function,x0,n):
        i = 0
        p = 0
        while i <= n:
            p = p + (function.diff(x,i).subs(x,x0))/(factorial(i))*(x-x0)**i
            i += 1
        return p
    for i in [0,1,2,3]:
        ts = collect(expand(taylor(x**0.5,0.5,i)),x)
        print("Talyor Expansion Order {}".format(i))
        print(sack(ts))

#algebras()

ops = 0


def add(X,Y):
    global ops
    ops +=1
    return (X[0]+Y[0],X[1]+Y[1])

def sub(X,Y):
    global ops
    ops +=1
    return (X[0]-Y[0],X[1]-Y[1])

def dot(X,Y):
    global ops
    ops +=2
    return X[0]*Y[0]+X[1]*Y[1]

def det(A,B): 
    global ops
    ops +=3
    return A[0]*B[1]-B[0]*A[1]
    
def norm2(X):
    return dot(X,X)
    
def mul(k,X):
    global ops
    ops +=1
    return (X[0]*k,X[1]*k)

def clamp(x,lower,upper):
    global ops
    ops += 2
    return min(max(lower,x),upper)

def fma(X,y,Z):
    global ops
    ops +=1
    return add(mul(y,X),Z)

def lerp(A,B,t):
    global ops
    ops +=1
    return add(mul(1-t,A),mul(t,B))

def fusqrt(x):
    global ops
    ops += 3 # using fused multiply add
    # returns upper bound for sqrt
    return 0.220970869120796 + x*(1.32582521472478 + x*(-0.883883476483184 + x*0.353553390593274))

# http://www.pouet.net/topic.php?which=9119&page=1
# this seems like a good exact solution which is faster than 
# iteratively solving for the roots

# From microsoft research paper for comparison purposes
# Find vector v given pixel 𝑝=(0,0) and Bézier points b0,b1,b2
def get_distance_vector(p,b0,b1,b2):
    global ops
    a,b,d = det(b0,b2), 2*det(b1,b0), 2*det(b2,b1)          # float a=det(b0,b2), b=2*det(b1,b0), d=2*det(b2,b1);
    f=b*d-a*a                                               # float f=b*d-a*a;
    d21,d10,d20 = sub(b2,b1), sub(b1,b0), sub(b2,b0)        # float d21=b2-b1, d10=b1-b0, d20=b2-b0;
    gf=mul(2,add(add(mul(b,d21),mul(d,d10)),mul(a,d20)))    # float2 gf=2*(b*d21+d*d10+a*d20);
    gf=(gf[1],-gf[0])                                       # gf=float2(gf.y,-gf.x);
    pp=mul(-f/dot(gf,gf),gf)                                # float2 pp=-f*gf/dot(gf,gf);
    d0p=sub(b0,pp)                                          # float2 d0p=b0-pp;
    ap,bp=det(d0p,d20), 2*det(d10,d0p)                      # float ap=det(d0p,d20), bp=2*det(d10,d0p);
    # (note that 2*ap+bp+dp=2*a+b+d=4*area(b0,b1,b2))       # 
    t=clamp((ap+bp)/(2*a+b+d), 0,1)                         # float t=clamp((ap+bp)/(2*a+b+d), 0,1);
    v=lerp(lerp(b0,b1,t),lerp(b1,b2,t),t)                   # return lerp(lerp(b0,b1,t),lerp(b1,b2,t),t);
    #return math.sqrt(dot(v,v))
    ops += 14
    return math.sqrt(dot(v,v))



def halley(t0,g,gp,gpp,tol=0.01):
    global ops
    t = t0
    gt = tol + 1
    iterations = 0
    # criteria is modified to stop when g is close to zero.
    # Since we only need the distance to the curve at the closest point, 
    # it does not matter where precisely the closest point is on the curve
    # This makes a difference in the number of iterations for areas of the plane that are roughly equidistant to a region of the parabola.  
    # Specifically, the region of greatest curvature on the parabola
    while abs(gt) > tol:
        gt = g(t)
        gpt = gp(t)
        gppt = gpp(t)
        t = t - 2*gt*gpt/(2*gpt*gpt-gt*gppt)
        ops += 15
        iterations += 1
    print("Halley Iterations: {}".format(iterations))
    return t

def min_distance(A,B,C,D,tol = 0.01):
    global ops
    # 2nd order polynomial coefficients for the bezier curve Z
    Alpha = add(sub(A,mul(2,B)),C)
    Beta  = mul(2,sub(B,A))
    Gamma = sub(A,D)

    # 3rd order polynomial coefficients for derivative of distance to the bezier curve "g"
    h = 2*dot(Alpha,Alpha)
    i = 3*dot(Alpha,Beta)
    j = 2*dot(Alpha,Gamma)+dot(Beta,Beta)
    k = dot(Beta,Gamma)
    h3 = 3*h
    h6 = 6*h
    i2 = 2*i
    #these can be implemented with the fma instruction to make them more efficient
    Z   = lambda t: add(Gamma,mul(t,add(Beta,mul(t,Alpha))))
    g   = lambda t: k+t*(j+t*(i+t*h))
    gp  = lambda t: j+t*(i2+t*h3)
    gpp = lambda t: i2+t*h6
    
    # v is the vertex of the parabola Z
    v = -i/h3

    # symmetry about the vertex for zeros of the 2nd derivative of distance
    phi = v*v-j/h3
    
    g0 = k
    g1 = h+i+j+k
    gv = g(v)
    
    ops += 15

    # when phi > 0, there are local maxima and minima of g.  These could produce zeros.
    if phi > 0:# could also use g'(v) < 0, but why waste the clock cycles when I need phi anyway
        sqrtphi = math.sqrt(phi)
        
        l = v - sqrtphi
        r = v + sqrtphi
        gl = g(l)
        gr = g(r)
        
        lm = ((l > 0) and (gl > 0) and (g0 < 0) and (l < 1)) or ((l > 1) and (g0 < 0) and (g1 > 0))
        rm = ((r < 1) and (gr < 0) and (g1 > 0) and (r > 0)) or ((r < 0) and (g0 < 0) and (g1 > 0))
        ls = gv > 0
        rs = not ls
        
        ops += 18

        if ls and lm:
            t =  halley(v-2*sqrtphi,g,gp,gpp,tol)
        elif rs and rm:
            t =  halley(v+2*sqrtphi,g,gp,gpp,tol)
        elif ls and rm:
            m = halley(v+2*sqrtphi,g,gp,gpp,tol)
            Zm = Z(m)
            Z0 = sub(A,D)
            t =  0 if dot(Z0,Z0) < dot(Zm,Zm) else m
        elif rs and lm:
            m = halley(v-2*sqrtphi,g,gp,gpp,tol)
            Zm = Z(m)
            Z1 = sub(C,D)
            t =  1 if dot(Z1,Z1) < dot(Zm,Zm) else m
        else:
            Z0 = sub(A,D)
            Z1 = sub(C,D)
            t =  0 if dot(Z0,Z0) < dot(Z1,Z1) else 1
    else:
        if g0 > 0:
            t =  0
        elif g1 < 0:
            t =  1
        elif gv > 0:
            t =  halley(0,g,gp,gpp,tol)
        else:
            t =  halley(1,g,gp,gpp,tol)

    P = Z(t)
    return math.sqrt(dot(P,P))


def bisection_method(A,B,C,D):
    global ops
    # 2nd order polynomial coefficients for the bezier curve Z
    Alpha = add(sub(A,mul(2,B)),C)
    Beta  = mul(2,sub(B,A))
    Gamma = sub(A,D)

    # 3rd order polynomial coefficients for derivative of distance to the bezier curve "g"
    h = 2*dot(Alpha,Alpha)
    i = 3*dot(Alpha,Beta)
    j = 2*dot(Alpha,Gamma)+dot(Beta,Beta)
    k = dot(Beta,Gamma)
    h3 = 3*h
    h6 = 6*h
    i2 = 2*i
    #these can be implemented with the fma instruction to make them more efficient
    Z   = lambda t: add(Gamma,mul(t,add(Beta,mul(t,Alpha))))
    g   = lambda t: k+t*(j+t*(i+t*h))
    gp  = lambda t: j+t*(i2+t*h3)
    gpp = lambda t: i2+t*h6
    
    # v is the vertex of the parabola Z
    v = -i/h3

    # symmetry about the vertex for zeros of the 2nd derivative of distance
    phi = v*v-j/h3
    
    g0 = k
    g1 = h+i+j+k
    gv = g(v)
    
    ops += 15

    # when phi > 0, there are local maxima and minima of g.  These could produce zeros.
    if phi > 0:# could also use g'(v) < 0, but why waste the clock cycles when I need phi anyway
        sqrtphi = math.sqrt(phi)
        
        l = v - sqrtphi
        r = v + sqrtphi
        gl = g(l)
        gr = g(r)
        
        lm = ((l > 0) and (gl > 0) and (g0 < 0) and (l < 1)) or ((l > 1) and (g0 < 0) and (g1 > 0))
        rm = ((r < 1) and (gr < 0) and (g1 > 0) and (r > 0)) or ((r < 0) and (g0 < 0) and (g1 > 0))
        ls = gv > 0
        rs = not ls
        
        ops += 18

        if ls and lm:
            t =  halley(v-2*sqrtphi,g,gp,gpp,tol)
        elif rs and rm:
            t =  halley(v+2*sqrtphi,g,gp,gpp,tol)
        elif ls and rm:
            m = halley(v+2*sqrtphi,g,gp,gpp,tol)
            Zm = Z(m)
            Z0 = sub(A,D)
            t =  0 if dot(Z0,Z0) < dot(Zm,Zm) else m
        elif rs and lm:
            m = halley(v-2*sqrtphi,g,gp,gpp,tol)
            Zm = Z(m)
            Z1 = sub(C,D)
            t =  1 if dot(Z1,Z1) < dot(Zm,Zm) else m
        else:
            Z0 = sub(A,D)
            Z1 = sub(C,D)
            t =  0 if dot(Z0,Z0) < dot(Z1,Z1) else 1
    else:
        if g0 > 0:
            t =  0
        elif g1 < 0:
            t =  1
        elif gv > 0:
            t =  halley(0,g,gp,gpp,tol)
        else:
            t =  halley(1,g,gp,gpp,tol)

    P = Z(t)
    return math.sqrt(dot(P,P))



#A = (1.94,0.55)
#B = (2.32,0.6)
#C = (2.42,0.35)
#D = (2.1,0.3)
#P = minima(A,B,C,D)
#P = raphson(A,B,C,D)

#print("Halley Method")
#ops = 0
#P = min_distance(A,B,C,D)
#print("Min Distance: {}".format(P))
#print("Add/Multiplies: {}".format(ops))
#
#print("Inversion Method")
#ops = 0
#P = get_distance_vector(A,B,C,D)
#print("Min Distance: {}".format(P))
#print("Add/Multiplies: {}".format(ops))


def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

plt.clf()
plt.setp(plt.gca(), autoscale_on=False)

while True:
    pts = []
    while len(pts) < 3:
        tellme('Select control points with mouse')
        pts = np.asarray(plt.ginput(3))
    
    N = 100
    A,B,C = pts
    X = np.zeros(N)
    Y = np.zeros(N)
    Z = np.zeros((N,N))
    for i in range(0,N):
        for j in range(0,N):
            x,y = i/N,j/N
            X[j] = x
            Y[i] = y
            Z[j,i] = min_distance(A,B,C,(x,y))
    
    plt.imshow(
        Z, 
        extent=[0, 1, 0, 1], 
        origin='lower',
        cmap='RdGy'
    )
    plt.axis(aspect='image')

    tellme('Click mouse to start over.')

    if plt.waitforbuttonpress():
        break

    
    plt.cla()
