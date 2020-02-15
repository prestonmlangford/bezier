from sympy import *
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

def diff(t,alpha,beta,gamma): return gamma + t*(beta + alpha*t)


