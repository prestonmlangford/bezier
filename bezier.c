typedef vec2;
#define gpoly(t) fma(fma(fma(h,t,i),t,j),t,k)
#define zpoly(t) fma(fma(Alpha,t,Beta),t,Gamma)

float halley(float t, float h,float i,float j,float k){
    float gpp,gp,g = 1;
    float i2 = 2*i;
    float h3 = 3*h;
    float h6 = 2*h3;
    
    #define halley_step \
    g = fma(fma(fma(h,t,i),t,j),t,k);\
    gp = fma(fma(h3,t,i2),t,j);\
    gpp = fma(h6,t,i2);\
    t = t - g*gp/(gp*gp-0.5*g*gpp);\
    
    halley_step
    halley_step
    halley_step
    
    //while(abs(g) > 0.01){//Maybe try loop unrolling and always assume it finishes in 3 iterations
    //    halley_step
    //}
    return t;
}

float min_distance(vec2 A,vec2 B,vec2 C,vec2 D){
    vec2 Alpha = A-2*B+C;
    vec2 Beta  = 2*(B-A);
    vec2 Gamma = A-D;
    float h = 2*dot(Alpha,Alpha);
    float i = 3*dot(Alpha,Beta);
    float j = 2*dot(Alpha,Gamma)+dot(Beta,Beta);
    float k = dot(Beta,Gamma);
    float h3 = 3*h;
    float v = -i/h3;
    float phi = v*v-j/h3;
    float g0 = k;
    float g1 = h+i+j+k;
    float gv = gpoly(v);
    
    if (phi > 0){
        
        float sqrtphi = sqrt(phi);
        float l = v - sqrtphi;
        float r = v + sqrtphi;
        float gl = gpoly(l);
        float gr = gpoly(r);
            
        
        int lm = ((l > 0) && (gl > 0) && (g0 < 0) && (l < 1)) || ((l > 1) && (g0 < 0) && (g1 > 0));
        int rm = ((r < 1) && (gr < 0) && (g1 > 0) && (r > 0)) || ((r < 0) && (g0 < 0) && (g1 > 0));
        int ls = gv > 0;
        int rs = 1-ls;
        float t;
        
        if(ls && lm){
            float t = halley(v-2*sqrtphi,h,i,j,k);
            return length(zpoly(t));
        }
        else if(rs && rm){
            float t = halley(v+2*sqrtphi,h,i,j,k);
            return length(zpoly(t));
        }
        else if(ls && rm){
            t =  halley(v+2*sqrtphi,h,i,j,k);
            return min(length(zpoly(t)),length(Gamma));
        }
        else q1 `if rs && lm:
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
            t =  halley(1,g,gp,gpp,tol);

    return length(fma(fma(Alpha,t,Beta),t,Gamma)));
}