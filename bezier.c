float gpoly(float t,float h,float i,float j,float k){return k+t*(j+t*(i+t*h));}
vec2 zpoly(float t,vec2 a,vec2 b,vec2 c){return c+t*(b+t*a);}


float halley(float t, float h,float i,float j,float k){
    float gpp = 0.0;
    float gp  = 0.0;
    float g = 1.0;
    float i2 = 2.0*i;
    float h3 = 3.0*h;
    float h6 = 2.0*h3;
    
    
    //while(abs(g) > 0.01){//Maybe try loop unrolling and always assume it finishes in 3 iterations
    //    g = k+t*(j+t*(i+t*h));
    //	gp = j+t*(i2+t*h3);
    //	gpp = i2+t*h6;
    //	t = t - g*gp/(gp*gp-0.5*g*gpp);
    //}
    
    g = k+t*(j+t*(i+t*h));
    gp = j+t*(i2+t*h3);
    gpp = i2+t*h6;
    t = t - g*gp/(gp*gp-0.5*g*gpp);
    
    g = k+t*(j+t*(i+t*h));
    gp = j+t*(i2+t*h3);
    gpp = i2+t*h6;
    t = t - g*gp/(gp*gp-0.5*g*gpp);
    
    g = k+t*(j+t*(i+t*h));
    gp = j+t*(i2+t*h3);
    gpp = i2+t*h6;
    t = t - g*gp/(gp*gp-0.5*g*gpp);
    
    g = k+t*(j+t*(i+t*h));
    gp = j+t*(i2+t*h3);
    gpp = i2+t*h6;
    t = t - g*gp/(gp*gp-0.5*g*gpp);
    
    return t;
}

float min_distance(vec2 A,vec2 B,vec2 C,vec2 D){
    vec2 Alpha = A-2.0*B+C;
    vec2 Beta  = 2.0*(B-A);
    vec2 Gamma = A-D;
    float h = 2.0*dot(Alpha,Alpha);
    float i = 3.0*dot(Alpha,Beta);
    float j = 2.0*dot(Alpha,Gamma)+dot(Beta,Beta);
    float k = dot(Beta,Gamma);
    float h3 = 3.0*h;
    float v = -i/h3;
    float phi = v*v-j/h3;
    float g0 = k;
    float g1 = h+i+j+k;
    float gv = gpoly(v,h,i,j,k);
    
    if (phi > 0.0){
        
        float sqrtphi = sqrt(phi);
        float l = v - sqrtphi;
        float r = v + sqrtphi;
        float gl = gpoly(l,h,i,j,k);
        float gr = gpoly(r,h,i,j,k);
            
        
        bool lm = ((l > 0.0) && (gl > 0.0) && (g0 < 0.0) && (l < 1.0)) || ((l > 1.0) && (g0 < 0.0) && (g1 > 0.0));
        bool rm = ((r < 1.0) && (gr < 0.0) && (g1 > 0.0) && (r > 0.0)) || ((r < 0.0) && (g0 < 0.0) && (g1 > 0.0));
        bool ls = gv > 0.0;
        bool rs = !ls;
        
        
        if(ls && lm){
            float t = halley(v-2.0*sqrtphi,h,i,j,k);
            return length(zpoly(t,Alpha,Beta,Gamma));
        }
        else if(rs && rm){
            float t = halley(v+2.0*sqrtphi,h,i,j,k);
            return length(zpoly(t,Alpha,Beta,Gamma));
        }
        else if(ls && rm){
            float t =  halley(v+2.0*sqrtphi,h,i,j,k);
            return min(length(zpoly(t,Alpha,Beta,Gamma)),length(Gamma));
        }
        else if(rs && lm){
            float m = halley(v-2.0*sqrtphi,h,i,j,k);
            return min(length(zpoly(m,Alpha,Beta,Gamma)),length(C-D));
        }
        else {
            return min(length(A-D),length(C-D));
        }

    } else {
        if(g0 > 0.0)
            return length(A-D);
        else if(g1 < 0.0)
            return length(C-D);
        else if(gv > 0.0){
            float t =  halley(0.0,h,i,j,k);
            return length(zpoly(t,Alpha,Beta,Gamma));
        }
        else {
            float t =  halley(1.0,h,i,j,k);
            return length(zpoly(t,Alpha,Beta,Gamma));
        }
    }
}

// http://www.pouet.net/topic.php?which=9119&page=1
int solveCubic(float a, float b, float c, float* r)
{
	float p = b - a*a / 3;
	float q = a * (2*a*a - 9*b) / 27 + c;
	float p3 = p*p*p;
	float d = q*q + 4*p3 / 27;
	float offset = -a / 3;
	if(d >= 0) { // Single solution
		float z = sqrtf(d);
		float u = (-q + z) / 2;
		float v = (-q - z) / 2;
		u = cuberoot(u);
		v = cuberoot(v);
		r[0] = offset + u + v;
		return 1;
	}
	float u = sqrtf(-p / 3);
	float v = acos(-sqrtf( -27 / p3) * q / 2) / 3;
	float m = cos(v), n = sin(v)*1.732050808f;
	r[0] = offset + u * (m + m);
	r[1] = offset - u * (n + m);
	r[2] = offset + u * (n - m);
	return 3;
}


float SignedDistanceSquared(point2D p, QBSpline2D s)
{
    float minDis = 1e20f;

	float res[3];
		
	point2D D = s.p0 - p;

	float a = s.A;
	float b = s.B;
	float c = s.C + dot(D, s.c);
	float d = dot(D, s.d);

	int n = solveCubic(b*a, c*a, d*a, res);

	for(int j=0; j<n; j++) {
		float t = Clamp(res[j]);
		point2D d = s.p0 + (s.b + s.c*t)*t - p;
		minDis = min(minDis, dot(d,d));
	}

    return minDis;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 xy = fragCoord.xy;
	
	vec2 b0 = vec2(0.25, .5) * iResolution.xy;
	// vec2 b1 = vec2(0.5, .75 + .1*sin(iTime)) * iResolution.xy;
	vec2 b1 = iMouse.xy;
	vec2 b2 = vec2(.75, .5) * iResolution.xy;
	vec2 mid = .5*(b0+b2) + vec2(0.0,0.01);
	
	float d = min_distance(b0, b1, b2, xy);
	float thickness = 1.0;
	
	float a;
	
	if(d < thickness) {
	  a = 1.0;
	} else {
	  // Anti-alias the edge.
	  a = 1.0 - smoothstep(d, thickness, thickness+1.0);
	}
	
	//fragColor = vec4(a,1.0,1.0, 1.0);
	
	fragColor = vec4(mod(d/50.0, 1.0),a,a, 1.0);
}