//bezier_sd and inv are from https://www.shadertoy.com/view/4sySDK
//area formula from https://www.geeksforgeeks.org/program-check-three-points-collinear/
//mainImage from https://www.shadertoy.com/view/XsX3zf

float area(vec2 a, vec2 b, vec2 c) 
{ 
    return abs(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y));
} 

float linear_sd(vec2 uv, vec2 p0, vec2 p1){
    vec2 v = normalize(p0-p1);
    vec2 n = vec2(v.y,-v.x);
    return dot(n,uv-p0);
}

mat2 inv(mat2 m){
	return mat2(m[1][1],-m[0][1],-m[1][0],m[0][0])/(m[0][0]*m[1][1]-m[1][0]*m[0][1]);
}


float bezier_sd(vec2 uv, vec2 p0, vec2 p1, vec2 p2){

	const mat2 trf1 = mat2(-1, 2, 1, 2);
	mat2 trf2 = inv(mat2(p0-p1, p2-p1));
	mat2 trf=trf1*trf2;

	uv-=p1;
	vec2 xy=trf*uv;
	xy.y-=1.;

	vec2 gradient;
	gradient.x=2.*trf[0][0]*(trf[0][0]*uv.x+trf[1][0]*uv.y)-trf[0][1];
	gradient.y=2.*trf[1][0]*(trf[0][0]*uv.x+trf[1][0]*uv.y)-trf[1][1];

	return (xy.x*xy.x-xy.y)/length(gradient);
}



float corrected_sd(vec2 uv, vec2 p0, vec2 p1, vec2 p2){
    float b = bezier_sd(uv,p0,p1,p2);
    float l = linear_sd(uv,p0,p2);
    float t = tanh(abs(2.1*l));
    return t*b+(1.0-t)*l;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 xy = fragCoord.xy;
	
	vec2 b0 = vec2(0.25, .5) * iResolution.xy;
	// vec2 b1 = vec2(0.5, .75 + .1*sin(iTime)) * iResolution.xy;
	vec2 b1 = iMouse.xy;
	vec2 b2 = vec2(.75, .5) * iResolution.xy;
	vec2 mid = .5*(b0+b2) + vec2(0.0,0.01);
	
	float d = abs(bezier_sd(xy, b0,b1,b2));
	float thickness = 1.0;
	
	float a;
	
	if(d < thickness) {
	  a = 1.0;
	} else {
	  // Anti-alias the edge.
	  a = 1.0 - smoothstep(d, thickness, thickness+1.0);
	}
	
	fragColor = vec4(a,0.0,0.0,0.0);
	
	//fragColor = vec4(mod(d/50.0, 1.0),a,a, 1.0);
}