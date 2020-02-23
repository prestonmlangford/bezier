//bezier_sd and inv are from https://www.shadertoy.com/view/4sySDK
//area formula from https://www.geeksforgeeks.org/program-check-three-points-collinear/
//mainImage from https://www.shadertoy.com/view/XsX3zf

float area(vec2 a, vec2 b, vec2 c) 
{ 
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y);
} 

float approx_distance(vec2 p, vec2 a, vec2 b, vec2 c) {
    float alpha = area(p,b,a);
    float beta  = area(p,c,b);
    float gamma = area(p,a,c);
    
    return 4.0*alpha*beta-gamma*gamma;
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