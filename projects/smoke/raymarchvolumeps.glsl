#version 120

uniform sampler3D g_densityTexture;
uniform sampler3D g_temperatureTexture;
uniform sampler2D g_blackBodyTexture;
uniform sampler3D g_noiseTexture;

uniform vec3 g_lightPos;
uniform vec3 g_lightIntensity;
uniform vec3 g_absorption;
uniform vec3 g_scatter;

/*
// interpolate between b and c using derivatives defined as (c-a) and (d-b)
float CubicInterpolate(float a, float b, float c, float d, float t)
{
	float tt = t*t;
	float ttt = tt*t;
	
	float dk1 = 0.5f*(c-a);
	float dk2 = 0.5f*(d-b);

	float r = (2.0f*ttt -3.0f*tt + 1.0f)*b + (ttt - 2.0f*tt + t)*dk1 + (-2.0f*ttt + 3.0f*tt)*c + (ttt-tt)*dk2;
		
	r = clamp(r, b, c);

	return r;
	
}

float texture3DCubic(sampler3D sampler, vec3 uvw)
{
	uvw *= 128.0;
	
	float i = floor(uvw.x);
	float j = floor(uvw.y);
	float k = floor(uvw.z);
	
	// tricubic interpolation
	float tx = (uvw.x-i);
	float ty = (uvw.y-j);
	float tz = (uvw.z-k);
	
	float a0 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k-1)/128.0).x, texture3D(sampler, vec3(i, j-1, k-1)/128.0).x, texture3D(sampler, vec3(i+1, j-1, k-1)/128.0).x, texture3D(sampler, vec3(i+2, j-1, k-1)/128.0).x, tx);				
	float a1 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k-1)/128.0).x, texture3D(sampler, vec3(i, j+0, k-1)/128.0).x, texture3D(sampler, vec3(i+1, j+0, k-1)/128.0).x, texture3D(sampler, vec3(i+2, j+0, k-1)/128.0).x, tx);				
	float a2 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k-1)/128.0).x, texture3D(sampler, vec3(i, j+1, k-1)/128.0).x, texture3D(sampler, vec3(i+1, j+1, k-1)/128.0).x, texture3D(sampler, vec3(i+2, j+1, k-1)/128.0).x, tx);				
	float a3 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k-1)/128.0).x, texture3D(sampler, vec3(i, j+2, k-1)/128.0).x, texture3D(sampler, vec3(i+1, j+2, k-1)/128.0).x, texture3D(sampler, vec3(i+2, j+2, k-1)/128.0).x, tx);				

	float a4 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+0)/128.0).x, texture3D(sampler, vec3(i, j-1, k+0)/128.0).x, texture3D(sampler, vec3(i+1, j-1, k+0)/128.0).x, texture3D(sampler, vec3(i+2, j-1, k+0)/128.0).x, tx);				
	float a5 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+0)/128.0).x, texture3D(sampler, vec3(i, j+0, k+0)/128.0).x, texture3D(sampler, vec3(i+1, j+0, k+0)/128.0).x, texture3D(sampler, vec3(i+2, j+0, k+0)/128.0).x, tx);				
	float a6 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+0)/128.0).x, texture3D(sampler, vec3(i, j+1, k+0)/128.0).x, texture3D(sampler, vec3(i+1, j+1, k+0)/128.0).x, texture3D(sampler, vec3(i+2, j+1, k+0)/128.0).x, tx);				
	float a7 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+0)/128.0).x, texture3D(sampler, vec3(i, j+2, k+0)/128.0).x, texture3D(sampler, vec3(i+1, j+2, k+0)/128.0).x, texture3D(sampler, vec3(i+2, j+2, k+0)/128.0).x, tx);				

	float a8 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+1)/128.0).x, texture3D(sampler, vec3(i, j-1, k+1)/128.0).x, texture3D(sampler, vec3(i+1, j-1, k+1)/128.0).x, texture3D(sampler, vec3(i+2, j-1, k+1)/128.0).x, tx);				
	float a9 =  CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+1)/128.0).x, texture3D(sampler, vec3(i, j+0, k+1)/128.0).x, texture3D(sampler, vec3(i+1, j+0, k+1)/128.0).x, texture3D(sampler, vec3(i+2, j+0, k+1)/128.0).x, tx);				
	float a10 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+1)/128.0).x, texture3D(sampler, vec3(i, j+1, k+1)/128.0).x, texture3D(sampler, vec3(i+1, j+1, k+1)/128.0).x, texture3D(sampler, vec3(i+2, j+1, k+1)/128.0).x, tx);				
	float a11 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+1)/128.0).x, texture3D(sampler, vec3(i, j+2, k+1)/128.0).x, texture3D(sampler, vec3(i+1, j+2, k+1)/128.0).x, texture3D(sampler, vec3(i+2, j+2, k+1)/128.0).x, tx);				

	float a12 = CubicInterpolate(texture3D(sampler, vec3(i-1, j-1, k+2)/128.0).x, texture3D(sampler, vec3(i, j-1, k+2)/128.0).x, texture3D(sampler, vec3(i+1, j-1, k+2)/128.0).x, texture3D(sampler, vec3(i+2, j-1, k+2)/128.0).x, tx);				
	float a13 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+0, k+2)/128.0).x, texture3D(sampler, vec3(i, j+0, k+2)/128.0).x, texture3D(sampler, vec3(i+1, j+0, k+2)/128.0).x, texture3D(sampler, vec3(i+2, j+0, k+2)/128.0).x, tx);				
	float a14 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+1, k+2)/128.0).x, texture3D(sampler, vec3(i, j+1, k+2)/128.0).x, texture3D(sampler, vec3(i+1, j+1, k+2)/128.0).x, texture3D(sampler, vec3(i+2, j+1, k+2)/128.0).x, tx);				
	float a15 = CubicInterpolate(texture3D(sampler, vec3(i-1, j+2, k+2)/128.0).x, texture3D(sampler, vec3(i, j+2, k+2)/128.0).x, texture3D(sampler, vec3(i+1, j+2, k+2)/128.0).x, texture3D(sampler, vec3(i+2, j+2, k+2)/128.0).x, tx);				

	float b0 = CubicInterpolate(a0, a1, a2, a3, ty);
	float b1 = CubicInterpolate(a4, a5, a6, a7, ty);
	float b2 = CubicInterpolate(a8, a9, a10, a11, ty);
	float b3 = CubicInterpolate(a12, a13, a14, a15, ty);
	
	float c0 = CubicInterpolate(b0, b1, b2, b3, tz);

	return c0;
}
*/
void main()
{
	// diagonal of the cube (with sides length 2)
	const float maxDist = sqrt(3.0)*2.0;

	const int numSamples = 512;
	const float scale = maxDist/float(numSamples);
	
	const int numLightSamples = 128;
	const float lscale = maxDist / float(numLightSamples);
	
	// convert to texture space
	vec3 worldpos = gl_TexCoord[0].xyz;
	 
	vec3 pos = worldpos*0.5 + vec3(0.5);
	vec3 eyePos = gl_ModelViewMatrixInverse[3].xyz*0.5 + vec3(0.5);
	vec3 eyeDir = normalize(pos-eyePos)*scale*0.5;	

	float jitter = texture3D(g_noiseTexture, pos*20).x;
	pos += eyeDir*jitter;

	vec3 lightPos = g_lightPos*0.5 + vec3(0.5);

	vec3 T = vec3(1.0);	  // transmittance
	vec3 Lo = vec3(0.0);  // inscattered radiance

	for (int i=0; i < numSamples; ++i)
	{	
		// sample density
		float density = texture3D(g_densityTexture, pos).x;

		// skip empty space
		if (density > 0.0)
		{			

			// point light dir in texture space
			vec3 lightDir = normalize(lightPos - pos)*lscale*0.5;

			// sample light
			vec3 Tl = vec3(1.0);	// transmittance along light ray
			vec3 lpos = pos + lightDir;
			
			for (int s=0; s < numLightSamples; ++s)
			{
				Tl *= 1.0-g_scatter*lscale*texture3D(g_densityTexture, lpos).x;

				if (Tl.x <= 0.01)
					break;
				
				lpos += lightDir;
			}
		
			vec3 Li = g_lightIntensity*Tl;
			vec3 Le = vec3(0.0);//texture2D(g_blackBodyTexture, vec2(texture3D(g_temperatureTexture, pos).x, 0.0)).xyz;
			//vec3 Le = texture3D(g_temperatureTexture, pos).xyz*20.0;
			
			Lo += (Le+Li)*T*density*scale;
			
			// attenuate ray-throughput
			T *= 1.0-density*scale*g_absorption;
			
			// early out
			if (T.x <= 0.01)
				break;

		}

		pos += eyeDir;
	}

	gl_FragColor.xyz = Lo;
	gl_FragColor.w = max(0.0, 1.0-dot(T, vec3(0.33, 0.33, 0.33)));	
}