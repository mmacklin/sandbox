#version 120

uniform sampler3D g_densityTexture;
uniform vec3 g_lightPos;
uniform vec3 g_lightIntensity;
uniform vec3 g_absorption;
uniform vec3 g_scatter;

void main()
{
	const int numLightSamples = 256;
	const float lscale = 8.0 / float(numLightSamples);
	
	vec3 pos = gl_TexCoord[0].xyz;
	vec3 eyePos = gl_ModelViewMatrixInverse[3].xyz;
	
	// attenuation
	float falloff = 1.0 / dot(pos-g_lightPos, pos-g_lightPos);
	
	// convert to texture space
	pos = pos*0.5 + vec3(0.5, 0.5, 0.5);
	vec3 lightPos = g_lightPos*0.5 +vec3(0.5, 0.5, 0.5);
	
	// convert to texture space (cube dimensions are -1, 1 in world space)	
	vec3 transmittance = vec3(1.0, 1.0, 1.0);
	vec3 Lo = vec3(0.0, 0.0, 0.0);

	vec3 lightDir = normalize(lightPos - pos)*lscale*0.5;

	// sample light
	vec3 ltransmittance = vec3(1.0, 1.0, 1.0);
	vec3 lpos = pos + lightDir;
	
	for (int s=0; s < numLightSamples; ++s)
	{
		ltransmittance *= 1.0-g_scatter*lscale*texture3D(g_densityTexture, lpos).x;			

		if (ltransmittance.x <= 0.01)
			break;
		
		lpos += lightDir;
	}
	
	vec3 Li = ltransmittance;
	Lo += Li*vec3(0.52f, 0.46f, 0.4f) + vec3(0.05, 0.05, 0.05);


	gl_FragColor.xyz = Lo*falloff*g_lightIntensity;
	gl_FragColor.w = 0.0f;//max(0.0, 1.0-dot(transmittance, vec3(0.33, 0.33, 0.33)));
}