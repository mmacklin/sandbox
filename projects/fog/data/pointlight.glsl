#include "common.glsl"

uniform mat4 g_lightToWorld;
uniform mat4 g_worldToLight;

uniform vec4 g_lightPos;
uniform vec4 g_lightCol;
uniform float g_scatteringCoefficient;

uniform sampler3D g_noiseTexture;

void main()
{
	vec3 surfacePos = gl_TexCoord[0].xyz;
	vec3 surfaceNormal = gl_TexCoord[1].xyz;
	vec3 cameraPos = gl_ModelViewMatrixInverse[3].xyz;
	vec3 lightPos = g_lightToWorld[3].xyz;

	//vec3 noise = texture3D(g_noiseTexture, 0.1*(surfacePos + (surfacePos-cameraPos)*0.5)).xyz;

	vec3 dir = surfacePos - cameraPos;// + noise;
	float l = length(dir);
	dir /= l;

	// calculate in-scattering contribution
	vec3 scatter = g_lightCol.xyz * vec3(0.2, 0.5, 0.8) * InScatter(cameraPos, dir, lightPos, l) * g_scatteringCoefficient;
	
	// calculate contribution from diffuse reflection
	float diffuse;
	float specular;

	float specularExponent = 15.0;
	float specularIntensity = 0.02;

	Li(surfacePos, surfaceNormal, lightPos, cameraPos, specularExponent, diffuse, specular);

	vec3 r = LinearToSrgb(g_lightCol.xyz * (diffuse + specular * specularIntensity) + scatter);

	gl_FragColor = vec4(r, 1.0);
	//gl_FragColor = vec4(texture3D( g_noiseTexture, surfacePos * 0.1).xyz, 1.0);
}