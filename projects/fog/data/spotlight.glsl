#include "common.glsl"

uniform sampler2D g_texDiffuse;
uniform sampler2D g_texNormal;
uniform sampler2DShadow g_texDepth; 
uniform mat4 g_lightToWorld;
uniform mat4 g_worldToLight;

uniform vec4 g_lightPos;
uniform vec4 g_lightCol;
uniform float g_scatteringCoefficient;

uniform sampler3D g_noiseTexture;

void IntersectCone(vec3 rayOrigin, vec3 rayDir, mat4 invConeTransform, float aperture, float height, out float minT, out float maxT)
{
	vec4 localOrigin = invConeTransform * vec4(rayOrigin, 1.0);
	vec4 localDir = invConeTransform * vec4(rayDir, 0.0);

	// could perform this on the cpu
	float tanTheta = tan(aperture);
	tanTheta *= tanTheta;

	float a = localDir.x*localDir.x + localDir.z*localDir.z - localDir.y*localDir.y*tanTheta;
	float b = 2.0*(localOrigin.x*localDir.x + localOrigin.z*localDir.z - localOrigin.y*localDir.y*tanTheta);
	float c = localOrigin.x*localOrigin.x + localOrigin.z*localOrigin.z - localOrigin.y*localOrigin.y*tanTheta;

	SolveQuadratic(a, b, c, minT, maxT);

	float y1 = localOrigin.y + localDir.y*minT;
	float y2 = localOrigin.y + localDir.y*maxT;

	// should be possible to simplify these branches if the compiler isn't already doing it
	
	if (y1 > 0.0 && y2 > 0.0)
	{
		// both intersections are in the reflected cone so return degenerate value
		minT = 0.0;
		maxT = -1.0;
	}
	else if (y1 > 0.0 && y2 < 0.0)
	{
		// closest t on the wrong side, furthest on the right side => ray enters volume but doesn't leave it (so set maxT arbitrarily large)
		minT = maxT;
		maxT = 10000.0;
	}
	else if (y1 < 0.0 && y2 > 0.0)
	{
		// closest t on the right side, largest on the wrong side => ray starts in volume and exits once
		maxT = minT;
		minT = 0.0;		
	}	
}

float SpotFalloff(vec3 surfacePos, float aperture, vec3 lightDir, vec3 lightPos)
{
	vec3 l = lightPos-surfacePos;
	float d = length(l);
	l /= d;

	float a = saturate(1.0 - acos(saturate(dot(l, lightDir))) / aperture);

	return a;
}

void main()
{
	vec3 surfacePos = gl_TexCoord[0].xyz;
	vec3 surfaceNormal = gl_TexCoord[1].xyz;
	vec3 cameraPos = gl_ModelViewMatrixInverse[3].xyz;
	vec3 lightPos = g_lightToWorld[3].xyz;
	vec3 lightDir = g_lightToWorld[1].xyz;

	vec3 dir = surfacePos - cameraPos;
	float l = length(dir);
	dir /= l;
	
	// hard-coded light shape
	float aperture = 0.4;
	float height = 30.0;
	float minT = 0.0;
	float maxT = 0.0;
	
	IntersectCone(cameraPos, dir, g_worldToLight, aperture, height, minT, maxT);
		
	// clamp bounds to scene geometry / camera
	maxT = clamp(maxT, 0.0, l);
	minT = max(0.0, minT);

	float t = max(0.0, maxT - minT);
	
	float innerMin = (minT + maxT)*0.5 - t * 0.6 * 0.5;

	//vec3 noise = texture3D(g_noiseTexture, 0.9*(cameraPos + (surfacePos-cameraPos)*0.1)).xyz;
	//dir  = normalize(dir + noise);

	vec3 scatter = g_lightCol.xyz * vec3(0.2, 0.5, 0.8) * InScatter(cameraPos + dir*minT, dir, lightPos, t)  * g_scatteringCoefficient;
	
	// calculate contribution from diffuse reflection
	float diffuse;
	float specular;
	float falloff = SpotFalloff(surfacePos, aperture, lightDir, lightPos);
	
	float specularExponent = 15.0;
	float specularIntensity = 0.02;
	
	Li(surfacePos, surfaceNormal, lightPos, cameraPos, specularExponent, diffuse, specular);
	
	vec3 r = LinearToSrgb(g_lightCol.xyz * (falloff * (diffuse + specular * specularIntensity)) + scatter);

	gl_FragColor = vec4(r, 1.0);

}