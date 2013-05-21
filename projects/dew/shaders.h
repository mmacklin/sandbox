#pragma once

// spherical harmonic coefficients for the 'beach' light probe
Vec3 gShDiffuse[] = 
{
	Vec3(1.51985, 1.53785, 1.56834),
	Vec3(-0.814902, -0.948101, -1.13014),
	Vec3(-0.443242, -0.441047, -0.421306),
	Vec3(1.16161, 1.07284, 0.881858),
	Vec3(-0.36858, -0.37136, -0.332637),
	Vec3(0.178697, 0.200577, 0.219209),
	Vec3(-0.0204381, -0.0136351, -0.00920174),
	Vec3(-0.349078, -0.292836, -0.214752),
	Vec3(0.399496, 0.334641, 0.219389),
};

// vertex shader
const char* vertexShader = STRINGIFY
(
	uniform mat4 lightTransform; 
	uniform mat4 worldTransform;

	void main()
	{
		vec3 n = gl_Normal;//normalize(gl_Normal);

		gl_Position = gl_ModelViewProjectionMatrix*vec4(gl_Vertex.xyz, 1.0);		
		gl_TexCoord[0] = worldTransform*vec4(n, 0.0);
		gl_TexCoord[1] = vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[2] = lightTransform*vec4(gl_Vertex.xyz+n, 1.0);
		gl_TexCoord[3].xyz = gl_Color.xyz;
	}
);

// pixel shader
const char* fragmentShaderMain = STRINGIFY
(
	uniform vec3 shDiffuse[9];
	uniform vec3 color;

	uniform vec3 lightDir;
	uniform vec3 lightPos;

	uniform sampler2DShadow shadowTex;
	uniform vec2 shadowTaps[12];

	// evaluate spherical harmonic function
	vec3 shEval(vec3 dir, vec3 sh[9])
	{
		// evaluate irradiance
		vec3 e = sh[0];

		e += -dir.y*sh[1];
		e +=  dir.z*sh[2];
		e += -dir.x*sh[3];

		e +=  dir.x*dir.y*sh[4];
		e += -dir.y*dir.z*sh[5];
		e += -dir.x*dir.z*sh[7];

		e += (3.0*dir.z*dir.z-1.0)*sh[6];
		e += (dir.x*dir.x - dir.y*dir.y)*sh[8];

		return max(e, vec3(0.0));
	}

	// sample shadow map
	float shadowSample()
	{
		vec3 pos = vec3(gl_TexCoord[2].xyz/gl_TexCoord[2].w);
		vec3 uvw = (pos.xyz*0.5)+vec3(0.5);

		// user clip
		if (uvw.x  < 0.0 || uvw.x > 1.0)
			return 0.0;
		if (uvw.y < 0.0 || uvw.y > 1.0)
			return 0.0;
		
		float s = 0.0;
		float radius = 0.003;

		for (int i=0; i < 12; i++)
		{
			s += shadow2D(shadowTex, vec3(uvw.xy + shadowTaps[i]*radius, uvw.z)).r;
		}

		s /= 12.0;
		return s;
	}

	void main()
	{
		vec3 n = gl_TexCoord[0].xyz;
		vec3 shadePos = gl_TexCoord[1].xyz;
		vec3 eyePos = gl_ModelViewMatrixInverse[3].xyz;
		vec3 eyeDir = normalize(eyePos-shadePos);
		vec3 vcolor = gl_TexCoord[3].xyz;
	
		vec3 lightCol = vec3(1.0, 1.0, 1.0)*0.8; 

		// SH-ambient	
		float ambientExposure = 0.04;
		vec3 ambient = shEval(n, shDiffuse)*ambientExposure;

		// wrapped spot light 
		float w = 0.1;
		float s = shadowSample();
		vec3 direct = clamp((dot(n, -lightDir) + w) / (1.0 + w), 0.0, 1.0)*lightCol*smoothstep(0.9, 1.0, dot(lightDir, normalize(shadePos-lightPos))); 

		vec3 l = (ambient + s*direct)*vcolor;//*(n*vec3(0.5) + vec3(0.5));//color;

		// convert from linear light space to SRGB
		gl_FragColor = vec4(pow(l, vec3(0.5)), 1.0);
	}
);

// pixel shader
const char* fragmentShaderDebug = STRINGIFY
(
	void main()
	{
		vec3 n = gl_TexCoord[0].xyz;
		gl_FragColor = vec4(n*0.5 + vec3(0.5), 1.0);
	}
);

