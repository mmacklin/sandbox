#include "shaders.h"

#include "core/shader.h"

// vertex shader
const char *vertexPointShader = STRINGIFY(

uniform float pointRadius;  // point size in world space
uniform float pointScale;   // scale to calculate size in pixels

uniform mat4 lightTransform; 
uniform vec3 lightDir;
uniform vec3 lightDirView;

void main()
{
    // calculate window-space point size
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
	gl_PointSize = pointRadius * (pointScale / gl_Position.w);    

	gl_TexCoord[0] = gl_MultiTexCoord0;    
	gl_TexCoord[1] = lightTransform*vec4(gl_Vertex.xyz-lightDir*pointScale*2.0, 1.0);
	gl_TexCoord[2] = gl_ModelViewMatrix*vec4(lightDir, 0.0);
	gl_TexCoord[3] = gl_Color;
}
);

// pixel shader for rendering points as shaded spheres
const char *fragmentPointShader = STRINGIFY(

uniform vec3 lightDir;
uniform vec3 lightPos;

uniform sampler2DShadow shadowTex;
uniform vec2 shadowTaps[12];


// sample shadow map
float shadowSample()
{
	vec3 pos = vec3(gl_TexCoord[1].xyz/gl_TexCoord[1].w);
	vec3 uvw = (pos.xyz*0.5)+vec3(0.5);

	// user clip
	if (uvw.x  < 0.0 || uvw.x > 1.0)
		return 0.0;
	if (uvw.y < 0.0 || uvw.y > 1.0)
		return 0.0;
	
	float s = 0.0;
	float radius = 0.08;

	for (int i=0; i < 8; i++)
	{
		s += shadow2D(shadowTex, vec3(uvw.xy + shadowTaps[i]*radius, uvw.z)).r;
	}

	s /= 8.0;
	return s;
}

void main()
{
    // calculate normal from texture coordinates
    vec3 normal;
    normal.xy = gl_TexCoord[0].xy*vec2(2.0, -2.0) + vec2(-1.0, 1.0);
    float mag = dot(normal.xy, normal.xy);
    if (mag > 1.0) discard;   // kill pixels outside circle
   	normal.z = sqrt(1.0-mag);

    // calculate lighting
	float shadow = shadowSample();
    vec3 diffuse = vec3(1.0, 1.0, 1.0)*max(0.0, -dot(gl_TexCoord[2].xyz, normal))*shadow;

    gl_FragColor = vec4(diffuse, 1.0);
//	gl_FragColor = vec4(N.x, N.y, N.z, 1.0);
//	gl_FragColor = vec4(dFdx(pos.x), dFdx(pos.y), 0.0, 1.0);
//	gl_FragColor.xyz = (gl_TexCoord[1].xyz/gl_TexCoord[1].w)*0.5 + 0.5;
//	gl_FragColor.w = 1.0;
//	gl_FragColor.z = 0.0;

}


);

// vertex shader
const char *vertexShader = STRINGIFY(

uniform mat4 lightTransform; 
uniform vec3 lightDir;

void main()
{
    // calculate window-space point size
	gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

	gl_TexCoord[0].xyz = gl_Normal;    
	gl_TexCoord[1] = lightTransform*vec4(gl_Vertex.xyz, 1.0);
	gl_TexCoord[2] = gl_ModelViewMatrix*vec4(lightDir, 0.0);
}
);

// pixel shader for rendering points as shaded spheres
const char *fragmentShader = STRINGIFY(

uniform vec3 lightDir;
uniform vec3 lightPos;

uniform sampler2DShadow shadowTex;
uniform vec2 shadowTaps[12];


// sample shadow map
float shadowSample()
{
	vec3 pos = vec3(gl_TexCoord[1].xyz/gl_TexCoord[1].w);
	vec3 uvw = (pos.xyz*0.5)+vec3(0.5);

	// user clip
	if (uvw.x  < 0.0 || uvw.x > 1.0)
		return 1.0;
	if (uvw.y < 0.0 || uvw.y > 1.0)
		return 1.0;
	
	float s = 0.0;
	float radius = 0.002;

	for (int i=0; i < 8; i++)
	{
		s += shadow2D(shadowTex, vec3(uvw.xy + shadowTaps[i]*radius, uvw.z)).r;
	}

	s /= 8.0;
	return s;
}

void main()
{
    // calculate lighting
	float shadow = shadowSample();
    float diffuse = max(0.0, -dot(lightDir, gl_TexCoord[0].xyz))*shadow;

    gl_FragColor = vec4(diffuse, diffuse, diffuse, 1.0);
}


);

void ShadowApply(GLint sprogram, Vec3 lightPos, Vec3 lightTarget, Matrix44 lightTransform, GLuint shadowTex)
{
	GLint uLightTransform = glGetUniformLocation(sprogram, "lightTransform");
	glUniformMatrix4fv(uLightTransform, 1, false, lightTransform);

	GLint uLightPos = glGetUniformLocation(sprogram, "lightPos");
	glUniform3fv(uLightPos, 1, lightPos);
	
	GLint uLightDir = glGetUniformLocation(sprogram, "lightDir");
	glUniform3fv(uLightDir, 1, Normalize(lightTarget-lightPos));

	const Vec2 taps[] = 
	{ 
		Vec2(-0.326212,-0.40581),Vec2(-0.840144,-0.07358),
		Vec2(-0.695914,0.457137),Vec2(-0.203345,0.620716),
		Vec2(0.96234,-0.194983),Vec2(0.473434,-0.480026),
		Vec2(0.519456,0.767022),Vec2(0.185461,-0.893124),
		Vec2(0.507431,0.064425),Vec2(0.89642,0.412458),
		Vec2(-0.32194,-0.932615),Vec2(-0.791559,-0.59771) 
	};
	
	GLint uShadowTaps = glGetUniformLocation(sprogram, "shadowTaps");
	glUniform2fv(uShadowTaps, 12, &taps[0].x);
	
	glEnable(GL_TEXTURE_2D);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, shadowTex);

}


void DrawPoints(float* positions, int n, float radius, float screenWidth, float screenAspect, float fov, Vec3 lightPos, Vec3 lightTarget, Matrix44 lightTransform, GLuint shadowTex)
{
	static int sprogram = -1;
	if (sprogram == -1)
		sprogram = CompileProgram(vertexPointShader, fragmentPointShader);

	if (sprogram)
	{
		glEnable(GL_POINT_SPRITE);
		glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glDepthMask(GL_TRUE);
		glEnable(GL_DEPTH_TEST);		

		glUseProgram(sprogram);
		glUniform1f( glGetUniformLocation(sprogram, "pointScale"), radius*2.0f);
		glUniform1f( glGetUniformLocation(sprogram, "pointRadius"), screenWidth / (2.0f*screenAspect*tanf(fov*0.5f)));

		// set shadow parameters
		ShadowApply(sprogram, lightPos, lightTarget, lightTransform, shadowTex);	

		// diffuse 
		glColor3f(0.5, 1, 0.5);

		glEnableClientState(GL_VERTEX_ARRAY);			
		glVertexPointer(3, GL_FLOAT, sizeof(float)*3, positions);

		glDrawArrays(GL_POINTS, 0, n);

		glUseProgram(0);
		glDisableClientState(GL_VERTEX_ARRAY);	
		glDisable(GL_POINT_SPRITE_ARB);
	}
}

void DrawPlane(const Vec4& p);

void DrawPlanes(Vec4* planes, int n, Vec3 lightPos, Vec3 lightTarget, Matrix44 lightTransform, GLuint shadowTex)
{
	static int sprogram = -1;
	if (sprogram == -1)
		sprogram = CompileProgram(vertexShader, fragmentShader);

	if (sprogram)
	{
		glDepthMask(GL_TRUE);
		glEnable(GL_DEPTH_TEST);		

		glUseProgram(sprogram);

		// set shadow parameters
		ShadowApply(sprogram, lightPos, lightTarget, lightTransform, shadowTex);	

		// diffuse 
		glColor3f(1, 1, 1);

		for (int i=0; i < n; ++i)
			DrawPlane(planes[i]);
	}
}
	
void ShadowCreate(GLuint& texture, GLuint& frameBuffer)
{
	glVerify(glGenFramebuffers(1, &frameBuffer));
	glVerify(glGenTextures(1, &texture));
	glVerify(glBindTexture(GL_TEXTURE_2D, texture));

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
	 
	// This is to allow usage of shadow2DProj function in the shader 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_R_TO_TEXTURE); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL); 
	glTexParameteri(GL_TEXTURE_2D, GL_DEPTH_TEXTURE_MODE, GL_INTENSITY); 

	glVerify(glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 1024, 1024, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL));

}

void ShadowBegin(GLuint texture, GLuint frameBuffer)
{
	glVerify(glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer));
	glVerify(glDrawBuffer(GL_NONE));
	glVerify(glReadBuffer(GL_NONE));
	glVerify(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, texture, 0));

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, 1024, 1024);
}

void ShadowEnd()
{
	glVerify(glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, 0));
	glVerify(glBindFramebuffer(GL_FRAMEBUFFER, 0));
}
