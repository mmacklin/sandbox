
uniform sampler2D texPosition;
uniform sampler2D texNormal;
uniform sampler2D texIrradiance;
uniform sampler2D texAlbedo;
uniform sampler2D texIndices;

float kPi = 3.1415926535897932384626433832795;
float kInvPi = 1.0f / kPi;

float CalculateClippedArea(vec3 v, vec3 nr, vec3 ne, float a)
{
	float r2 = abs(a) / kPi;
	float r = sqrt(r2);
	float h = max(dot(v, nr), 0.0);
	float d = h / (1.0-abs(dot(nr,ne)));

	//float c = sqrt(r2-d*d);
	float theta = 2.0*acos(clamp(d/r,0.0,1.0));
	float stheta = sin(theta);
	
	float cseg = 0.5*r2*(theta-stheta);
	return abs(a)-cseg;	
	//return a;
}

void main()
{
	vec4 receiverPosition = texture2D(texPosition, gl_TexCoord[0].xy);

	//if(isnan(receiverPosition.x)) 
		//discard;

	vec4 receiverNormal = texture2D(texNormal, gl_TexCoord[0].xy);
	vec4 receiverIrradiance = vec4(0.0); 
	float receiverArea = receiverPosition.w;

	//gl_FragColor = receiverNormal;

	
	// if we are not a leaf surfel
	if (receiverNormal.w > 0.0)
	{
		// first child
		vec2 childIndex = texture2D(texIndices, gl_TexCoord[0].xy).zw;
		
		float totalArea = 0.0;

		// iterate over children gathering radiance
		int numChildren = int(receiverNormal.w);
		for (int i=0; i < numChildren; i++)
		{
			vec4 childPosition = texture2D(texPosition, childIndex);
			vec4 childIrradiance = texture2D(texIrradiance, childIndex);
			
			receiverIrradiance += childIrradiance*abs(childPosition.w);
			totalArea += abs(childPosition.w);

			// go to next emitter
			childIndex = texture2D(texIndices, childIndex).xy;
		}
		
		receiverIrradiance /= totalArea;
	}
	else
	{
		vec2 emitterIndex = vec2(0.0, 0.0);

		// have to double up the loops as there seems to be a 255 iteration limit per loop
		while (emitterIndex.y != 1.0)
		while (emitterIndex.y != 1.0)
		{
			vec4 emitterPosition = texture2D(texPosition, emitterIndex);
			vec4 emitterNext = texture2D(texIndices, emitterIndex);
			float emitterArea = emitterPosition.w;
				
			// vector to emitter
			vec3 v = emitterPosition.xyz - receiverPosition.xyz;
			float dSq = dot(v,v) + 1e-16;		

			// if close enough then traverse children
			if (dSq < 4.0*-emitterArea*kInvPi)
			{
				emitterNext.xy = emitterNext.zw;
				emitterArea = 0.0;
			}

			vec4 emitterNormal = texture2D(texNormal, emitterIndex);
			vec4 emitterIrradiance = texture2D(texIrradiance, emitterIndex);
			vec4 emitterColour = texture2D(texAlbedo, emitterIndex);

			// calc form factor	
			vec3 dir = v * inversesqrt(dSq);
			
			float area = abs(emitterArea);
			
			float ev = dot(dir,emitterNormal.xyz);
			float rv = dot(dir,receiverNormal.xyz);
			//float ff = abs((area*clamp(rv, 0.0, 1.0)*ev)/(kPi*dSq+area));
			float ff = 2.0*kPi*abs(ev*clamp(rv,0.0, 1.0)*(1.0 - inversesqrt(1.0 + ( area / ( dSq*kPi )  ) ) ) );
			
			// add bounced light and glow sources
			if (emitterColour.a > 0.0)
			{
				receiverIrradiance.xyz += ff*emitterColour.xyz*emitterColour.a;
			}
			else if (ev < 0.0)
			{
				receiverIrradiance.xyz += ff*(emitterColour.xyz*kInvPi*emitterIrradiance.xyz);				
			}
			else
			{
				receiverIrradiance.xyz -= ff*emitterIrradiance.xyz*0.8f;
			}
				
			emitterIndex = emitterNext.xy;
		}
	}
	
	gl_FragColor = max(receiverIrradiance, vec4(0.0));
	
}


