
void main()
{

	gl_Position = gl_ModelViewProjectionMatrix*gl_Vertex;
	gl_TexCoord[0] = gl_Vertex;
	gl_TexCoord[1].xyz = gl_Normal.xyz;
	gl_TexCoord[2].xy = gl_MultiTexCoord0.xy;
}

