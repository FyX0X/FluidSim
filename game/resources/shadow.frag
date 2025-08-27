#version 330 core
out vec4 FragColor;

uniform vec3 uShadowColor;

void main()
{
	//FragColor = vec4(uShadowColor, 1.0); // Default to black for the shadow mask
	
	// Set shadow color black for now
	FragColor = vec4(0, 0, 0, 1.0);
}