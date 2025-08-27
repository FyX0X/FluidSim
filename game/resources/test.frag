#version 330 core

out vec4 FragColor;
in vec2 TexCoords; // Texture coordinates for the fragment

void main()
{
	// Set the fragment color to red for testing
	FragColor = vec4(1.0, 0.0, 1.0, 1.0); // PINK color
}