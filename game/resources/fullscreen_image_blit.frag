#version 330 core

out vec4 FragColor;
in vec2 TexCoords; // Texture coordinates for the fragment

uniform sampler2D uTexture; // Texture to sample from

void main()
{
	FragColor = texture(uTexture, TexCoords); // Sample the texture at the given coordinates
}