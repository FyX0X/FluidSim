
#version 330 core
out vec4 FragColor;

in vec2 vTexCoord;
in vec4 vColor;
in float vTexIndex;

uniform sampler2D uTextures[16];

void main()
{
    int index = int(vTexIndex);
    FragColor = vColor * texture(uTextures[index], vTexCoord);
}