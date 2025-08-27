// gaussian_blur_13tap.frag
#version 330 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D image;
uniform bool horizontal;
uniform float texelSize; // 1 / texture width (if horizontal) or 1 / height (if vertical)

// 13-tap Gaussian weights (precomputed with sigma ? 4.0)
float weights[7] = float[](
    0.196482, // weight[0]
    0.176032, // weight[1]
    0.121003, // weight[2]
    0.064759, // weight[3]
    0.027995, // weight[4]
    0.0093,   // weight[5]
    0.002216  // weight[6]
);

void main()
{
    vec2 texOffset = horizontal ? vec2(texelSize, 0.0) : vec2(0.0, texelSize);
    vec3 result = texture(image, TexCoords).rgb * weights[0];

    for (int i = 1; i < 7; ++i)
    {
        result += texture(image, TexCoords + texOffset * i).rgb * weights[i];
        result += texture(image, TexCoords - texOffset * i).rgb * weights[i];
    }

    FragColor = vec4(result, 1.0);
}
