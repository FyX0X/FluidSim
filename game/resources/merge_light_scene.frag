#version 330 core
in vec2 TexCoords;
out vec4 FragColor;

uniform sampler2D sceneTexture;
uniform sampler2D lightTexture;

void main()
{
    vec3 sceneColor = texture(sceneTexture, TexCoords).rgb;
    vec3 lightColor = texture(lightTexture, TexCoords).rgb;

    // Additive blending of light on scene
    vec3 finalColor = sceneColor * lightColor;

    // Optional: tone mapping or clamp
    //finalColor = finalColor / (finalColor + vec3(1.0));
    //finalColor = pow(finalColor, vec3(1.0 / 2.2)); // gamma correction

    finalColor = clamp(finalColor, 0.0, 1.0);

    FragColor = vec4(finalColor, 1.0);
}
