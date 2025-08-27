#version 330 core

out vec4 FragColor;


uniform mat4 invProj;
uniform mat4 invView;

uniform vec2 uLightPos;       // Light position in world coords
uniform vec3 uLightColor;     // Light RGB color
uniform float uLightRadius;   // Light radius (max influence)
uniform float uLightIntensity; // Light intensity scalar


uniform vec2 uScreenSize;      // Needed to convert gl_FragCoord to NDC

void main()
{
    // Get NDC coords from screen coords
    vec2 ndc = (gl_FragCoord.xy / uScreenSize) * 2.0 - 1.0;
    vec4 clipSpace = vec4(ndc, 0.0, 1.0);

    // Convert to world space
    vec4 worldSpace = invView * invProj * clipSpace;
    vec2 worldPos = worldSpace.xy / worldSpace.w;


    // Calculate distance from current fragment to light center
    float dist = length(worldPos - uLightPos);

    // If outside light radius, no contribution
    if (dist > uLightRadius)
        discard;

    // Compute attenuation (inverse quadratic falloff)
    float attenuation = 1.0 - (dist / uLightRadius);
    attenuation = attenuation * attenuation; // quadratic falloff

    // Final light color with intensity and attenuation
    vec3 light = uLightColor * attenuation * uLightIntensity;

    FragColor = vec4(light, 1.0);
}
