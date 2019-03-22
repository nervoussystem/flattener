#version 120

precision mediump float;

varying vec3 vTransformedNormal;
varying vec3 vColor;
varying vec4 vPosition;

uniform samplerCube envMap;
    

uniform vec3 ambientLightingColor;
uniform vec3 directionalDiffuseColor;
uniform vec3 lightingDirection;
uniform float reflectivity;
uniform float materialShininess;
        
void main(void) {
  vec3 ambientLightWeighting = ambientLightingColor;
  vec3 normal = normalize(vTransformedNormal);
  float diffuseLightBrightness = max(dot(normal,lightingDirection),0.0);
  vec3 diffuseLightWeighting = diffuseLightBrightness*directionalDiffuseColor;
  vec3 r = -reflect(lightingDirection, normal);
  vec3 v = -vPosition.xyz;
  v = normalize(v);
  vec3 reflectEye = (-reflect(v, normal));
  vec3 cameraEye = vec3(0,1,0);
  
  vec3 envColor = textureCube(envMap, reflectEye).rgb;
  vec4 spec = vec4(.3,.3,.3,1.0)*pow(max(0.0,dot(r, v)), materialShininess);
  gl_FragColor = vec4(ambientLightWeighting * vColor+ envColor*reflectivity
                      + vColor * diffuseLightWeighting,
                      1.0)+spec;
}