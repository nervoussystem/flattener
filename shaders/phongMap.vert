#version 120

precision mediump float;

varying vec3 vTransformedNormal;
varying vec4 vPosition;
varying vec3 vColor;

void main(void) {
   gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
   vColor = vec3(gl_Color);
   vPosition = gl_ModelViewMatrix * gl_Vertex;
   vTransformedNormal = gl_NormalMatrix * gl_Normal;
}