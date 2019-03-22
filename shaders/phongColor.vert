precision mediump float;

attribute vec3 vertexPosition;
attribute vec3 vertexNormal;
attribute vec4 vertexColor;
    
uniform mat4 mvMatrix;
uniform mat4 pMatrix;
uniform mat3 nMatrix;

varying vec3 vTransformedNormal;
varying vec4 vPosition;
varying vec4 matColor;

void main(void) {
   gl_Position = pMatrix * mvMatrix * vec4(vertexPosition, 1.0);
   vPosition = mvMatrix * vec4(vertexPosition, 1.0);
   vTransformedNormal = nMatrix * vertexNormal;
   matColor = vertexColor;
}
