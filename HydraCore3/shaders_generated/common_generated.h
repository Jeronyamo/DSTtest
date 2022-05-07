/////////////////////////////////////////////////////////////////////
/////////////  Required  Shader Features ////////////////////////////
/////////////////////////////////////////////////////////////////////
#extension GL_EXT_shader_explicit_arithmetic_types_int8: require

/////////////////////////////////////////////////////////////////////
/////////////////// include files ///////////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
/////////////////// declarations in class ///////////////////////////
/////////////////////////////////////////////////////////////////////
#ifndef uint32_t
#define uint32_t uint
#endif
#define MAXFLOAT 1e37f
#define MINFLOAT 1e37f
#define FLT_MAX 1e37f
#define FLT_MIN -1e37f
#define FLT_EPSILON 1e-6f;
struct Lite_HitT
{
  float t;
  int   primId; 
  int   instId;
  int   geomId;
};
#define Lite_Hit Lite_HitT
struct SurfaceHitT
{
  vec3 pos;
  vec3 norm;
  vec2 uv;
};
#define SurfaceHit SurfaceHitT
struct RectLightSource
{
  vec4 pos;
  vec4 intensity;
  vec4 norm;
  vec2 size;
  vec2 dummy;
};
const float GEPSILON = 5e-6f;
const float DEPSILON = 1e-20f;
struct MisDataT
{
  float matSamplePdf; ///< previous angle pdf (pdfW) that were used for sampling material. if < 0, then material sample was pure specular 
  float cosTheta;     ///< previous dot(matSam.direction, hit.norm)
};
#define MisData MisDataT
struct RandomGenT
{
  uvec2 state;

};
#define RandomGen RandomGenT
struct BsdfSample
{
  vec3 color;
  vec3 direction;
  float  pdf; 
  int    flags;
};
struct BsdfEval
{
  vec3 color;
  float  pdf; 
};
const uint BRDF_TYPE_LAMBERT = 1;
const uint BRDF_TYPE_GGX = 2;
const uint BRDF_TYPE_GLTF = 5;
const uint BRDF_TYPE_GLASS = 6;
const uint BRDF_TYPE_MIRROR = 7;
const uint BRDF_TYPE_LIGHT_SOURCE = 0xEFFFFFFF;
struct GLTFMaterial
{
  vec4 row0[1];     ///< texture matrix
  vec4 row1[1];     ///< texture matrix
  uint   texId[4];    ///< texture id

  vec4 baseColor;   ///< color for both lambert and emissive lights; baseColor.w store emission
  vec4 metalColor;  ///< in our implementation we allow different color for metals and diffuse
  vec4 coatColor;   ///< in our implementation we allow different color for coating (fresnel) and diffuse
  uint  brdfType;     ///<
  uint  lightId;      ///< identifier of light if this material is light  
  float alpha;        ///< blend factor between lambert and reflection : alpha*baseColor + (1.0f-alpha)*baseColor
  float glosiness;    ///< material glosiness or intensity for lights, take color from baseColor
};
struct CRT_Hit 
{
  float    t;         ///< intersection distance from ray origin to object
  uint primId; 
  uint instId;
  uint geomId;    ///< use 4 most significant bits for geometry type; thay are zero for triangles 
  float    coords[4]; ///< custom intersection data; for triangles coords[0] and coords[1] stores baricentric coords (u,v)
};
const uint INTEGRATOR_STUPID_PT = 0;
const uint INTEGRATOR_SHADOW_PT = 1;
const uint INTEGRATOR_MIS_PT = 2;
const uint RAY_FLAG_IS_DEAD = 0x80000000;
const uint RAY_FLAG_OUT_OF_SCENE = 0x40000000;
const uint RAY_FLAG_HIT_LIGHT = 0x20000000;
const uint RAY_FLAG_HAS_NON_SPEC = 0x10000000;

#include "include/Integrator_ubo.h"

/////////////////////////////////////////////////////////////////////
/////////////////// local functions /////////////////////////////////
/////////////////////////////////////////////////////////////////////
bool isfinite(float x) { return !isinf(x); }

void CoordinateSystem(vec3 v1, inout vec3 v2, inout vec3 v3) {
  float invLen = 1.0f;

  if (abs(v1.x) > abs(v1.y))
  {
    invLen = 1.0f / sqrt(v1.x*v1.x + v1.z*v1.z);
    (v2)  = vec3((-1.0f) * v1.z * invLen,0.0f,v1.x * invLen);
  }
  else
  {
    invLen = 1.0f / sqrt(v1.y * v1.y + v1.z * v1.z);
    (v2)  = vec3(0.0f,v1.z * invLen,(-1.0f) * v1.y * invLen);
  }

  (v3) = cross(v1, (v2));
}

float GGX_Distribution(const float cosThetaNH, const float alpha) {
  const float alpha2  = alpha * alpha;
  const float NH_sqr  = clamp(cosThetaNH * cosThetaNH, 0.0f, 1.0f);
  const float den     = NH_sqr * alpha2 + (1.0f - NH_sqr);
  return alpha2 / max(float((M_PI)) * den * den, 1e-6f);
}

float GGX_GeomShadMask(const float cosThetaN, const float alpha) {
  // Height - Correlated G.
  //const float tanNV      = sqrt(1.0f - cosThetaN * cosThetaN) / cosThetaN;
  //const float a          = 1.0f / (alpha * tanNV);
  //const float lambda     = (-1.0f + sqrt(1.0f + 1.0f / (a*a))) / 2.0f;
  //const float G          = 1.0f / (1.0f + lambda);

  // Optimized and equal to the commented-out formulas on top.
  const float cosTheta_sqr = clamp(cosThetaN*cosThetaN, 0.0f, 1.0f);
  const float tan2         = (1.0f - cosTheta_sqr) / max(cosTheta_sqr, 1e-6f);
  const float GP           = 2.0f / (1.0f + sqrt(1.0f + alpha * alpha * tan2));
  return GP;
}

vec3 MapSampleToCosineDistribution(float r1, float r2, vec3 direction, vec3 hit_norm, float power) {
  if(power >= 1e6f)
    return direction;

  const float sin_phi = sin(M_TWOPI * r1);
  const float cos_phi = cos(M_TWOPI * r1);

  //sincos(2.0f*r1*3.141592654f, &sin_phi, &cos_phi);

  const float cos_theta = pow(1.0f - r2, 1.0f / (power + 1.0f));
  const float sin_theta = sqrt(1.0f - cos_theta*cos_theta);

  vec3 deviation;
  deviation.x = sin_theta*cos_phi;
  deviation.y = sin_theta*sin_phi;
  deviation.z = cos_theta;

  vec3 ny = direction,  nx,  nz;
  CoordinateSystem(ny, nx, nz);

  {
    vec3 temp = ny;
    ny = nz;
    nz = temp;
  }

  vec3 res = nx*deviation.x + ny*deviation.y + nz*deviation.z;

  float invSign = dot(direction, hit_norm) > 0.0f ? 1.0f : -1.0f;

  if (invSign*dot(res, hit_norm) < 0.0f) // reflected ray is below surface #CHECK_THIS
  {
    res = (-1.0f)*nx*deviation.x + ny*deviation.y - nz*deviation.z;
    //belowSurface = true;
  }

  return res;
}

vec3 SphericalDirectionPBRT(const float sintheta, const float costheta, const float phi) { 
  return vec3(sintheta * cos(phi),sintheta * sin(phi),costheta); 
}

float misHeuristicPower1(float p) { return isfinite(p) ? abs(p) : 0.0f; }

vec2 mulRows2x4(const vec4 row0, const vec4 row1, vec2 v) {
  vec2 res;
  res.x = row0.x*v.x + row0.y*v.y + row0.w;
  res.y = row1.x*v.x + row1.y*v.y + row1.w;
  return res;
}

vec3 ggxSample(vec2 rands, vec3 v, vec3 n, float roughness) {
  const float roughSqr = roughness * roughness;
    
  vec3 nx,  ny, nz = n;
  CoordinateSystem(nz, nx, ny);
    
  const vec3 wo = vec3(dot(v, nx),dot(v, ny),dot(v, nz));
  const float phi       = rands.x * M_TWOPI;
  const float cosTheta  = clamp(sqrt((1.0f - rands.y) / (1.0f + roughSqr * roughSqr * rands.y - rands.y)), 0.0f, 1.0f);
  const float sinTheta  = sqrt(1.0f - cosTheta * cosTheta);
  const vec3 wh = SphericalDirectionPBRT(sinTheta, cosTheta, phi);
    
  const vec3 wi = 2.0f * dot(wo, wh) * wh - wo;      // Compute incident direction by reflecting about wm  
  return normalize(wi.x * nx + wi.y * ny + wi.z * nz); // back to normal coordinate system
}

float ggxEvalBSDF(vec3 l, vec3 v, vec3 n, float roughness) {
  if(abs(dot(l, n)) < 1e-5f)
    return 0.0f; 
 
  const float dotNV = dot(n, v);  
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 0.0f; 

  const float  roughSqr = roughness * roughness;
  const vec3 h = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  const float G     = GGX_GeomShadMask(dotNV, roughSqr)*GGX_GeomShadMask(dotNL, roughSqr);      

  return (D * G / max(4.0f * dotNV * dotNL, 1e-6f));  // Pass single-scattering
}

float epsilonOfPos(vec3 hitPos) { return max(max(abs(hitPos.x), max(abs(hitPos.y), abs(hitPos.z))), 2.0f*GEPSILON)*GEPSILON; }

vec3 lambertSample(vec2 rands, vec3 v, vec3 n) {
   return MapSampleToCosineDistribution(rands.x, rands.y, n, n, 1.0f);
}

float PdfAtoW(const float aPdfA, const float aDist, const float aCosThere) {
  return (aPdfA*aDist*aDist) / max(aCosThere, 1e-30f);
}

float lambertEvalPDF(vec3 l, vec3 v, vec3 n) { 
  return abs(dot(l, n)) * INV_PI;
}

float ggxEvalPDF(vec3 l, vec3 v, vec3 n, float roughness) { 
  const float dotNV = dot(n, v);
  const float dotNL = dot(n, l);
  if (dotNV < 1e-6f || dotNL < 1e-6f)
    return 1.0f;

  const float  roughSqr  = roughness * roughness;
    
  const vec3 h = normalize(v + l); // half vector.
  const float dotNH = dot(n, h);
  const float dotHV = dot(h, v);
  const float D     = GGX_Distribution(dotNH, roughSqr);
  return  D * dotNH / (4.0f * max(dotHV, 1e-6f));
}

vec3 mul4x3(mat4 m, vec3 v) {
  return (m*vec4(v, 1.0f)).xyz;
}

vec3 gltfConductorFresnel(vec3 f0, float VdotH) {
  const float tmp = 1.0f - abs(VdotH);
  return f0 + (vec3(1.0f,1.0f,1.0f) - f0) * (tmp*tmp*tmp*tmp*tmp);
}

uint NextState(inout RandomGen gen) {
  const uint x = (gen.state).x * 17 + (gen.state).y * 13123;
  (gen.state).x = (x << 13) ^ x;
  (gen.state).y ^= (x << 7);
  return x;
}

float gltfFresnelMix2(float VdotH) {
  //return 0.25f;
  //const float f1  = (1.0f-ior)/(1+ior);
  //const float f0  = f1*f1;
  // Note that the dielectric index of refraction ior = 1.5 is now f0 = 0.04
  const float tmp = 1.0f - abs(VdotH);
  return 0.04f + 0.96f*(tmp*tmp*tmp*tmp*tmp);
}

float lambertEvalBSDF(vec3 l, vec3 v, vec3 n) {
  return INV_PI;
}

vec3 EyeRayDir(float x, float y, float w, float h, mat4 a_mViewProjInv) {
  vec4 pos = vec4(2.0f * (x + 0.5f) / w - 1.0f, -2.0f * (y + 0.5f) / h + 1.0f, 0.0f, 1.0f);

  pos = (a_mViewProjInv*pos);
  pos /= pos.w;

  pos.y *= (-1.0f);

  return normalize(pos.xyz);
}

vec2 rndFloat2_Pseudo(inout RandomGen gen) {
  uint x = NextState(gen); 

  const uint x1 = (x * (x * x * 15731 + 74323) + 871483);
  const uint y1 = (x * (x * x * 13734 + 37828) + 234234);

  const float scale     = (1.0f / 4294967296.0f);

  return vec2(float((x1)), float((y1)))*scale;
}

vec3 mul3x3(mat4 m, vec3 v) { 
  return (m*vec4(v, 0.0f)).xyz;
}

void transform_ray3f(mat4 a_mWorldViewInv, inout vec3 ray_pos, inout vec3 ray_dir) {
  vec3 pos = mul4x3(a_mWorldViewInv, (ray_pos));
  vec3 pos2 = mul4x3(a_mWorldViewInv, ((ray_pos) + 100.0f*(ray_dir)));

  vec3 diff = pos2 - pos;

  (ray_pos)  = pos;
  (ray_dir)  = normalize(diff);
}

float misWeightHeuristic(float a, float b) {
  const float w = misHeuristicPower1(a) / max(misHeuristicPower1(a) + misHeuristicPower1(b), 1e-30f);
  return isfinite(w) ? w : 0.0f;
}

float PdfWtoA(const float aPdfW, const float aDist, const float aCosThere) {
  return aPdfW * abs(aCosThere) / max(aDist*aDist, 1e-30f);
}

vec4 rndFloat4_Pseudo(inout RandomGen gen) {
  uint x = NextState(gen);

  const uint x1 = (x * (x * x * 15731 + 74323) + 871483);
  const uint y1 = (x * (x * x * 13734 + 37828) + 234234);
  const uint z1 = (x * (x * x * 11687 + 26461) + 137589);
  const uint w1 = (x * (x * x * 15707 + 789221) + 1376312589);

  const float scale = (1.0f / 4294967296.0f);

  return vec4(float((x1)), float((y1)), float((z1)), float((w1)))*scale;
}

vec3 OffsRayPos(const vec3 a_hitPos, const vec3 a_surfaceNorm, const vec3 a_sampleDir) {
  const float signOfNormal2 = dot(a_sampleDir, a_surfaceNorm) < 0.0f ? -1.0f : 1.0f;
  const float offsetEps     = epsilonOfPos(a_hitPos);
  return a_hitPos + signOfNormal2*offsetEps*a_surfaceNorm;
}

float maxcomp(vec3 v) { return max(v.x, max(v.y, v.z)); }

uint RealColorToUint32_f3(vec3 real_color) {
  float  r = real_color.x*255.0f;
  float  g = real_color.y*255.0f;
  float  b = real_color.z*255.0f;
  uint8_t red = uint8_t(r), green = uint8_t(g), blue = uint8_t(b);
  return int(red) | (int(green) << 8) | (int(blue) << 16) | 0xFF000000;
}

uint fakeOffset(uint x, uint y, uint pitch) { return y*pitch + x; }  // RTV pattern, for 2D threading

#define KGEN_FLAG_RETURN            1
#define KGEN_FLAG_BREAK             2
#define KGEN_FLAG_DONT_SET_EXIT     4
#define KGEN_FLAG_SET_EXIT_NEGATIVE 8
#define KGEN_REDUCTION_LAST_STEP    16

