////////////////////////////////////////////////////
//// input file: /mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/HydraCore3/integrator_pt1.cpp
////////////////////////////////////////////////////
#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

void Integrator::InitRandomGens(int a_maxThreads)
{
  m_randomGens.resize(a_maxThreads);
  #pragma omp parallel for default(shared)
  for(int i=0;i<a_maxThreads;i++)
    m_randomGens[i] = RandomGenInit(i);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::kernel_InitEyeRay(uint tid, const uint* packedXY, float4* rayPosAndNear, float4* rayDirAndFar) // (tid,tidX,tidY,tidZ) are SPECIAL PREDEFINED NAMES!!!
{
  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  float3 rayDir = EyeRayDir(x, y, m_winWidth, m_winHeight, m_projInv); 
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, 
                  &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, MAXFLOAT);
}

void Integrator::kernel_InitEyeRay2(uint tid, const uint* packedXY, 
                                   float4* rayPosAndNear, float4* rayDirAndFar,
                                   float4* accumColor,    float4* accumuThoroughput,
                                   RandomGen* gen, uint* rayFlags) // 
{
  *accumColor        = make_float4(0,0,0,1);
  *accumuThoroughput = make_float4(1,1,1,1);
  RandomGen genLocal = m_randomGens[tid];
  *rayFlags          = 0;

  const uint XY = packedXY[tid];

  const uint x = (XY & 0x0000FFFF);
  const uint y = (XY & 0xFFFF0000) >> 16;

  const float2 pixelOffsets = rndFloat2_Pseudo(&genLocal) - float2(0.5f);
  const float fx = float(x) + pixelOffsets.x;
  const float fy = float(y) + pixelOffsets.y;

  float3 rayDir = EyeRayDir(fx, fy, m_winWidth, m_winHeight, m_projInv); 
  float3 rayPos = float3(0,0,0);

  transform_ray3f(m_worldViewInv, &rayPos, &rayDir);
  
  *rayPosAndNear = to_float4(rayPos, 0.0f);
  *rayDirAndFar  = to_float4(rayDir, MAXFLOAT);
  *gen           = genLocal;
}

//uint* temp_out_color;
bool Integrator::kernel_RayTrace(uint tid, const float4* rayPosAndNear, float4* rayDirAndFar,
                                Lite_Hit* out_hit, float2* out_bars)
{
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  CRT_Hit hit = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);
  Lite_Hit res;
  res.primId = hit.primId;
  res.instId = hit.instId;
  res.geomId = hit.geomId;
  res.t      = hit.t;
  
  float2 baricentrics = float2(hit.coords[0], hit.coords[1]);

  //float3 color = float3(hit.coords[2], hit.coords[3], 0.f);
  //if (color.x > 1.0f) color.x = 1.0f;
  //if (color.y > 1.0f) color.y = 1.0f;
 // if (color.z > 1.0f) color.z = 1.0f;
 // const uint XY = m_packedXY[tid];
 // const uint x = (XY & 0x0000FFFF);
 // const uint y = (XY & 0xFFFF0000) >> 16;

  //temp_out_color[y * m_winWidth + x] = RealColorToUint32_f3(color);

  *out_hit  = res;
  *out_bars = baricentrics;
  return (res.primId != -1);
}

void Integrator::kernel_RayTrace2(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar,
                                 float4* out_hit1, float4* out_hit2, uint* rayFlags)
{
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const float4 rayPos = *rayPosAndNear;
  const float4 rayDir = *rayDirAndFar ;

  const CRT_Hit hit   = m_pAccelStruct->RayQuery_NearestHit(rayPos, rayDir);

  if(hit.geomId != uint32_t(-1))
  {
    const float2 uv       = float2(hit.coords[0], hit.coords[1]);
    const float3 hitPos   = to_float3(rayPos) + hit.t*to_float3(rayDir);

    const uint triOffset  = m_matIdOffsets[hit.geomId];
    const uint vertOffset = m_vertOffset  [hit.geomId];
  
    const uint A = m_triIndices[(triOffset + hit.primId)*3 + 0];
    const uint B = m_triIndices[(triOffset + hit.primId)*3 + 1];
    const uint C = m_triIndices[(triOffset + hit.primId)*3 + 2];
  
    const float3 A_norm = to_float3(m_vNorm4f[A + vertOffset]);
    const float3 B_norm = to_float3(m_vNorm4f[B + vertOffset]);
    const float3 C_norm = to_float3(m_vNorm4f[C + vertOffset]);

    const float2 A_texc = m_vTexc2f[A + vertOffset];
    const float2 B_texc = m_vTexc2f[B + vertOffset];
    const float2 C_texc = m_vTexc2f[C + vertOffset];
      
    float3 hitNorm     = (1.0f - uv.x - uv.y)*A_norm + uv.y*B_norm + uv.x*C_norm;
    float2 hitTexCoord = (1.0f - uv.x - uv.y)*A_texc + uv.y*B_texc + uv.x*C_texc;
  
    // transform surface point with matrix and flip normal if needed
    //
    hitNorm = normalize(mul3x3(m_normMatrices[hit.instId], hitNorm));
    const float flipNorm = dot(to_float3(rayDir), hitNorm) > 0.001f ? -1.0f : 1.0f;
    hitNorm = flipNorm*hitNorm;
  
    *rayFlags  = packMatId(currRayFlags, m_matIdByPrimId[m_matIdOffsets[hit.geomId] + hit.primId]);
    *out_hit1  = to_float4(hitPos,  hitTexCoord.x); 
    *out_hit2  = to_float4(hitNorm, hitTexCoord.y); 
  }
  else
    *rayFlags = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_OUT_OF_SCENE) ;
}

void Integrator::kernel_PackXY(uint tidX, uint tidY, uint* out_pakedXY)
{
  const uint inBlockIdX = tidX % 8; // 8x8 blocks
  const uint inBlockIdY = tidY % 8; // 8x8 blocks
 
  const uint localIndex = inBlockIdY*8 + inBlockIdX;
  const uint wBlocks    = m_winWidth/8;

  const uint blockX     = tidX/8;
  const uint blockY     = tidY/8;
  const uint offset     = (blockX + blockY*wBlocks)*8*8 + localIndex;

  out_pakedXY[offset] = ((tidY << 16) & 0xFFFF0000) | (tidX & 0x0000FFFF);
}

void Integrator::kernel_RealColorToUint32(uint tid, float4* a_accumColor, uint* out_color)
{
  out_color[tid] = RealColorToUint32(*a_accumColor);
}

void Integrator::kernel_GetRayColor(uint tid, const Lite_Hit* in_hit, const uint* in_pakedXY, uint* out_color)
{ 
  const Lite_Hit lhit = *in_hit;
  if(lhit.geomId == -1)
  {
    out_color[tid] = 0;
    return;
  }

  const uint32_t matId = m_matIdByPrimId[m_matIdOffsets[lhit.geomId] + lhit.primId];
  const float4 mdata   = m_materials[matId].baseColor;
  //const float3 color = mdata.w > 0.0f ? clamp(float3(mdata.w, mdata.w, mdata.w), 0.0f, 1.0f) : to_float3(mdata);
  float3 temp_color = mdata.w > 0.0f ? clamp(float3(mdata.w,mdata.w,mdata.w), 0.0f, 1.0f) : to_float3(mdata);
  if (lhit.instId == 0)
      temp_color = float3(0.1f + 0.9f * (lhit.instId & 1), 0.1f + 0.8f * (lhit.instId & 2), 0.1f + 0.8f * (lhit.instId & 4));
      //if (tid == 7111u || tid == 7119u)
      //    temp_color = float3(1.f, 1.f, 1.f);
  const float3 color = temp_color;
  //

  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;

  out_color[y*m_winWidth+x] += RealColorToUint32_f3(color);   ////////
}



void Integrator::kernel_SampleLightSource(uint tid, const float4* rayPosAndNear, const float4* rayDirAndFar, 
                                         const float4* in_hitPart1, const float4* in_hitPart2, 
                                         const uint* rayFlags,  
                                         RandomGen* a_gen, float4* out_shadeColor)
{
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);
  const float3 ray_dir = to_float3(*rayDirAndFar);
  
  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;

  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.uv   = float2(data1.w, data2.w);

  const float2 uv = rndFloat2_Pseudo(a_gen);
  
  const float2 sampleOff = 2.0f*(float2(-0.5f,-0.5f) + uv)*m_light.size;
  const float3 samplePos = to_float3(m_light.pos) + float3(sampleOff.x, -1e-5f, sampleOff.y);
  const float  hitDist   = sqrt(dot(hit.pos - samplePos, hit.pos - samplePos));

  const float3 shadowRayDir = normalize(samplePos - hit.pos); // explicitSam.direction;
  const float3 shadowRayPos = hit.pos + hit.norm*std::max(maxcomp(hit.pos), 1.0f)*5e-6f; // TODO: see Ray Tracing Gems

  const bool inShadow = m_pAccelStruct->RayQuery_AnyHit(to_float4(shadowRayPos, 0.0f), to_float4(shadowRayDir, hitDist*0.9995f));
  
  if(!inShadow && dot(shadowRayDir, to_float3(m_light.norm)) < 0.0f)
  {
    const float lightPickProb = 1.0f;
    const float pdfA          = 1.0f / (4.0f*m_light.size.x*m_light.size.y);
    const float cosVal        = std::max(dot(shadowRayDir, (-1.0f)*to_float3(m_light.norm)), 0.0f);
    const float lgtPdfW       = PdfAtoW(pdfA, hitDist, cosVal);
    const BsdfEval bsdfV      = MaterialEval(matId, shadowRayDir, (-1.0f)*ray_dir, hit.norm, hit.uv);
    const float cosThetaOut   = std::max(dot(shadowRayDir, hit.norm), 0.0f);

    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT)
      misWeight = misWeightHeuristic(lgtPdfW, bsdfV.pdf);

    if(cosVal <= 0.0f)
      *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 0.0f);
    else
      *out_shadeColor = to_float4((1.0f/lightPickProb)*(to_float3(m_light.intensity)*bsdfV.color/lgtPdfW)*cosThetaOut*misWeight, 0.0f);
  }
  else
    *out_shadeColor = float4(0.0f, 0.0f, 0.0f, 1.0f);
}

void Integrator::kernel_NextBounce(uint tid, uint bounce, const float4* in_hitPart1, const float4* in_hitPart2, const float4* in_shadeColor, 
                                  float4* rayPosAndNear, float4* rayDirAndFar, float4* accumColor, float4* accumThoroughput, RandomGen* a_gen, MisData* misPrev, uint* rayFlags)
{
  const uint currRayFlags = *rayFlags;
  if(isDeadRay(currRayFlags))
    return;
    
  const uint32_t matId = extractMatId(currRayFlags);

  // process surcase hit case
  //
  const float3 ray_dir = to_float3(*rayDirAndFar);
  const float3 ray_pos = to_float3(*rayPosAndNear);
  
  const float4 data1 = *in_hitPart1;
  const float4 data2 = *in_hitPart2;
  
  SurfaceHit hit;
  hit.pos  = to_float3(data1);
  hit.norm = to_float3(data2);
  hit.uv   = float2(data1.w, data2.w);
  
  const MisData prevBounce = *misPrev;
  const float   prevPdfW   = prevBounce.matSamplePdf;
  const float   prevPdfA   = (prevPdfW >= 0.0f) ? PdfWtoA(prevPdfW, length(ray_pos - hit.norm), prevBounce.cosTheta) : 1.0f;

  // process light hit case
  //
  if(m_materials[matId].brdfType == BRDF_TYPE_LIGHT_SOURCE)
  {
    const float lightIntensity = m_materials[matId].baseColor.w;
    float misWeight = 1.0f;
    if(m_intergatorType == INTEGRATOR_MIS_PT) 
    {
      if(bounce > 0)
      {
        const int lightId   = 0; // #TODO: get light id from material info
        const float lgtPdf  = LightPdfSelectRev(lightId)*LightEvalPDF(lightId, ray_pos, ray_dir, &hit);
        misWeight           = misWeightHeuristic(prevPdfW, lgtPdf);
        if (prevPdfW <= 0.0f) // specular bounce
          misWeight = 1.0f;
      }
    }
    else if(m_intergatorType == INTEGRATOR_SHADOW_PT && hasNonSpecular(currRayFlags))
      misWeight = 0.0f;
    
    float4 currAccumColor       = *accumColor;
    float4 currAccumThoroughput = *accumThoroughput;
    
    const float lightDirectionAtten = dot(to_float3(*rayDirAndFar), float3(0,-1,0)) < 0.0f ? 1.0f : 0.0f;
    
    currAccumColor.x += currAccumThoroughput.x*lightIntensity*misWeight*lightDirectionAtten;
    currAccumColor.y += currAccumThoroughput.y*lightIntensity*misWeight*lightDirectionAtten;
    currAccumColor.z += currAccumThoroughput.z*lightIntensity*misWeight*lightDirectionAtten;
    if(bounce > 0)
      currAccumColor.w *= prevPdfA;
    
    *accumColor = currAccumColor;
    *rayFlags   = currRayFlags | (RAY_FLAG_IS_DEAD | RAY_FLAG_HIT_LIGHT);
    return;
  }
  
  const float4 uv         = rndFloat4_Pseudo(a_gen);
  const BsdfSample matSam = MaterialSampleAndEval(matId, uv, (-1.0f)*ray_dir, hit.norm, hit.uv);
  const float3 bxdfVal    = matSam.color * (1.0f / std::max(matSam.pdf, 1e-10f));
  const float  cosTheta   = dot(matSam.direction, hit.norm);

  MisData nextBounceData;                   // remember current pdfW for next bounce
  nextBounceData.matSamplePdf = matSam.pdf; //
  nextBounceData.cosTheta     = cosTheta;   //
  *misPrev = nextBounceData;                //

  if(m_intergatorType == INTEGRATOR_STUPID_PT)
  {
    *accumThoroughput *= cosTheta*to_float4(bxdfVal, 0.0f); 
  }
  else if(m_intergatorType == INTEGRATOR_SHADOW_PT || m_intergatorType == INTEGRATOR_MIS_PT)
  {
    const float4 currThoroughput = *accumThoroughput;
    const float4 shadeColor      = *in_shadeColor;
    float4 currAccumColor        = *accumColor;

    currAccumColor.x += currThoroughput.x * shadeColor.x;
    currAccumColor.y += currThoroughput.y * shadeColor.y;
    currAccumColor.z += currThoroughput.z * shadeColor.z;
    if(bounce > 0)
      currAccumColor.w *= prevPdfA;

    *accumColor       = currAccumColor;
    *accumThoroughput = currThoroughput*cosTheta*to_float4(bxdfVal, 0.0f); 
  }

  *rayPosAndNear = to_float4(OffsRayPos(hit.pos, hit.norm, matSam.direction), 0.0f);
  *rayDirAndFar  = to_float4(matSam.direction, MAXFLOAT);
  *rayFlags      = currRayFlags | matSam.flags;
}

void Integrator::kernel_ContributeToImage(uint tid, const float4* a_accumColor, const RandomGen* gen, const uint* in_pakedXY, float4* out_color)
{
  const uint XY = in_pakedXY[tid];
  const uint x  = (XY & 0x0000FFFF);
  const uint y  = (XY & 0xFFFF0000) >> 16;
 
  out_color[y*m_winWidth+x] += *a_accumColor;
  m_randomGens[tid] = *gen;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::PackXY(uint tidX, uint tidY)
{
  kernel_PackXY(tidX, tidY, m_packedXY.data());
}

void Integrator::CastSingleRay(uint tid, uint* out_color)
{
  float4 rayPosAndNear, rayDirAndFar;
  kernel_InitEyeRay(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar);
  //temp_out_color = out_color;
  Lite_Hit hit; 
  float2   baricentrics; 
  if(!kernel_RayTrace(tid, &rayPosAndNear, &rayDirAndFar, &hit, &baricentrics))
    return;

  kernel_GetRayColor(tid, &hit, m_packedXY.data(), out_color);
}

void Integrator::NaivePathTrace(uint tid, uint a_maxDepth, float4* out_color)
{
  float4 accumColor, accumThoroughput;
  float4 rayPosAndNear, rayDirAndFar;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThoroughput, &gen, &rayFlags);

  for(int depth = 0; depth < a_maxDepth; depth++) 
  {
    float4   shadeColor, hitPart1, hitPart2;
    kernel_RayTrace2(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThoroughput, &gen, &mis, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
  }

  //kernel_HitEnvironment(tid, &rayFlags, &rayDirAndFar, &accumColor);

  kernel_ContributeToImage(tid, &accumColor, &gen, m_packedXY.data(), 
                           out_color);
}

void Integrator::PathTrace(uint tid, uint a_maxDepth, float4* out_color)
{
  float4 accumColor, accumThoroughput;
  float4 rayPosAndNear, rayDirAndFar;
  RandomGen gen; 
  MisData   mis;
  uint      rayFlags;
  kernel_InitEyeRay2(tid, m_packedXY.data(), &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThoroughput, &gen, &rayFlags);

  for(int depth = 0; depth < a_maxDepth; depth++) 
  {
    float4   shadeColor, hitPart1, hitPart2;
    kernel_RayTrace2(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
    
    kernel_SampleLightSource(tid, &rayPosAndNear, &rayDirAndFar, &hitPart1, &hitPart2, &rayFlags, 
                             &gen, &shadeColor);

    kernel_NextBounce(tid, depth, &hitPart1, &hitPart2, &shadeColor,
                      &rayPosAndNear, &rayDirAndFar, &accumColor, &accumThoroughput, &gen, &mis, &rayFlags);
    if(isDeadRay(rayFlags))
      break;
  }

  kernel_ContributeToImage(tid, &accumColor, &gen, m_packedXY.data(), 
                           out_color);
                           
}
////////////////////////////////////////////////////
//// input file: /mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/HydraCore3/integrator_pt2.cpp
////////////////////////////////////////////////////
#include "integrator_pt.h"
#include "include/crandom.h"

#include <chrono>
#include <string>

float Integrator::LightPdfSelectRev(int a_lightId) 
{ 
  return 1.0f; 
}

float Integrator::LightEvalPDF(int a_lightId, float3 illuminationPoint, float3 ray_dir, const SurfaceHit* pSurfaceHit)
{
  const float3 lpos   = pSurfaceHit->pos;
  const float3 lnorm  = pSurfaceHit->norm;
  const float hitDist = length(illuminationPoint - lpos);
  const float pdfA    = 1.0f / (4.0f*m_light.size.x*m_light.size.y);
  const float cosVal  = std::max(dot(ray_dir, -1.0f*lnorm), 0.0f);
  return PdfAtoW(pdfA, hitDist, cosVal);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BsdfSample Integrator::MaterialSampleAndEval(int a_materialId, float4 rands, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ m_materials[a_materialId].texId[0] ]->sample(texCoordT));

  const uint   type      = m_materials[a_materialId].brdfType;
  const float3 color     = to_float3(m_materials[a_materialId].baseColor)*texColor;
  const float3 specular  = to_float3(m_materials[a_materialId].metalColor);
  const float3 coat      = to_float3(m_materials[a_materialId].coatColor);
  const float  roughness = 1.0f - m_materials[a_materialId].glosiness;
  float  alpha           = m_materials[a_materialId].alpha;
  
  // TODO: read color     from texture
  // TODO: read roughness from texture
  // TODO: read alpha     from texture

  // TODO: check if glosiness in 1 (roughness is 0), use special case mirror brdf
  //if(roughness == 0.0f && type == BRDF_TYPE_GGX)
  //  type = BRDF_TYPE_MIRROR;

  BsdfSample res;
  switch(type)
  {
    case BRDF_TYPE_GLTF:
    case BRDF_TYPE_GGX:
    case BRDF_TYPE_LAMBERT:
    default:
    {
      const float3 ggxDir = ggxSample(float2(rands.x, rands.y), v, n, roughness);
      const float  ggxPdf = ggxEvalPDF (ggxDir, v, n, roughness); 
      const float  ggxVal = ggxEvalBSDF(ggxDir, v, n, roughness);
      
      const float3 lambertDir = lambertSample(float2(rands.x, rands.y), v, n);
      const float  lambertPdf = lambertEvalPDF(lambertDir, v, n);
      const float  lambertVal = lambertEvalBSDF(lambertDir, v, n);

      const float3 h = normalize(v - ggxDir); // half vector.
    
      if(type == BRDF_TYPE_GGX)
        alpha = 1.0f;

      // (1) select between metal and dielectric via rands.z
      //
      float pdfSelect = 1.0f;
      if(rands.z < alpha) // select metall
      {
        pdfSelect *= alpha;
        const float3 F = gltfConductorFresnel(specular, dot(h,v));
        res.direction = ggxDir;
        res.color     = ggxVal*F*alpha;
        res.pdf       = ggxPdf;
      }
      else                // select dielectric
      {
        pdfSelect *= 1.0f - alpha;
        
        // (2) now select between specular and diffise via rands.w
        //
        float fDielectric = gltfFresnelMix2(dot(h,v));
        if(type == BRDF_TYPE_LAMBERT)
          fDielectric = 0.0f;

        if(rands.w < fDielectric) // specular
        {
          pdfSelect *= fDielectric;
          res.direction = ggxDir;
          res.color     = ggxVal*coat*fDielectric*(1.0f - alpha);
          res.pdf       = ggxPdf;
        } 
        else
        {
          pdfSelect *= (1.0f-fDielectric); // lambert
          res.direction = lambertDir;
          res.color     = lambertVal*color*(1.0f-fDielectric)*(1.0f - alpha);
          res.pdf       = lambertPdf;
        }
      }
      
      res.pdf *= pdfSelect;
      res.flags = RAY_FLAG_HAS_NON_SPEC;
    }
    break;
    case BRDF_TYPE_MIRROR:
    {
      res.direction = reflect(v, n);
      // BSDF is multiplied (outside) by cosThetaOut.
      // For mirrors this shouldn't be done, so we pre-divide here instead.
      //
      const float cosThetaOut = dot(res.direction, n);
      res.color     = specular*(1.0f/std::max(cosThetaOut, 1e-6f));
      res.pdf       = 1.0f;
      res.flags     = 0;
    }
    break;
  }

  return res;
}

BsdfEval Integrator::MaterialEval(int a_materialId, float3 l, float3 v, float3 n, float2 tc)
{
  const float2 texCoordT = mulRows2x4(m_materials[a_materialId].row0[0], m_materials[a_materialId].row1[0], tc);
  const float3 texColor  = to_float3(m_textures[ m_materials[a_materialId].texId[0] ]->sample(texCoordT));

  const uint type       = m_materials[a_materialId].brdfType;
  const float3 color    = to_float3(m_materials[a_materialId].baseColor)*texColor;
  const float3 specular = to_float3(m_materials[a_materialId].metalColor);
  const float3 coat     = to_float3(m_materials[a_materialId].coatColor);
  const float roughness = 1.0f - m_materials[a_materialId].glosiness;
        float  alpha    = m_materials[a_materialId].alpha;

  // TODO: read color     from texture
  // TODO: read roughness from texture
  // TODO: read alpha     from texture
 
  // TODO: check if glosiness in 1 (roughness is 0), use special case mirror brdf
  //if(roughness == 0.0f && type == BRDF_TYPE_GGX)
  //  type = BRDF_TYPE_MIRROR;


  BsdfEval res;
  switch(type)
  {
    case BRDF_TYPE_GLTF:
    case BRDF_TYPE_GGX:
    case BRDF_TYPE_LAMBERT:
    default:
    {
      if(type == BRDF_TYPE_GGX)
        alpha = 1.0f;
        
      const float ggxVal = ggxEvalBSDF(l, v, n, roughness);
      const float ggxPdf = ggxEvalPDF (l, v, n, roughness);
      
      const float lambertVal = lambertEvalBSDF(l, v, n);
      const float lambertPdf = lambertEvalPDF (l, v, n);
      
      const float3 h = normalize(v + l);
      const float3 F = gltfConductorFresnel(specular, dot(h,v));

      const float3 specularColor = ggxVal*F;                  // (1) eval metal and (same) specular component
      float  fDielectric         = gltfFresnelMix2(dot(h,v)); // (2) eval dielectric component
      if(type == BRDF_TYPE_LAMBERT)
        fDielectric = 0.0f;
      const float  dielectricPdf = (1.0f-fDielectric)*lambertPdf       + fDielectric*ggxPdf;
      const float3 dielectricVal = (1.0f-fDielectric)*lambertVal*color + fDielectric*ggxVal*coat;

      res.color = alpha*specularColor + (1.0f - alpha)*dielectricVal; // (3) accumulate final color and pdf
      res.pdf   = alpha*ggxPdf        + (1.0f - alpha)*dielectricPdf; // (3) accumulate final color and pdf
    }
    break;
    case BRDF_TYPE_MIRROR:
    {
      res.color = float3(0,0,0);
      res.pdf   = 0.0f;
    }
    break;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Integrator::PackXYBlock(uint tidX, uint tidY, uint a_passNum)
{
  //#pragma omp parallel for default(shared)
  for(int y=0;y<tidY;y++)
    for(int x=0;x<tidX;x++)
      PackXY(x, y);
}

void Integrator::CastSingleRayBlock(uint tid, uint* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
 // std::cout << "START " << start. << std::end;
  #pragma omp parallel for default(shared)
    for (long i = 0; i < tid; i++) {
        //auto start2 = std::chrono::high_resolution_clock::now();
        //if (i >= 100420u && i <= 100421u)
        CastSingleRay(i, out_color);
        //std::cout << "Indiv " << i << std::endl;
        //std::cout << "Indiv " << i << " ray time: " << (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start2).count() / 1000.f) << std::endl;

  }
  naivePtTime = static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()) / 1000.f;
}

void Integrator::NaivePathTraceBlock(uint tid, uint a_maxDepth, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for default(shared)
  for(long i=0;i<tid;i++)
    for(int j=0;j<a_passNum;j++)
      NaivePathTrace(i, a_maxDepth, out_color);
  naivePtTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::PathTraceBlock(uint tid, uint a_maxDepth, float4* out_color, uint a_passNum)
{
  auto start = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for default(shared)
  for(long i=0;i<tid;i++)
    for(int j=0;j<a_passNum;j++)
      PathTrace(i, a_maxDepth, out_color);
  shadowPtTime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count()/1000.f;
}

void Integrator::GetExecutionTime(const char* a_funcName, float a_out[4])
{
  if(std::string(a_funcName) == "NaivePathTrace" || std::string(a_funcName) == "NaivePathTraceBlock")
    a_out[0] = naivePtTime;
  else if(std::string(a_funcName) == "PathTrace" || std::string(a_funcName) == "PathTraceBlock")
    a_out[0] = shadowPtTime;
}
