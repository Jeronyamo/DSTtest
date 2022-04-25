#include <iostream>
#include <fstream>

#include "integrator_pt.h"
#include "Bitmap.h"
#include "ArgParser.h"

//#include "vk_context.h"
//std::shared_ptr<Integrator> CreateIntegrator_Generated(int a_maxThreads, vk_utils::VulkanContext a_ctx, size_t a_maxThreadsGenerated);

int main(int argc, const char** argv)
{
  #ifndef NDEBUG
  bool enableValidationLayers = true;
  #else
  bool enableValidationLayers = false;
  #endif

  int WIN_WIDTH  = 512;
  int WIN_HEIGHT = 512;

  std::vector<uint32_t> pixelData(WIN_WIDTH*WIN_HEIGHT);
  std::vector<float4>   realColor(WIN_WIDTH*WIN_HEIGHT);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  std::shared_ptr<Integrator> pImpl = nullptr;
  ArgParser args(argc, argv);
  
  std::string scenePath = "../resources/HydraCore/hydra_app/tests/02_cry_sponza/statex_00001.xml";
  if(args.hasOption("-in"))
    scenePath = args.getOptionValue<std::string>("-in");

  std::string imageOut = "z_out.bmp";
  if(args.hasOption("-out"))
    imageOut = args.getOptionValue<std::string>("-out");

  const std::string imageOutClean = imageOut.substr(0, imageOut.find_last_of("."));

  std::string integratorType = "mispt";
  if(args.hasOption("-integrator"))
    integratorType = args.getOptionValue<std::string>("-integrator");

  int PASS_NUMBER = 100;
  if(args.hasOption("-spp"))
    PASS_NUMBER = args.getOptionValue<int>("-spp");
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  bool onGPU = args.hasOption("--gpu");
  if(onGPU)
  {
    unsigned int a_preferredDeviceId = args.getOptionValue<int>("--gpu_id", 0);
  }
  else
    pImpl = std::make_shared<Integrator>(WIN_WIDTH*WIN_HEIGHT);
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////

  pImpl->SetViewport(0,0,WIN_WIDTH,WIN_HEIGHT);
  std::cout << "[main]: Loading scene ... " << scenePath.c_str() << std::endl;
  pImpl->LoadScene(scenePath.c_str());
  pImpl->CommitDeviceData();

  // remember (x,y) coords for each thread to make our threading 1D
  //

  ISceneObject* treePtr = dynamic_cast<ISceneObject*>(pImpl->m_pAccelStruct.get());
  std::cout << "[main]: PackXYBlock() ... " << std::endl;
  pImpl->PackXYBlock(WIN_WIDTH, WIN_HEIGHT, 1);

  const float normConst = 1.0f/float(PASS_NUMBER);
  const float invGamma  = 1.0f/2.2f;
  
  // now test path tracing
  //
  float minValPdf = 1e30f;
  float maxValPdf = -1e29f;
  if(integratorType == "naivept" || integratorType == "all")
  {
    const int NAIVE_PT_REPEAT = 1;
  
    std::cout << "[main]: NaivePathTraceBlock() ... " << std::endl;
    memset(realColor.data(), 0, sizeof(float)*4*realColor.size());
    pImpl->SetIntegratorType(Integrator::INTEGRATOR_STUPID_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->NaivePathTraceBlock(WIN_HEIGHT*WIN_HEIGHT, 6, realColor.data(), PASS_NUMBER*NAIVE_PT_REPEAT);
    
    for(int i=0;i<WIN_HEIGHT*WIN_HEIGHT;i++)
    {
      float4 color = realColor[i]*normConst*(1.0f/float(NAIVE_PT_REPEAT));
      if(std::isfinite(color.w))
      {
        minValPdf  = std::min(minValPdf, color.w);
        maxValPdf  = std::max(maxValPdf, color.w);
      }
      color.x      = std::pow(color.x, invGamma);
      color.y      = std::pow(color.y, invGamma);
      color.z      = std::pow(color.z, invGamma);
      color.w      = 1.0f;
      pixelData[i] = RealColorToUint32(clamp(color, 0.0f, 1.0f));
    }
    
    const std::string outName = (integratorType == "naivept") ? imageOut : imageOutClean + "_naivept.bmp"; 
    SaveBMP(outName.c_str(), pixelData.data(), WIN_WIDTH, WIN_HEIGHT);
  }

  //std::cout << "minValPdf = " << minValPdf << std::endl;
  //std::cout << "maxValPdf = " << maxValPdf << std::endl;
  //const float normInv = 1.0f / (maxValPdf - minValPdf);
  //for(int i=0;i<WIN_HEIGHT*WIN_HEIGHT;i++)
  //{
  //  const float pdf = std::isfinite(realColor[i].w) ? (realColor[i].w*normConst*(1.0f/float(NAIVE_PT_REPEAT)) - minValPdf)*normInv : 0.0f;
  //  float4 color;
  //  color.x      = pdf;
  //  color.y      = std::sqrt(pdf);
  //  color.z      = std::sqrt(pdf);
  //  color.w      = 1.0f;
  //  pixelData[i] = RealColorToUint32(clamp(color, 0.0f, 1.0f));
  //}
  //
  //if(onGPU)
  //  SaveBMP("zout_gpu1.bmp", pixelData.data(), WIN_WIDTH, WIN_HEIGHT);
  //else
  //  SaveBMP("zout_cpu1.bmp", pixelData.data(), WIN_WIDTH, WIN_HEIGHT);

  // -------------------------------------------------------------------------------
  if(integratorType == "shadowpt" || integratorType == "all")
  {
    std::cout << "[main]: PathTraceBlock(Shadow-PT) ... " << std::endl;
    memset(realColor.data(), 0, sizeof(float)*4*realColor.size());
    pImpl->SetIntegratorType(Integrator::INTEGRATOR_SHADOW_PT);
    pImpl->UpdateMembersPlainData();
    pImpl->PathTraceBlock(WIN_HEIGHT*WIN_HEIGHT, 6, realColor.data(), PASS_NUMBER);
  
    for(int i=0;i<WIN_HEIGHT*WIN_HEIGHT;i++)
    {
      float4 color = realColor[i]*normConst;
      color.x      = std::pow(color.x, invGamma);
      color.y      = std::pow(color.y, invGamma);
      color.z      = std::pow(color.z, invGamma);
      color.w      = 1.0f;
      pixelData[i] = RealColorToUint32(clamp(color, 0.0f, 1.0f));
    }
  
    const std::string outName = (integratorType == "shadowpt") ? imageOut : imageOutClean + "_shadowpt.bmp"; 
    SaveBMP(outName.c_str(), pixelData.data(), WIN_WIDTH, WIN_HEIGHT);
  }

  // -------------------------------------------------------------------------------
  if(integratorType == "mispt" || integratorType == "all")
  {
    std::cout << "[main]: PathTraceBlock(MIS-PT) ... " << std::endl;
    //memset(realColor.data(), 0, sizeof(float)*4*realColor.size());
    memset(pixelData.data(), 0u, sizeof(unsigned) * pixelData.size());
    pImpl->SetIntegratorType(Integrator::INTEGRATOR_MIS_PT);
    pImpl->UpdateMembersPlainData();
    //pImpl->PathTraceBlock(WIN_HEIGHT*WIN_HEIGHT, 6, realColor.data(), PASS_NUMBER);
    pImpl->CastSingleRayBlock(WIN_HEIGHT * WIN_HEIGHT, pixelData.data(), PASS_NUMBER);
  /*
    for(int i=0;i<WIN_HEIGHT*WIN_HEIGHT;i++)
    {
      float4 color = realColor[i]*normConst;
      color.x      = std::pow(color.x, invGamma);
      color.y      = std::pow(color.y, invGamma);
      color.z      = std::pow(color.z, invGamma);
      color.w      = 1.0f;
      pixelData[i] = RealColorToUint32(clamp(color, 0.0f, 1.0f));
    }
  */
    const std::string outName = (integratorType == "mispt") ? imageOut : imageOutClean + "_mispt.bmp"; 
    SaveBMP(outName.c_str(), pixelData.data(), WIN_WIDTH, WIN_HEIGHT);
  }
  // -------------------------------------------------------------------------------

  std::cout << std::endl;
  float timings[4] = {0,0,0,0};
  pImpl->GetExecutionTime("NaivePathTraceBlock", timings);
  std::cout << "NaivePathTraceBlock(exec)  = " << timings[0]              << " ms " << std::endl;
  std::cout << "NaivePathTraceBlock(copy)  = " << timings[1] + timings[2] << " ms " << std::endl;
  std::cout << "NaivePathTraceBlock(ovrh)  = " << timings[3]              << " ms " << std::endl;
  std::cout << std::endl;
  pImpl->GetExecutionTime("PathTraceBlock", timings);
  std::cout << "PathTraceBlock(exec) = " << timings[0]              << " ms " << std::endl;
  std::cout << "PathTraceBlock(copy) = " << timings[1] + timings[2] << " ms " << std::endl;
  std::cout << "PathTraceBlock(ovrh) = " << timings[3]              << " ms " << std::endl;
  return 0;
}

