#!/bin/sh
glslangValidator -V PackXYMega.comp -o PackXYMega.comp.spv -DGLSL -I.. -I/mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/LiteMath 
glslangValidator -V CastSingleRayMega.comp -o CastSingleRayMega.comp.spv -DGLSL -I.. -I/mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/LiteMath 
glslangValidator -V NaivePathTraceMega.comp -o NaivePathTraceMega.comp.spv -DGLSL -I.. -I/mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/LiteMath 
glslangValidator -V PathTraceMega.comp -o PathTraceMega.comp.spv -DGLSL -I.. -I/mnt/c/Users/Blueberry_iScream/source/repos/kernel_slicer/apps/LiteMath 
