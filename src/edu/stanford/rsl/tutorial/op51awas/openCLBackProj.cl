typedef float TvoxelValue;
typedef float Tcoord_dev;

// texture sampling
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void openCLBackProj(__read_only image2d_t clSinoTex, __global TvoxelValue* gRes, __constant Tcoord_dev* gVolumeSize) {
	
	
	
}