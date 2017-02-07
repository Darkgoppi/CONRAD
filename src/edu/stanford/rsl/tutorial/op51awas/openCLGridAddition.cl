typedef float TvoxelValue;
typedef float Tcoord_dev;

// texture sampling
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void openCLGridAdd(__read_only image2d_t clGrid1Tex, __read_only image2d_t clGrid2Tex, __global TvoxelValue* gRes, __constant Tcoord_dev* gVolumeSize) {
	int gidx = get_group_id(0);
	int gidy = get_group_id(1);
	int lidx = get_local_id(0);
	int lidy = get_local_id(1);
	
	int locSizeX = get_local_size(0);
	int locSizeY = get_local_size(1);
	
	int x = mad24(gidx, locSizeX, lidx);
	int y = mad24(gidy, locSizeY, lidy);
	
	unsigned int yStride = gVolumeSize[0];
	
	if (x >= gVolumeSize[0] || y >= gVolumeSize[1]) {
		return;
	}
	
	unsigned long idx = y * yStride + x;
	
	float val1 = read_imagef(clGrid1Tex, sampler, (float2)(x+0.5f, y+0.5f)).x;
	float val2 = read_imagef(clGrid2Tex, sampler, (float2)(x+0.5f, y+0.5f)).x;
	
	gRes[idx] = val1 + val2;
	
	return;
}	