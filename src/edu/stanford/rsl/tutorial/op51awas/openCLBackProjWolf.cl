typedef float TvoxelValue;
typedef float Tcoord_dev;

__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void openCLBackProj(__read_only image2d_t clSinoTex, __global TvoxelValue* gRes, __constant float* resSpac, __constant float* resOrigin, __constant Tcoord_dev* gImgSize) 
{
	int gidx = get_group_id(0);
	int gidy = get_group_id(1);
	int lidx = get_local_id(0);
	int lidy = get_local_id(1);
	
	int locSizex = get_local_size(0);
	int locSizey = get_local_size(1);
	
	int x = mad24(gidx, locSizex, lidx);
	int y = mad24(gidy, locSizey, lidy);
	
	float resSpacX = resSpac[0];
	float resSpacY = resSpac[1];
	float resOrig = resOrigin[0];
	
	// compute world coordinates of current result image position
	float worldX = x * resSpacX + resOrig;
	float worldY = y * resSpacX + resOrig;
	
	unsigned int yStride = gImgSize[0];
	unsigned int ySize = gImgSize[1];
	
	float detectorSpacing = resSpacX;
	float angleSpacing = resSpacY;
	
	float detectorOrig = resOrig;
	
	if (x >= gImgSize[0] || y >= gImgSize[0]) {
		return;
	}
	
	unsigned long idx = y * yStride + x;
	float pixVal = 0;
	
	for (int i = 0; i < ySize; i++) {
	
		// compute angle of detector
		float realAngle = i * angleSpacing;
		float angleX = cos(realAngle);
		float angleY = sin(realAngle);
		
		// compute inner product of world coordinates and angle position of detector
		float2 worldCoord = (float2)(worldX, worldY);
		float2 normVec = (float2)(angleX, angleY);
		float innerPro = dot(worldCoord, normVec);
	
		float dist = (innerPro - detectorOrig) / detectorSpacing;
		
		float val = read_imagef(clSinoTex, sampler, (float2)(dist+0.5f, (float)i+0.5f)).x;
		
		pixVal += val;
	
	}
	
	gRes[idx] = pixVal;
	return;
	
}