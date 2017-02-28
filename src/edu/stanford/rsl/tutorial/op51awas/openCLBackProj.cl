typedef float TvoxelValue;
typedef float Tcoord_dev;

// texture sampling
__constant sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void openCLBackProj(__read_only image2d_t clSinoTex, __global TvoxelValue* gRes, __constant float* resSpac, __constant float* resOrig, __constant Tcoord_dev* gImgSize) {
	
	int x = get_global_id(0);
	int y = get_global_id(1);
	
	float spacingX = resSpac[0];
	float spacingY = resSpac[1];
	
	float origX = resOrig[0];
	float origY = resOrig[1];
	
	// compute world coordinates of current result image position
	float worldX = x * spacingX + origX;
	float worldY = y * spacingY + origY;
	
	printf("%f \n", spacingX);
	
	unsigned int yStride = gImgSize[0];
	unsigned int ySize = gImgSize[1];
	
	float detectorSpacing = resSpac[2];
	float angleSpacing = resSpac[3];
	
	float detectorOrig = resOrig[2];
	
	if (x >= gImgSize[0] || y >= gImgSize[0]) {
		return;
	}
	
//	printf("%d \n", ySize);
	
	unsigned long idx = y * yStride + x;
	float pixVal = 0;
	
	bool inFor = false;
	
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
		
		if (val > 0) {
			inFor = true;
		}
		
		pixVal += val;
	
	}
	
	if (inFor) {
		printf("Ich habe fertig \n");
	}
	gRes[idx] = pixVal;
//	gRes[idx] /= (gImgSize[1] / M_PI);
//	gRes[idx] = 5;
	
	return;
	
}