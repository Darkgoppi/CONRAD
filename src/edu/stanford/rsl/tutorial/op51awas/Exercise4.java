package edu.stanford.rsl.tutorial.op51awas;

import java.io.IOException;
import java.io.InputStream;
import java.nio.FloatBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLImage2d;
import com.jogamp.opencl.CLImageFormat;
import com.jogamp.opencl.CLImageFormat.ChannelOrder;
import com.jogamp.opencl.CLImageFormat.ChannelType;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLMemory.Mem;
import com.jogamp.opencl.CLProgram;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGridOperators;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;
import ij.ImageJ;

public class Exercise4 {

	//--------------Add-Shleife-----------------------------
	private static OpenCLGrid2D addLoop(Phantom phantom){
		OpenCLGrid2D testOpenCL = new OpenCLGrid2D(phantom);
		OpenCLGrid2D testOpenCL2 = new OpenCLGrid2D(phantom);
		
		OpenCLGridOperators clop = OpenCLGridOperators.getInstance();
		for (int i = 0; i < 1000; i++) {
			clop.addBy(testOpenCL, testOpenCL2);
		} 
		return testOpenCL;
	}
	
	
	//---------------Add-Kernel-----------------------------
	private static Grid2D addKernel(Phantom phantom) throws IOException {
		
		OpenCLGrid2D testOpenCL3 = new OpenCLGrid2D(phantom);
		OpenCLGrid2D testOpenCL4 = new OpenCLGrid2D(phantom);
		OpenCLGrid2D testOpenCL5 = new OpenCLGrid2D(new Grid2D(phantom.getWidth(), phantom.getHeight()));
		
		int[] imgSize = new int[] {phantom.getWidth(), phantom.getHeight()};
		
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice device = context.getMaxFlopsDevice();
		CLCommandQueue commandQueue = device.createCommandQueue();
		
		InputStream is = Exercise4.class.getResourceAsStream("openCLGridAddWolf.cl");
		CLProgram program = context.createProgram(is).build();
		CLKernel kernelFunction = program.createCLKernel("addKernel");
		
		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(imgSize.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(new float[] {imgSize[0], imgSize[1]});
		gImgSize.getBuffer().rewind();
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		//Grid1
		testOpenCL3.getDelegate().prepareForDeviceOperation();
		testOpenCL3.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> g1Tex = null;
		g1Tex = context.createImage2d(testOpenCL3.getDelegate().getCLBuffer().getBuffer(), imgSize[0], imgSize[1], format, Mem.READ_ONLY);
		testOpenCL3.getDelegate().release();
		
		//Grid2
		testOpenCL4.getDelegate().prepareForDeviceOperation();
		testOpenCL4.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> g2Tex = null;
		g2Tex = context.createImage2d(testOpenCL4.getDelegate().getCLBuffer().getBuffer(), imgSize[0], imgSize[1], format, Mem.READ_ONLY);
		testOpenCL4.getDelegate().release();
		
		//Grid3
		testOpenCL5.getDelegate().prepareForDeviceOperation();
		commandQueue.putWriteImage(g1Tex, true).putWriteImage(g2Tex, true)
		.putWriteBuffer(testOpenCL5.getDelegate().getCLBuffer(), true).putWriteBuffer(gImgSize, true).finish();
		kernelFunction.rewind();
		kernelFunction.putArg(g1Tex).putArg(g2Tex).putArg(testOpenCL5.getDelegate().getCLBuffer()).putArg(gImgSize);
		
		//Grafikkarten-Settings
		int bpBlockSize[] = {32, 32};
		int maxWorkGroupSize = device.getMaxWorkGroupSize();
		int[] realLocalSize = new int[]{
				Math.min((int)Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[0]), 
				Math.min((int)Math.pow(maxWorkGroupSize, 1/2.0), bpBlockSize[1])};
		
		int[] globalWorksize = new int[]{OpenCLUtil.roundUp(realLocalSize[0], imgSize[0]), 
						OpenCLUtil.roundUp(realLocalSize[1], imgSize[1])};
		
		commandQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorksize[0], globalWorksize[1], realLocalSize[0], realLocalSize[1]).finish();
		testOpenCL5.getDelegate().notifyDeviceChange();
		
		Grid2D result = new Grid2D(testOpenCL5);
		return result;
	}
	
	//------------Back-Projector---------------------------
	private static Grid2D backProj(Grid2D sinogram) throws IOException {
		
		OpenCLGrid2D sinogramCL = new OpenCLGrid2D(sinogram);
		OpenCLGrid2D testOpenCL6 = new OpenCLGrid2D(new Grid2D(sinogram.getWidth(), sinogram.getWidth()));
		testOpenCL6.setSpacing(sinogram.getSpacing());
		testOpenCL6.setOrigin(sinogram.getOrigin()[0], sinogram.getOrigin()[0]);
		
		int[] imgSize = new int[] {sinogram.getWidth(), sinogram.getHeight()};
		
		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice dev = context.getMaxFlopsDevice();
		CLCommandQueue comQueue = dev.createCommandQueue();

		CLProgram prog = context.createProgram(Exercise4.class.getResourceAsStream("OpenCLBackProjWolf.cl")).build();
		CLKernel kernelFunction = prog.createCLKernel("openCLBackProj");

		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(imgSize.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(new float[] { imgSize[0], imgSize[1] });
		gImgSize.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> resSpac = context.createFloatBuffer(sinogramCL.getSpacing().length, Mem.READ_ONLY);
		resSpac.getBuffer().put(new float[] { (float)sinogramCL.getSpacing()[0], (float)sinogramCL.getSpacing()[1] });
		resSpac.getBuffer().rewind();
		
		CLBuffer<FloatBuffer> resOrig = context.createFloatBuffer(1, Mem.READ_ONLY);
		resOrig.getBuffer().put(new float[] {(float)sinogramCL.getOrigin()[0]});
		resOrig.getBuffer().rewind();
		
		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY, ChannelType.FLOAT);
		
		//Input-Sinogram
		sinogramCL.getDelegate().prepareForDeviceOperation();
		sinogramCL.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> clSinoTex = context.createImage2d(sinogramCL
				.getDelegate().getCLBuffer().getBuffer(), imgSize[0], imgSize[1],
				format, Mem.READ_ONLY);
		sinogramCL.getDelegate().release();
		
		//Ergebnis-BackProjection
		testOpenCL6.getDelegate().prepareForDeviceOperation();
		comQueue.putWriteImage(clSinoTex, true).putWriteBuffer(testOpenCL6.getDelegate().getCLBuffer(), true).putWriteBuffer(resSpac, true)
											   .putWriteBuffer(resOrig, true).putWriteBuffer(gImgSize, true).finish();
		kernelFunction.rewind();
		kernelFunction.putArg(clSinoTex).putArg(testOpenCL6.getDelegate().getCLBuffer()).putArg(resSpac).putArg(resOrig).putArg(gImgSize);

		//Grafikkarten-Settings
		int bpBlockSize[] = { 32, 32 };
		int maxWorkGroupSize = dev.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),	bpBlockSize[0]),
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),	bpBlockSize[1])};
		
		int[] globalWorkSize = new int[]{OpenCLUtil.roundUp(realLocalSize[0], imgSize[0]), 
				OpenCLUtil.roundUp(realLocalSize[1], imgSize[0])};
		
		comQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
		testOpenCL6.getDelegate().notifyDeviceChange();
		
		return (new Grid2D(testOpenCL6));
	}
	
	//------------------------------------------------------
	//---------------MAIN-----------------------------------
	//------------------------------------------------------
	public static void main(String[] args) throws IOException{
				
		new ImageJ();
		
		//--------Variablen--------
		//Phantom
		int width = 200;
		int height = 200;
		double spacingX = 1.0;		//0.5
		double spacingY = 1.0;		//0.5
		
		//BackProj
		int numProj = 200;
		int detSize = 250;
		double detSpacing = 1.0;		//0.75, 1.25 
		double angularRange = Math.PI;
		double angularInc = angularRange/numProj;
		double samplingRate = 0.05;
		
		//------Aufgaben-----------
		Phantom phantom = new Phantom(width, height, spacingX, spacingY);
		phantom.show("Phantom");
		
		addLoop(phantom).show("Add-Schleife");
		
		addKernel(phantom).show("Add-Kernel");
		
		Grid2D sinogram = Exercise2.radonTrans(phantom, numProj, detSize, detSpacing, angularRange, angularInc, samplingRate);
		sinogram.show("Sinogram");
		long start = System.currentTimeMillis();
		backProj(sinogram).show("Back-Projection");
		long end = System.currentTimeMillis();
		System.out.println("Zeit f�r BackProj: " + (end - start) + "ms --> Vergleich Exercise 2: 3733ms; Vergleich Exercise 3: 16084ms");
		
	}
}
