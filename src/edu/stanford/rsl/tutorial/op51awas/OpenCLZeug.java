package edu.stanford.rsl.tutorial.op51awas;

import java.io.IOException;
import java.nio.FloatBuffer;
import java.util.Arrays;

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

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGridOperators;
import edu.stanford.rsl.conrad.opencl.OpenCLUtil;

public class OpenCLZeug {

	// OpenCL CheatSheet

	public static void main(String[] args) {
//		addGrid2DOpenCL();
		try {
			addGrid2DHostCode();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void addGrid2DOpenCL() {
		SimplePhantom phant = new SimplePhantom(250, 250, new double[] { 0.5,
				0.5 });

		OpenCLGrid2D clGrid1 = new OpenCLGrid2D(phant);
		OpenCLGrid2D clGrid2 = new OpenCLGrid2D(phant);
		OpenCLGridOperators clOp = OpenCLGridOperators.getInstance();

		for (int i = 0; i < 1000000; i++) {
			clOp.addBy(clGrid1, clGrid2);
		}

		new ImageJ();
		clGrid1.show();
	}

	public static void addGrid2DHostCode() throws IOException {
		int[] size = new int[] { 128, 128 };

		OpenCLGrid2D clGrid1 = new OpenCLGrid2D(new Grid2D(size[0], size[1]));
		Arrays.fill(clGrid1.getBuffer(), 1);

		OpenCLGrid2D clGrid2 = new OpenCLGrid2D(new Grid2D(size[0], size[1]));
		Arrays.fill(clGrid2.getBuffer(), 2);

		OpenCLGrid2D clGridRes = new OpenCLGrid2D(new Grid2D(size[0], size[1]));

		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice dev = context.getMaxFlopsDevice();
		CLCommandQueue comQueue = dev.createCommandQueue();

		CLProgram prog = context.createProgram(OpenCLZeug.class.getResourceAsStream("openCLGridAddition.cl")).build();
		CLKernel kernelFunction = prog.createCLKernel("openCLGridAdd");

		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(size.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(new float[] { size[0], size[1] });
		gImgSize.getBuffer().rewind();

		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY,
				ChannelType.FLOAT);

		clGrid1.getDelegate().prepareForDeviceOperation();
		clGrid1.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> clGrid1Tex = context.createImage2d(clGrid1
				.getDelegate().getCLBuffer().getBuffer(), size[0], size[1],
				format, Mem.READ_ONLY);
		clGrid1.getDelegate().release();

		clGrid2.getDelegate().prepareForDeviceOperation();
		clGrid2.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> clGrid2Tex = context.createImage2d(clGrid2
				.getDelegate().getCLBuffer().getBuffer(), size[0], size[1],
				format, Mem.READ_ONLY);
		clGrid2.getDelegate().release();

		clGridRes.getDelegate().prepareForDeviceOperation();

		comQueue.putWriteImage(clGrid1Tex, true)
				.putWriteImage(clGrid2Tex, true)
				.putWriteBuffer(clGridRes.getDelegate().getCLBuffer(), true)
				.putWriteBuffer(gImgSize, true).finish();

		kernelFunction.rewind();
		kernelFunction.putArg(clGrid1Tex).putArg(clGrid2Tex).putArg(clGridRes.getDelegate().getCLBuffer()).putArg(gImgSize);

		int bpBlockSize[] = { 32, 32 };
		int maxWorkGroupSize = dev.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),
						bpBlockSize[0]),
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),
						bpBlockSize[1]) };
		
		int[] globalWorkSize = new int[]{OpenCLUtil.roundUp(realLocalSize[0], size[0]), OpenCLUtil.roundUp(realLocalSize[1], size[1])};
		
		comQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
		clGridRes.getDelegate().notifyDeviceChange();
		
		new ImageJ();
		new Grid2D(clGridRes).show();
	}
	
	public Grid2D backProjOpenCL(Grid2D sinogram, double reconSpacingX, double reconSpacingY) throws IOException {
		int[] size = sinogram.getSize();

		OpenCLGrid2D clSino = new OpenCLGrid2D(sinogram);

		OpenCLGrid2D clGridRes = new OpenCLGrid2D(new Grid2D(size[0], size[0]));
		clGridRes.setSpacing(reconSpacingX, reconSpacingY);
		clGridRes.setOrigin(-(size[0]-1.0)/2*reconSpacingX, -(size[0]-1.0)/2*reconSpacingY);

		CLContext context = OpenCLUtil.getStaticContext();
		CLDevice dev = context.getMaxFlopsDevice();
		CLCommandQueue comQueue = dev.createCommandQueue();

		CLProgram prog = context.createProgram(OpenCLZeug.class.getResourceAsStream("openCLBackProj.cl")).build();
		CLKernel kernelFunction = prog.createCLKernel("openCLBackProj");

		CLBuffer<FloatBuffer> gImgSize = context.createFloatBuffer(size.length, Mem.READ_ONLY);
		gImgSize.getBuffer().put(new float[] { size[0], size[1] });
		gImgSize.getBuffer().rewind();

		CLImageFormat format = new CLImageFormat(ChannelOrder.INTENSITY,
				ChannelType.FLOAT);

		clSino.getDelegate().prepareForDeviceOperation();
		clSino.getDelegate().getCLBuffer().getBuffer().rewind();
		CLImage2d<FloatBuffer> clGrid1Tex = context.createImage2d(clSino
				.getDelegate().getCLBuffer().getBuffer(), size[0], size[1],
				format, Mem.READ_ONLY);
		clSino.getDelegate().release();

		clGridRes.getDelegate().prepareForDeviceOperation();

		comQueue.putWriteImage(clGrid1Tex, true)
				.putWriteBuffer(clGridRes.getDelegate().getCLBuffer(), true)
				.putWriteBuffer(gImgSize, true).finish();

		kernelFunction.rewind();
		kernelFunction.putArg(clGrid1Tex).putArg(clGridRes.getDelegate().getCLBuffer()).putArg(gImgSize);

		int bpBlockSize[] = { 32, 32 };
		int maxWorkGroupSize = dev.getMaxWorkGroupSize();
		int[] realLocalSize = new int[] {
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),
						bpBlockSize[0]),
				Math.min((int) Math.pow(maxWorkGroupSize, 1 / 2.0),
						bpBlockSize[1]) };
		
		int[] globalWorkSize = new int[]{OpenCLUtil.roundUp(realLocalSize[0], size[0]), OpenCLUtil.roundUp(realLocalSize[1], size[1])};
		
		comQueue.put2DRangeKernel(kernelFunction, 0, 0, globalWorkSize[0], globalWorkSize[1], realLocalSize[0], realLocalSize[1]).finish();
		clGridRes.getDelegate().notifyDeviceChange();
		
		return (new Grid2D(clGridRes));
	}
	
	

}
