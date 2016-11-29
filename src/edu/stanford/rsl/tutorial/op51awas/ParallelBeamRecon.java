package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class ParallelBeamRecon {
	
	private int numOfProjections;
	private double angularIncrement;
	private int detectorSize;
//	private double detectorSpacing;
	
	private Box box;
	
	private Grid2D sinogram;
	
	public static void main(String[] args) {
		ParallelBeamRecon recon = new ParallelBeamRecon(250, 500, 1.0);
		Grid2D phantom = new SimplePhantom(500, 500, new double[]{0.5, 0.5});
		phantom.show("phantom");
		recon.computeSinogram(phantom, 1.0);
		recon.getSinogram().show("sinogram");
		
		recon.backProj(recon.getSinogram(), phantom).show("Backproj");
	}
	
	public ParallelBeamRecon(int numOfProjections, int detectorSize, double detectorSpacing) {
		this.numOfProjections = numOfProjections;
		this.angularIncrement = Math.PI/numOfProjections;
		this.detectorSize = detectorSize;
//		this.detectorSpacing = detectorSpacing;
		
		this.sinogram = new Grid2D(detectorSize, numOfProjections);
		this.sinogram.setSpacing(detectorSpacing, angularIncrement);
//		double[] origin = new double[2];
		this.sinogram.setOrigin(-(detectorSize-1)*(detectorSpacing/2.0), 0);
	}
	
	/**
	 * generates the sinogram from a given 2D image
	 * @param image - given image
	 * @param samplingRate - sampling rate along the ray through the image
	 * @return sinogram
	 */
	public Grid2D computeSinogram(Grid2D image, double samplingRate) {
		
		// compute physical size of input image
		double width = image.getSize()[0] * image.getSpacing()[0];
		double height = image.getSize()[1] * image.getSpacing()[1];
	
		// generate translation to shift box into physical origin of input image
        Translation trans = new Translation(image.getOrigin()[0], image.getOrigin()[1], -1);
		
        // box for overlaying image to get hits of a line at the image borders
		box = new Box(width, height, 2);
		box.applyTransform(trans);
		
		// initial angle of detector
		double angle = 0.0;
		
		
		for (int i = 0; i < numOfProjections; i++) {
			angle = i*angularIncrement;
			
			for (int j = 0; j < detectorSize; j++) {
				
				double detectorPosition = sinogram.indexToPhysical(j, i)[0];
				double detectorX = Math.sin(angle)*detectorPosition;
				double detectorY = Math.cos(angle)*detectorPosition;
				PointND detectorPoint = new PointND(detectorX, detectorY,0);
				
				double dirX = Math.cos(angle);
				double dirY = Math.sin(angle);
				PointND dirPoint = new PointND((-dirX + detectorX), (dirY + detectorY),0);
				
				StraightLine intLine = new StraightLine(detectorPoint, dirPoint);
				ArrayList<PointND> borderPoints = box.intersect(intLine);
				
				if(borderPoints.size() == 2){
					PointND start = borderPoints.get(0);
					PointND end = borderPoints.get(1);
					SimpleVector dir = end.getAbstractVector().clone();
					dir.subtract(start.getAbstractVector());
					double length = dir.normL2();
					
					dir.divideBy(length / samplingRate);
					
					float sum = 0.0f;
					for (int t = 0; t < (length/samplingRate); t++) {
						PointND cur = new PointND(start);
						cur.getAbstractVector().add(dir.multipliedBy(t));
						double x = cur.getCoordinates()[0] * image.getSpacing()[0];
						double y = cur.getCoordinates()[1] * image.getSpacing()[1];
						
						double[] coords = image.physicalToIndex(x, y);
						
						sum += InterpolationOperators.interpolateLinear(image, coords[0], coords[1]);
						
					}
					
					sum /= samplingRate;
					sinogram.setAtIndex(j, i, sum);
					
				}
				
			}
		}
		
		return sinogram;
	}
	
	public Grid2D backProj(Grid2D sinogram, Grid2D phantom) {
		
		int sinowidth = sinogram.getWidth();
		double widthSpacing = sinogram.getSpacing()[0];
		int sinoheight = sinogram.getHeight();
		double heightSpacing = sinogram.getSpacing()[1];
		int width = phantom.getWidth();
		int height = phantom.getHeight();
		double[] spacing = phantom.getSpacing();
		
		Grid2D backProj = new Grid2D(width, height);
		backProj.setSpacing(spacing);
		backProj.setOrigin(-(width-1)/2*spacing[0], -(height-1)/2*spacing[1]);
		
		for (int i = 0; i < sinoheight; i++) {
			double realAngle = i * heightSpacing;
			double xValue = Math.cos(realAngle);
			double yValue = Math.sin(realAngle);
			
			SimpleVector normVek = new SimpleVector(yValue, xValue);
			
			for (int x = 0; x < backProj.getSize()[0]; x++) {
				for (int y = 0; y < backProj.getSize()[1]; y++) {
					double[] coords = backProj.indexToPhysical(x, y);
					SimpleVector pix = new SimpleVector(coords[0], coords[1]);
					double innerPro = SimpleOperators.multiplyInnerProd(normVek, pix);
					double dist = innerPro + sinowidth/2;
					double index = dist / widthSpacing;
					Grid1D sub = new Grid1D(sinogram.getSubGrid(i));
					if(sub.getSize()[0] <= index+1 || index < 0){
						continue;
					}
					float intens = InterpolationOperators.interpolateLinear(sub, index);
					backProj.addAtIndex(x, y, intens);
				}
			}
		}
		NumericPointwiseOperators.divideBy(backProj, (float)(sinoheight / Math.PI));
		return backProj;
	}
	
	public Grid2D getSinogram() {
		return this.sinogram;
	}
}