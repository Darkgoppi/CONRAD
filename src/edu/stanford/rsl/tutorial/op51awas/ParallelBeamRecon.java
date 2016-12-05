package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import ij.ImageJ;

public class ParallelBeamRecon {
	
	public enum FilterType {NONE, RAMLAK, RAMP};
	
	private int numOfProjections;
	private double angularIncrement;
	private int detectorSize;
	private double detectorSpacing;
	
	private Box box;
	
	private Grid2D sinogram;
	
	public static void main(String[] args) {
		new ImageJ();
		
		ParallelBeamRecon recon = new ParallelBeamRecon(250, 500, 1.0);
		Grid2D phantom = new SimplePhantom(500, 500, new double[]{0.5, 0.5});
		phantom.show("phantom");
		
		Grid2D sinogram = recon.computeSinogram(phantom, 1.0);
		sinogram.show("sinogram");
		recon.backProj(sinogram, 500, 500, 1.0, 1.0).show("Backproj");
		
		Grid2D fSino = recon.filterSino(sinogram, FilterType.RAMLAK);
		fSino.show("RamLak");
		recon.backProj(fSino, 500, 500, 0.5, 0.5).show("RamLak");
	}
	
	public ParallelBeamRecon(int numOfProjections, int detectorSize, double detectorSpacing) {
		this.numOfProjections = numOfProjections;
		this.angularIncrement = Math.PI/numOfProjections;
		this.detectorSize = detectorSize;
		this.detectorSpacing = detectorSpacing;
		
//		this.sinogram = new Grid2D(detectorSize, numOfProjections);
//		this.sinogram.setSpacing(detectorSpacing, angularIncrement);
////		double[] origin = new double[2];
//		this.sinogram.setOrigin(-(detectorSize-1.0)*(detectorSpacing/2.0), 0);
	}
	
	/**
	 * generates the sinogram from a given 2D image
	 * @param image - given image
	 * @param samplingRate - sampling rate along the ray through the image
	 * @return sinogram
	 */
	public Grid2D computeSinogram(Grid2D image, double samplingRate) {
		
		// initialize container for sinogram
		Grid2D sinogram = new Grid2D(detectorSize, numOfProjections);
		sinogram.setSpacing(detectorSpacing, angularIncrement);
		sinogram.setOrigin(-(detectorSize-1.0)*(detectorSpacing/2.0), 0);
		
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
	
	/**
	 * 
	 * @param sinogram
	 * @param width
	 * @param height
	 * @param reconSpacingX
	 * @param reconSpacingY
	 * @return
	 */
	public Grid2D backProj(Grid2D sinogram, int width, int height, double reconSpacingX, double reconSpacingY) {
		
		// get required information for reconstruction
		int sinoWidth = sinogram.getWidth();
		double widthSpacing = sinogram.getSpacing()[0];
		int sinoHeight = sinogram.getHeight();
		double heightSpacing = sinogram.getSpacing()[1];
		
		// create new Grid2D to hold the back projection
		Grid2D backProj = new Grid2D(width, height);
		backProj.setSpacing(reconSpacingX, reconSpacingY);
		backProj.setOrigin(-(width-1.0)/2*reconSpacingX, -(height-1.0)/2*reconSpacingY);
		
		// walk over all lines in sinogram is equal to number of projections in sinogram
		for (int i = 0; i < sinoHeight; i++) {
			double realAngle = i * heightSpacing;
			double xValue = Math.cos(realAngle);
			double yValue = Math.sin(realAngle);
			
			SimpleVector normVek = new SimpleVector(yValue, xValue);
			
			for (int x = 0; x < backProj.getWidth(); x++) {
				for (int y = 0; y < backProj.getHeight(); y++) {
					double[] coords = backProj.indexToPhysical(x, y);
					SimpleVector pix = new SimpleVector(coords[0], coords[1]);
					double innerPro = SimpleOperators.multiplyInnerProd(normVek, pix);
					double dist = innerPro + (sinoWidth-1)/2;
					double index = dist / widthSpacing;
					Grid1D sub = sinogram.getSubGrid(i);
					if(sub.getSize()[0] <= index+1 || index < 0){
						continue;
					}
					float intens = InterpolationOperators.interpolateLinear(sub, index);
					backProj.addAtIndex(x, y, intens);
				}
			}
		}
		NumericPointwiseOperators.divideBy(backProj, (float)(sinoHeight / Math.PI));
		return backProj;
	}
	
	/**
	 * 
	 * @param sinogram
	 * @param filterType
	 * @return
	 */
	public Grid2D filterSino(Grid2D sinogram, FilterType filterType) {
		
		// initialize filter grid
		Grid1DComplex filter = new Grid1DComplex(sinogram.getSize()[0]);
		int filterSize = filter.getSize()[0];
		
		// definition of filters
		if (filterType == FilterType.NONE) {
			// no filter to apply
			return sinogram;
		} else if (filterType == FilterType.RAMLAK) {
			/*
			 * =========================================================================== 
			 * generate ramlak filter in spatial domain and convert it to frequency domain
			 * ===========================================================================
			 */
			// apply definition of ramlak filter to Grid1D
			filter.setAtIndex(0, 0.25f);
			
			float factorOdd = -1.0f/(float)Math.pow(Math.PI, 2);
			for (int i = 1; i < filterSize/2; i++) {
				if (i%2 == 1) {
					filter.setAtIndex(i, factorOdd/(float)Math.pow(i,2));
				} else {
					filter.setAtIndex(i, 0.0f);
				}
			}
			
			float tmp;
			for (int i = filterSize/2; i < filterSize; i++) {
				tmp = filterSize - i;
                if((i%2) == 1){
                    filter.setAtIndex(i, factorOdd / (float) Math.pow(tmp, 2));
                }
			}
			
			// convert filter into frequency domain
//			filter.show();
			filter.transformForward();
//			filter.show();
			
		} else {
			// generate ramp filter directly in frequency domain no conversion required
			
			filter.setRealAtIndex(0, 0.0f);
			filter.setImagAtIndex(0, 0.0f);
			
			for (int i = 1; i < filterSize/2; i++) {
				filter.setRealAtIndex(i, i);
				filter.setImagAtIndex(i, 0.0f);
			}
			float tmp;
			for (int i = filterSize/2; i < filterSize; i++) {
				tmp = filterSize - i;
				filter.setRealAtIndex(i, tmp);
				filter.setImagAtIndex(i, 0.0f);
			}
			
		}
		
		Grid2D filteredSino = new Grid2D(sinogram.getWidth(), sinogram.getHeight());
		
		// walk over all lines, therefore get subgrid line per line and apply filter to line
		for (int i = 0; i < sinogram.getHeight(); i++) {
			Grid1DComplex sinof = new Grid1DComplex(sinogram.getSubGrid(i), true);
			sinof.transformForward();
			
			for (int p = 0; p < sinof.getSize()[0]; p++) {
				sinof.multiplyAtIndex((p),  filter.getRealAtIndex(p), filter.getImagAtIndex(p));
		    }
			sinof.transformInverse();
			
			Grid1D ret = new Grid1D(sinogram.getWidth());
	        ret = sinof.getRealSubGrid(0, sinogram.getSize()[0]);
			for(int p = 0; p < ret.getSize()[0]; p++) {
                filteredSino.putPixelValue(p, i, ret.getAtIndex(p));
            }
		}
		
		return filteredSino;
	}
	
	public Grid2D getSinogram() {
		return this.sinogram;
	}
}