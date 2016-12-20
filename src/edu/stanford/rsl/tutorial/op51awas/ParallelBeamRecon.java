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
	
	public static void main(String[] args) {
		new ImageJ();
		
		ParallelBeamRecon recon = new ParallelBeamRecon(200, 500, 1.0);
		Grid2D phantom = new SimplePhantom(500, 500, new double[]{0.5, 0.5}, 50.0);
		phantom.show("phantom");
		
		Grid2D sinogram = recon.computeSinogram(phantom, 0.05);
		sinogram.show("sinogram");
		recon.backProj(sinogram, 500, 500, 1.0, 1.0).show("Backproj");
		
		Grid2D lakSino = recon.filterSino(sinogram, FilterType.RAMLAK);
		lakSino.show("RamLak");
		recon.backProj(lakSino, 500, 500, 1.0, 1.0).show("RamLak");
		
		Grid2D ramSino = recon.filterSino(sinogram, FilterType.RAMP);
		ramSino.show("Ramp");
		recon.backProj(ramSino, 500, 500, 1.0, 1.0).show("Ramp");
	}
	
	public ParallelBeamRecon() {
		
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
		Box box = new Box(width, height, 2);
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
		int sinoHeight = sinogram.getHeight();
		double heightSpacing = sinogram.getSpacing()[1];
		
		// create new Grid2D to hold the back projection
		Grid2D backProj = new Grid2D(width, height);
		backProj.setSpacing(reconSpacingX, reconSpacingY);
		backProj.setOrigin(-(width-1.0)/2*reconSpacingX, -(height-1.0)/2*reconSpacingY);
		
		// walk over all lines in sinogram which is equal to number of projections in sinogram
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
					double[] dist = sinogram.physicalToIndex(innerPro, 0);
					Grid1D sub = sinogram.getSubGrid(i);
					float intens = InterpolationOperators.interpolateLinear(sub, dist[0]);
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
			/*
			 * ========================================================================
			 * generate ramp filter directly in frequency domain no conversion required
			 * ======================================================================== 
			 */
			filter.setAtIndex(0, 0.0f);
			
			double deltaf = 1.0/(sinogram.getSpacing()[0] * filterSize);
			for (int j = 1; j < filterSize/2; j++) {
				filter.setAtIndex(j, (float)(Math.abs(j * deltaf)));
			}
			for (int j = filterSize/2; j < filterSize; j++) {
				filter.setAtIndex(j,(float)(Math.abs(((filterSize/2 - 1.0) *  deltaf) - ((j - filterSize/2) * deltaf))));
			}
			
//			filter.show();
			
		}
		
		Grid2D filteredSino = new Grid2D(sinogram);
		
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
}