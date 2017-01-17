package edu.stanford.rsl.tutorial.op51awas;

import ij.ImageJ;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import edu.stanford.rsl.tutorial.op51awas.ParallelBeamRecon.FilterType;

public class FanBeamRecon {

	private int detectorSize;
	private double detectorSpacing;
	private int numOfProjections;
	private double angularRange;
	private double angularIncrement;
	private double sourceIsoDist;
	private double sourceDetDist;

	public static void main(String[] args) {
		
		SimplePhantom phant = new SimplePhantom(250, 250, new double[]{0.5, 0.5});
//		SimplePhantom phant = new SimplePhantom(250, 250, new double[]{0.5, 0.5}, 25.0);
		FanBeamRecon rec = new FanBeamRecon(250, 300, 0.7, 1200, 750);
		
		Grid2D fano = rec.computeFanogram(phant, 0.5);
		Grid2D sino = rec.performRebinning(fano);
		
		new ImageJ();
		phant.show("phantom");
		fano.show("fanogram");
		sino.show("sinogram");
		
		ParallelBeamRecon recon = new ParallelBeamRecon();
		
		recon.backProj(sino, 250, 250, 1.0, 1.0).show("#nofilter");
				
		Grid2D sinoRamLak = recon.filterSino(sino, FilterType.RAMLAK);
		recon.backProj(sinoRamLak, 250, 250, 1.0, 1.0).show("RamLak");
		
		Grid2D sinoRamp = recon.filterSino(sino, FilterType.RAMP);
		recon.backProj(sinoRamp, 250, 250, 1.0, 1.0).show("Ramp");
	}

	public FanBeamRecon(int numOfProjections, int detectorSize, double detectorSpacing, double dSD, double dSI) {
		this.detectorSize = detectorSize;
		this.detectorSpacing = detectorSpacing;
		this.numOfProjections = numOfProjections;
		this.sourceIsoDist = dSI;
		this.sourceDetDist = dSD;

		double fanAngle = Math.atan(((double)detectorSize)/(2 * dSD));
		this.angularRange = Math.PI + 2*fanAngle;
		this.angularIncrement = angularRange / numOfProjections;
	}

	public Grid2D computeFanogram(Grid2D image, double samplingRate) {

		// initialize container for sinogram
		Grid2D fanogram = new Grid2D(detectorSize, numOfProjections);
		fanogram.setSpacing(detectorSpacing, angularIncrement);
		fanogram.setOrigin(-(detectorSize-1.0)*(detectorSpacing/2.0), 0);

		// compute physical size of input image
		double width = image.getSize()[0] * image.getSpacing()[0];
		double height = image.getSize()[1] * image.getSpacing()[1];

		// generate translation to shift box into physical origin of input image
		Translation trans = new Translation(image.getOrigin()[0], image.getOrigin()[1], -1);

		// box for overlaying image to get hits of a line at the image borders
		Box box = new Box(width, height, 2);
		box.applyTransform(trans);

		// initial angle of detector
		double curAngle = 0.0;
		
		for (int i = 0; i < numOfProjections; i++) {
			curAngle = i*angularIncrement;
			
			// determine source position in world coordinates
			double xDir = Math.cos(curAngle);
			double yDir = Math.sin(curAngle);
			PointND sourcePos = new PointND(sourceIsoDist*xDir, sourceIsoDist*yDir, 0);
			
			// determine detector middle, thanks to symmetry just the negative of source position
			PointND detectorMiddle = new PointND((-sourceDetDist+sourceIsoDist)*xDir, (-sourceDetDist+sourceIsoDist)*yDir, 0);
			
			SimpleVector detectorMid = detectorMiddle.getAbstractVector().clone();
			
			// determine normal vector of conjunction between source position and detector middle
//			double normDirX = Math.cos(Math.PI - curAngle);
//			double normDirY = Math.sin(Math.PI - curAngle);
			double normDirX = yDir;
			double normDirY = -xDir;
			SimpleVector dir = new SimpleVector(normDirX, normDirY, 0);
			
			double detectorHalf = (detectorSize*detectorSpacing)/2;
			for (int s = 0; s < detectorSize; s++) {
				double step = 0.0;
				if (s < detectorSize/2) {
					step = -detectorHalf + s * detectorSpacing;
				} else {
					step = s * detectorSpacing - detectorHalf;
				}
				
				SimpleVector curStep = dir.multipliedBy(step).clone();
				curStep.add(detectorMid);
				PointND curPoint = new PointND(curStep);
				StraightLine line = new StraightLine(sourcePos, curPoint);
				ArrayList<PointND> hits = box.intersect(line);
				
				if (hits.size() == 2) {
					PointND start = hits.get(0);
					PointND end = hits.get(1);
					SimpleVector lineDir = end.getAbstractVector().clone();
					lineDir.subtract(start.getAbstractVector());
					double length = lineDir.normL2();
					
					lineDir.divideBy(length / samplingRate);
					
					float sum = 0.0f;
					for (int t = 0; t < (length/samplingRate); t++) {
						PointND cur = new PointND(start);
						cur.getAbstractVector().add(lineDir.multipliedBy(t));
						double curX = cur.getCoordinates()[0] * image.getSpacing()[0];
						double curY = cur.getCoordinates()[1] * image.getSpacing()[1];
						
						double[] coords = image.physicalToIndex(curX, curY);
						
						sum += InterpolationOperators.interpolateLinear(image, coords[0], coords[1]);
						
					}
					
					sum /= samplingRate;
					fanogram.setAtIndex(s, i, sum);
					
				}
			}
		}

		return fanogram;
	}
	
	public Grid2D performRebinning(Grid2D fanogram) {
		
		double[] spacing = fanogram.getSpacing();
		int height = 180;
		spacing[1] = Math.PI/height; // 1 .0 fuer Blume
		Grid2D sinogram = new Grid2D(fanogram.getWidth(), height);
		sinogram.setSpacing(spacing);
		sinogram.setOrigin(fanogram.getOrigin());
		
		for (int s = 0; s < sinogram.getWidth(); s++) {
			for (int theta = 0; theta < sinogram.getHeight(); theta++) {
				double[] physCoords = sinogram.indexToPhysical(s, theta);
				
				double gamma = Math.asin((physCoords[0]/sourceIsoDist));
				double beta = physCoords[1] - gamma;
				
				if (beta < 0) {
					gamma = -gamma;
					beta += Math.PI - 2 * gamma;
				}
				
				double t = Math.tan(gamma) * sourceDetDist;
				
				double[] fanoIdxCoords = fanogram.physicalToIndex(t, beta);
				
				
				
				// TODO needs some more tweeking: I know there is a problem with height indices
//				if (fanoIdxCoords[1] < 0) {
//					fanoIdxCoords = fanogram.physicalToIndex(-t, beta);
//					fanoIdxCoords[1] += (fanogram.getHeight()-1);
//				} else if (fanoIdxCoords[1] > (fanogram.getHeight())-1) {
//					fanoIdxCoords = fanogram.physicalToIndex(-t, beta);
//					fanoIdxCoords[1] -= (fanogram.getHeight()-1);
//				}
				
				sinogram.setAtIndex(s, theta, InterpolationOperators.interpolateLinear(fanogram, fanoIdxCoords[0], fanoIdxCoords[1]));
			}
		}
		
		return sinogram;
	}

	public void setAngularIncrement() {

	}

}
