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
		FanBeamRecon rec = new FanBeamRecon(250, 250, 1.0, 1200, 750);
		
		Grid2D sino = rec.computeFanogram(phant, 1.0);
		
		new ImageJ();
		phant.show();
		sino.show();
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
		double curAngle = 0.0;
		
		for (int i = 0; i < numOfProjections; i++) {
			curAngle = i*angularIncrement;
			
			// determine source position in world coordinates
			double xDir = Math.cos(curAngle);
			double yDir = Math.sin(curAngle);
			PointND sourcePos = new PointND(sourceIsoDist*xDir, sourceIsoDist*yDir, 0);
			
			// determine detector middle, thanks to symmetry just the negative of source position
			PointND detectorMiddle = new PointND((-sourceDetDist+sourceIsoDist)*xDir, (-sourceDetDist+sourceIsoDist)*yDir, 0);
			
			if (curAngle > Math.PI/2-angularIncrement && curAngle < Math.PI/2+angularIncrement) {
				System.out.println(sourcePos);
			}
			SimpleVector detectorMid = detectorMiddle.getAbstractVector().clone();
			
			// determine normal vector of conjunction between source position and detector middle
//			double normDirX = Math.cos(Math.PI - curAngle);
//			double normDirY = Math.sin(Math.PI - curAngle);
			double normDirX = xDir;
			double normDirY = -yDir;
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
				if (s == 0) {
					System.out.println(curStep.toString());
				}
				PointND curPoint = new PointND(curStep);
				StraightLine line = new StraightLine(sourcePos, curPoint);
				ArrayList<PointND> hits = box.intersect(line);
				
				if (hits.size() == 2) {
					PointND start = hits.get(0);
					PointND end = hits.get(1);
					SimpleVector lineDir = end.getAbstractVector().clone();
					dir.subtract(start.getAbstractVector());
					double length = lineDir.normL2();
					
					dir.divideBy(length / samplingRate);
					
					float sum = 0.0f;
					for (int t = 0; t < (length/samplingRate); t++) {
						PointND cur = new PointND(start);
						cur.getAbstractVector().add(dir.multipliedBy(t));
						double curX = cur.getCoordinates()[0] * image.getSpacing()[0];
						double curY = cur.getCoordinates()[1] * image.getSpacing()[1];
						
						double[] coords = image.physicalToIndex(curX, curY);
						
						sum += InterpolationOperators.interpolateLinear(image, coords[0], coords[1]);
						
					}
					
					sum /= samplingRate;
					sinogram.setAtIndex(s, i, sum);
					
				}
			}
		}

		return sinogram;
	}

	public void setAngularIncrement() {

	}

}
