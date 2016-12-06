package edu.stanford.rsl.tutorial.op51awas;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
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

	}

	public FanBeamRecon(int numOfProjections, int detectorSize, double detectorSpacing, double dSD) {
		this.detectorSize = detectorSize;
		this.detectorSpacing = detectorSpacing;
		this.numOfProjections = numOfProjections;
		this.sourceIsoDist = dSD/2;
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
			double x = Math.cos(curAngle) * sourceIsoDist;
			double y = Math.sin(curAngle) * sourceIsoDist;
			PointND sourcePos = new PointND(x, y);
			
			// determine detector middle, thanks to symmetry just the negative of source position
			PointND detectorMiddle = new PointND(-x, -y);
			
			// determine normal vector of conjunction between source position and detector middle
			double normDirX = Math.cos(Math.PI - curAngle);
			double normDirY = Math.sin(Math.PI - curAngle);
			SimpleVector dir = new SimpleVector(normDirX, normDirY);
						
			for (int s = 0; s < detectorSize; s++) {
				
			}
		}

		return null;
	}

	public void setAngularIncrement() {

	}

}
