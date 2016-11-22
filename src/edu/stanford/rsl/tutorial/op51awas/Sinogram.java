package edu.stanford.rsl.tutorial.op51awas;

import imagescience.transform.Transform;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class Sinogram {
	
	private int numOfProjections;
	private double angularIncrement;
	private int detectorSize;
	private double detectorSpacing;
	
	private Box box;
	
	private Grid2D sinogram;
	
	public static void main(String[] args) {
		Sinogram sinogram = new Sinogram(250, 500, 1.0);
		Grid2D phantom = new SimplePhantom(500, 500, new double[]{0.5, 0.5});
		phantom.show("phantom");
		sinogram.computeSinogram(phantom, 1.0);
		sinogram.getSinogram().show("sinogram");
	}
	
	public Sinogram(int numOfProjections, int detectorSize, double detectorSpacing) {
		this.numOfProjections = numOfProjections;
		this.angularIncrement = Math.PI/numOfProjections;
		this.detectorSize = detectorSize;
		this.detectorSpacing = detectorSpacing;
		
		this.sinogram = new Grid2D(detectorSize, numOfProjections);
		this.sinogram.setSpacing(detectorSpacing, angularIncrement);
		double[] origin = new double[2];
		this.sinogram.setOrigin(-(detectorSize-1)*(detectorSpacing/2.0), 0);
	}
	
	public Grid2D computeSinogram(Grid2D image, double samplingRate) {
		
		double width = image.getSize()[0] * image.getSpacing()[0];
		double height = image.getSize()[1] * image.getSpacing()[1];
	
        Translation trans = new Translation(image.getOrigin()[0], image.getOrigin()[1], -1);
		
		box = new Box(width, height, 2);
		box.applyTransform(trans);
		
		double angle = 0.0;
		
		for (int i = 0; i < numOfProjections; i++) {
			angle = i*angularIncrement;
			
			for (int j = 0; j < detectorSize; j++) {
				
				if (j == 450) {
					System.out.println("");
				}
				
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
	
	public Grid2D getSinogram() {
		return this.sinogram;
	}
}