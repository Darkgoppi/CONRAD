package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.numerics.SimpleVector;

public class Sinogram {
	
	private int numOfProjections;
	private double angularIncrement;
	private int detectorSize;
	private double detectorSpacing;
	
	private Box box;
	
	private Grid2D sinogram;
	
	public static void main(String[] args) {
		Sinogram sinogram = new Sinogram(100, 256, 1.0);
		Grid2D phantom = new SimplePhantom(500, 500, new double[]{0.5, 0.5});
		sinogram.computeSinogram(phantom, 1.0);
		sinogram.getSinogram().show("sinogram");
	}
	
	public Sinogram(int numOfProjections, int detectorSize, double detectorSpacing) {
		this.numOfProjections = numOfProjections;
		this.angularIncrement = 180/numOfProjections;
		this.detectorSize = detectorSize;
		this.detectorSpacing = detectorSpacing;
		
		this.sinogram = new Grid2D(detectorSize, numOfProjections);
		this.sinogram.setSpacing(detectorSpacing);
		double[] origin = new double[2];
		origin[0] = -(detectorSize-1)*(detectorSpacing/2);
		origin[1] = -(numOfProjections-1)*(detectorSpacing/2);
		this.sinogram.setOrigin(origin);
	}
	
	public void computeSinogram(Grid2D image, double samplingRate) {
		
		int width = image.getWidth();
		int height = image.getHeight();
		
		box = new Box(width, height, 1);
		
		double[] min = image.indexToPhysical(0, 0);
		double[] max = image.indexToPhysical(width, height);
		
		double angle = 0.0;
		
		double a = 0;
		double b = 0;
		
		
		for (int i = 0; i < numOfProjections; i++) {
			angle = i*angularIncrement * (Math.PI/180);
			
			for (int j = 0; j < detectorSize; j++) {
				
				double[] detectorPosition = sinogram.indexToPhysical(j, i);
				double detectorX = Math.sin(angle)*detectorPosition[0];
				double detectorY = Math.cos(angle)*detectorPosition[1];
				PointND detectorPoint = new PointND(detectorX, detectorY,0);
				
				double dirX = Math.sin(angle);
				double dirY = Math.cos(angle);
				PointND dirPoint = new PointND((-dirX + detectorX), (dirY + detectorY),0);
				
				StraightLine intLine = new StraightLine(detectorPoint, dirPoint);
				ArrayList<PointND> borderPoints = box.intersect(intLine);
				
				if(borderPoints.size() == 2){
					PointND start = borderPoints.get(0);
					SimpleVector dir = start.getAbstractVector().clone();
					dir.subtract(borderPoints.get(1).getAbstractVector());
					double length = dir.normL2();
					
					dir.divideBy(length / samplingRate);
					
					SimpleVector startVec = start.getAbstractVector();
					
					float sum = 0.0f;
					for (int t = 0; t < (length/samplingRate); t++) {
						startVec.add(dir);
						double x = startVec.getElement(0);
						double y = startVec.getElement(1);
						
						double[] coords = image.physicalToIndex(x, y);
						
						sum += InterpolationOperators.interpolateLinear(image, coords[0], coords[1]);
						
					}
					
					sinogram.setAtIndex(j, i, sum);
					
				}
				
			}
		}
	}
	
	public Grid2D getSinogram() {
		return this.sinogram;
	}
}