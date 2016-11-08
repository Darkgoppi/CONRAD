package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.geometry.AbstractCurve;
import edu.stanford.rsl.conrad.geometry.AbstractShape;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.transforms.Transform;

public class Sinogram {
	
	private int numOfProjections;
	private double angularIncrement;
	private int detectorSize;
	private double detectorSpacing;
	
	private Box box;
	
	private Grid2D sinogram;
	
	public static void main(String[] args) {
		
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
			angle = i*angularIncrement;
			
			for (int j = 0; j < detectorSize; j++) {
				
				double[] detectorPosition = sinogram.indexToPhysical(j, i);
				double x = Math.sin(angle)*detectorPosition[0];
				double y = Math.cos(angle)*detectorPosition[0];
				
//				AbstractCurve curve = new AbstractCurve()
//				box.getHitsOnBoundingBox(curve)
					
			}
		}
	}
	
	private double getPositionOnDetector(int idx) {
		int center = detectorSize/2;
		
		if (idx < center) {
			return - (center - idx) * detectorSpacing;
		} else {
			return (center - idx) * detectorSpacing;
		}
	}
}
