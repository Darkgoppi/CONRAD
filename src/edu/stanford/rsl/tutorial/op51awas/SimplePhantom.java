package edu.stanford.rsl.tutorial.op51awas;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;

public class SimplePhantom extends Grid2D{
	
	public static void main(String[] args) {
		
		new ImageJ();
		
		double[] spacing = new double[]{0.5, 0.5};
		
		SimplePhantom phant = new SimplePhantom(500, 500, spacing);
		phant.show();
		
		double[] origin = phant.getOrigin();
		for (int i = 0; i < origin.length; i++) {
			System.out.println(origin[i]);
		}
	}

	public SimplePhantom(int width, int height, double[] spacing) {
		// create Grid2D with given width and height by calling constructor of super class
		super(width, height);
		setSpacing(spacing);
		
		double[] origin = new double[2];
		origin[0] = -(width-1)*(spacing[0]/2);
		origin[1] = -(height-1)*(spacing[1]/2);
		
		setOrigin(origin);
		
		// add a quadratic shape to the phantom
		int boxSize = 0;
//		int rectangleWidth = (int)(width*0.3);
		if (width < height) {
			boxSize = (int)(width * 0.66);
		} else {
			boxSize = (int)(height * 0.66);
		}
		
		int startPointI = (width - boxSize)/3;
		int startPointJ = (height - boxSize)/2;
		
		for (int i = startPointI; i < startPointI+boxSize; i++) {
			for (int j = startPointJ; j < startPointJ+boxSize; j++) {
				this.setAtIndex(i, j, 0.4f);
			}
		}
		
		// add a circle shape to the phantom
		double radius = 0.0;
		if (width < height) {
			radius = width * 0.25;
		} else {
			radius = height * 0.25;
		}
		int[] center = new int[2];
		center[0] = width / 2;
		center[1] = height / 2;
		
		for (int i = center[0]-(int)radius; i <= center[0]+(int)radius; i++) {
			for (int j = center[1]-(int)radius; j <= center[1]+(int)radius; j++) {
				double tmpRadius = Math.sqrt(Math.pow((i-center[0]), 2) + Math.pow((j-center[1]), 2));
				if (tmpRadius <= radius) {
					this.setAtIndex(i, j, 0.8f);
				}
			}
		}
		

		
//		int recStartI = (height - boxSize)/4;
//		int recStartJ = (width - rectangleWidth)/2;
//		for (int i = recStartI; i < boxSize; i++) {
//			for (int j = recStartJ; j < rectangleWidth; j++) {
//				this.setAtIndex(i, j, 0.6f);
//			}
//		}
		
		
		// add a triangle to the phantom
//		int triaStart = height - boxSize;
//		int triaWStart = width - boxSize;
		int startPoint = (int)(boxSize*0.2);
		for (int i = startPoint; i < boxSize/2; i++) {
			for (int j = startPoint; j <= i; j++) {
				this.setAtIndex(i, j, 0.6f);
			}
		}
		
	}

}
