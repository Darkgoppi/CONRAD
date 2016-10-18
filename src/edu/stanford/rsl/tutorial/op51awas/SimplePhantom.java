package edu.stanford.rsl.tutorial.op51awas;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;

public class SimplePhantom extends Grid2D{
	
	public static void main(String[] args) {
		
		new ImageJ();
		
		double[] spacing = new double[]{0.5, 0.5};
		
		SimplePhantom phant = new SimplePhantom(256, 256, spacing);
		phant.show();
	}

	public SimplePhantom(int width, int height, double[] spacing) {
		// create Grid2D with given width and height by calling constructor of super class
		super(width, height);
		setSpacing(spacing);
		
		// add a quadratic shape to the phantom
		int boxSize = 0;
		if (width < height) {
			boxSize = (int)(width * 0.75);
		} else {
			boxSize = (int)(height * 0.75);
		}
		
		int startPointI = (width - boxSize)/2;
		int startPointJ = (height - boxSize)/2;
		
		for (int i = startPointI; i < boxSize; i++) {
			for (int j = startPointJ; j < boxSize; j++) {
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
		
		// add a triangle to the phantom
		int bottomWidth = (int)(width * 0.25);
		int triaHeight = (int)(height * 0.4);
		int ratio = triaHeight/bottomWidth;

		int idx = 0;
		for (int j = height - 1; j >= height - triaHeight; j--) {
			for (int i = width - 1; i >= width - bottomWidth - idx * ratio; i--) {
				this.setAtIndex(i, j, 0.6f);
			}
			idx++;
		}
		
	}

}
