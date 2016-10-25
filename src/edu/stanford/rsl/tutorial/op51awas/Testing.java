package edu.stanford.rsl.tutorial.op51awas;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;

public class Testing {
	
	public static void main(String[] args) {
		
		new ImageJ();
		
		SimplePhantom phant1 = new SimplePhantom(256, 256, new double[]{0.5, 0.5});
		SimplePhantom phant2 = new SimplePhantom(256, 256, new double[]{0.5, 0.5});
		
		Grid2D phant3 = (Grid2D)NumericPointwiseOperators.addedBy(phant1, phant2);
		phant3.show();
		
		NumericPointwiseOperators.subtractBy(phant2, phant3);
		phant2.show();
		
		NumericPointwiseOperators.multiplyBy(phant1, phant1);
		phant1.show();
		
		System.out.println(InterpolationOperators.interpolateLinear(phant2, 42, 42.5));
	}
	
}
