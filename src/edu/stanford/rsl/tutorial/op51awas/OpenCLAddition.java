package edu.stanford.rsl.tutorial.op51awas;

import ij.ImageJ;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGrid2D;
import edu.stanford.rsl.conrad.data.numeric.opencl.OpenCLGridOperators;

public class OpenCLAddition {
	
	// OpenCL CheatSheet

	public static void main(String[] args) {
		
		SimplePhantom phant = new SimplePhantom(250, 250, new double[]{0.5, 0.5});
		
		OpenCLGrid2D clGrid1 = new OpenCLGrid2D(phant);
		OpenCLGrid2D clGrid2 = new OpenCLGrid2D(phant);
		OpenCLGridOperators clOp = OpenCLGridOperators.getInstance();
		
		for (int i = 0; i < 1000000; i++) {
			clOp.addBy(clGrid1, clGrid2);
		}
		
		new ImageJ();
		clGrid1.show();
	}
	
}
