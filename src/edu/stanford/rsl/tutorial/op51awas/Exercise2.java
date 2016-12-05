package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid1D;
import edu.stanford.rsl.conrad.data.numeric.Grid1DComplex;
import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.data.numeric.NumericPointwiseOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleOperators;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import ij.ImageJ;

public class Exercise2 {

	//-----------------------------------------------------------------------------------------
	public static Grid2D radonTrans(SimplePhantom p, int numProj, int detSize, double detSpacing, double angularRange, double angularInc, double samplingRate) {
		
//		double angularSpacing = (angularRange*(360/(2*Math.PI)))/numProj; 
//		int detOrigin = detSize / 2;
		
		Grid2D sino = new Grid2D(new float[detSize * numProj], detSize, numProj);
		sino.setSpacing(detSpacing, angularInc);
		sino.setOrigin(-(detSize-1)*(detSpacing/2), 0);
		Box b = new Box(p.getSize()[0]*p.getSpacing()[0], p.getSize()[1]*p.getSpacing()[1], 1);
		Translation trans = new Translation(-((p.getSize()[0]-1)*p.getSpacing()[0]/2), -((p.getSize()[1]-1)*p.getSpacing()[1]/2), -1);
		b.applyTransform(trans);
		
		double actualAngle = 0;
		
		for (int i = 0; i < numProj; i++) {
			actualAngle = i * angularInc;
			
			for (int s = 0; s < detSize; s++) {
				double dist = sino.indexToPhysical(s, i)[0];
				double xValue = Math.cos(actualAngle)*dist;
				double yValue = Math.sin(actualAngle)*dist;
				PointND p1 = new PointND(xValue, yValue, 0);
				
				//Richtung der Drehung berechnen
				double pointX = Math.sin(actualAngle);
				double pointY = Math.cos(actualAngle);
				PointND p2 = new PointND(-pointX + xValue, yValue + pointY, 0);
				
				StraightLine strLi = new StraightLine(p1,p2);
				ArrayList<PointND> list = b.intersect(strLi);			//PointND als Datentyp, weil wir ja auch Punkte eingeben �ber die StraightLine
				float sum = 0f;
				if(list.size()==2){
					PointND end = list.get(1);
					PointND start = list.get(0);
					SimpleVector dir = end.getAbstractVector().clone();
					dir.subtract(list.get(0).getAbstractVector());
					double length = dir.normL2();
					dir.divideBy(length / samplingRate);				
//					PointND start1 = new PointND(start.getCoordinates()[0]+((p.getSize()[0]*p.getSpacing()[0]/2)),	//keine R�ckrechnung n�tig, da dieser Punkt
//							start.getCoordinates()[0]+((p.getSize()[1]*p.getSpacing()[1]/2)), 1);					//in der Interpolation direkt gemacht wird!
					for (double j = 0.0; j < length / samplingRate; j++) {
//						PointND current = new PointND(start1);
						PointND current = new PointND(start);
						current.getAbstractVector().add(dir.multipliedBy(j));
						double x = current.getCoordinates()[0] / p.getSpacing()[0];
						double y = current.getCoordinates()[1] / p.getSpacing()[1];
						sum += InterpolationOperators.interpolateLinear(p, p.physicalToIndex(x, y)[0], p.physicalToIndex(x, y)[1]);
					}					//weil hier auf die Funktion PhysicalToIndex zugegriffen wird, muss kein start-punkt aufwendig berechnet werden!
				}
				sum /= samplingRate;
				sino.setAtIndex(s, i, sum);
			}
		}
		return sino;
	}
	
	//-----------------------------------------------------------------------------------------
	private static Grid2D backProj(Grid2D sinogram, SimplePhantom phantom) {
		
		int sinowidth = sinogram.getWidth();
		double widthSpacing = sinogram.getSpacing()[0];
		int sinoheight = sinogram.getHeight();
		double heightSpacing = sinogram.getSpacing()[1];
		double[] spacing = phantom.getSpacing();
		
		Grid2D backProj = new Grid2D(phantom.getSize()[0], phantom.getSize()[1]);
		backProj.setSpacing(spacing);
		backProj.setOrigin(-(phantom.getSize()[0]-1.0)/2*spacing[0], -(phantom.getSize()[1]-1.0)/2*spacing[1]);
		
		for (int i = 0; i < sinoheight; i++) {
			double realAngle = i * heightSpacing;
			double xValue = Math.cos(realAngle);
			double yValue = Math.sin(realAngle);
			
			SimpleVector normVek = new SimpleVector(yValue, xValue);
			
			for (int x = 0; x < backProj.getSize()[0]; x++) {
				for (int y = 0; y < backProj.getSize()[1]; y++) {
					double[] coords = backProj.indexToPhysical(x, y);
					SimpleVector pix = new SimpleVector(coords[0], coords[1]);
					double innerPro = SimpleOperators.multiplyInnerProd(normVek, pix);
					double dist = innerPro + (sinowidth-1.0)/2;
					double index = dist / widthSpacing;
					Grid1D sub = new Grid1D(sinogram.getSubGrid(i));
					if(sub.getSize()[0] <= index+1 || index < 0){
						continue;
					}
					float intens = InterpolationOperators.interpolateLinear(sub, index);
					backProj.addAtIndex(x, y, intens);
				}
			}
		}
		NumericPointwiseOperators.divideBy(backProj, (float)(sinoheight / Math.PI));
		return backProj;
	}
	
	//-----------------------------------------------------------------------------------------
	private static Grid2D rampFilter(Grid2D sinogram) {
		//gleich im Frequenz-Bereich definieren
		
		double widthSpacing = sinogram.getSpacing()[0];
		int height = sinogram.getHeight();
		
//		int width = sinogram.getWidth();
//		double heightSpacing = sinogram.getSpacing()[1];
		
		Grid2D sinoRamp = new Grid2D(sinogram);
		
		for (int i = 0; i < height; i++) {
			
			Grid1D sub = sinogram.getSubGrid(i);
//			sub.show();
			Grid1DComplex subCom = new Grid1DComplex(sub);
			int k = subCom.getSize()[0];					//k sollte 2048 gro� sein (Det-Width ist 600)
			double deltaf = 1/(widthSpacing * k);
			for (int j = 0; j < k; j++) {
				subCom.setAtIndex(j, 0.0f);
			}
			
			for (int j = 1; j < k/2; j++) {
				subCom.setAtIndex(j, (float)(Math.abs(j * deltaf)));
			}
			for (int j = k/2; j < k; j++) {
				subCom.setAtIndex(j,(float)(Math.abs(((k/2 - 1) *  deltaf) - ((j - k/2) * deltaf))));		//((k/2 - 1) *  deltaf) - ((j - k/2) * deltaf)
			}	 
//			if(i == 0){
//				subCom.show();
//				i = 0;
//			}
			
			Grid1DComplex filteredSino = new Grid1DComplex(sub);
			filteredSino.transformForward();
			for (int j = 0; j < filteredSino.getSize()[0]; j++) {
				filteredSino.multiplyAtIndex(j, subCom.getRealAtIndex(j), subCom.getImagAtIndex(j));
			}
			filteredSino.transformInverse();
			Grid1D realSino = new Grid1D(sub.getSize()[0]);
			realSino = filteredSino.getRealSubGrid(0, sub.getSize()[0]);
			
			for (int j = 0; j < realSino.getSize()[0]; j++) {
				sinoRamp.setAtIndex(j, i, realSino.getAtIndex(j));
			}
			
		}
		
//		sinogram.show();
		return sinoRamp;
	}
	
	//-----------------------------------------------------------------------------------------
	private static Grid2D ramLakFilter(Grid2D sinogram) {
		//initialisieren in spatial domain
		
		double width = sinogram.getSpacing()[0];
		int height = sinogram.getSize()[1];
		
		Grid2D sinoRamLak = new Grid2D(sinogram);
		
		for (int i = 0; i < height; i++) {
			Grid1D sub = new Grid1D(sinogram.getSubGrid(i));
//			sub.show();
			Grid1DComplex subCom = new Grid1DComplex(sub);
			int k = subCom.getSize()[0];
			//diskrete Ausf�hrung
			for (int j = 0; j < subCom.getSize()[0]; j++) {
				subCom.setAtIndex(j, 0.0f);
			}
//			subCom.show();
			subCom.setAtIndex(0, 0.25f);	//aus VL-Folien
			for (int j = 1; j < k/2; j++) {
				if(j%2 == 1){
					subCom.setAtIndex(j, (float)(-1/(Math.pow(j, 2) * Math.pow(Math.PI, 2))));
				}
			}
			for (int j = k/2; j < k; j++) {
				if(j%2 == 1){
					subCom.setAtIndex(j, (float)(-1/(Math.pow((k - j), 2) * Math.pow(Math.PI, 2))));
				}
			}
//			if(i == 0){
//				subCom.show();
//				i = 0;
//			}
			
			subCom.transformForward();
//			subCom.show();
			Grid1DComplex filteredSino = new Grid1DComplex(sub);
			filteredSino.transformForward();
			for (int j = 0; j < filteredSino.getSize()[0]; j++) {
				filteredSino.multiplyAtIndex(j, subCom.getRealAtIndex(j), subCom.getImagAtIndex(j));
			}
			filteredSino.transformInverse();
			Grid1D realSino = new Grid1D(sub.getSize()[0]);
			realSino = filteredSino.getRealSubGrid(0, sub.getSize()[0]);
			
			for (int j = 0; j < realSino.getSize()[0]; j++) {
				sinoRamLak.setAtIndex(j, i, realSino.getAtIndex(j));
			}
		}
		
		return sinoRamLak;
	}
	
	//-----------------------------------------------------------------------------------------
	//-----------------------------------MAIN--------------------------------------------------
	//-----------------------------------------------------------------------------------------
	public static void main(String[] args) {
		
		ImageJ ij = new ImageJ();
		
		int width = 500;
		int height = 500;
		int numProj = 180;
		int detSize = 600;
		double detSpacing = 2.0;
		double angularRange = Math.PI;
		double angularInc = angularRange/numProj;
		double spacingX = 1.0;
		double spacingY = 1.0;
		double samplingRate = 0.5;
		
		SimplePhantom phantom = new SimplePhantom(width, height, new double[]{1.0, 1.0});
		phantom.show("Phantom");
		System.out.println("Phantom-Origin: " + phantom.getOrigin()[0] + " " + phantom.getOrigin()[1]);
		System.out.println("Phantom-Spacing: " + phantom.getSpacing()[0] + " " + phantom.getSpacing()[1]);
		
		Grid2D sinogram = radonTrans(phantom, numProj, detSize, detSpacing, angularRange, angularInc, samplingRate);
		sinogram.show("Sinogram");
		System.out.println("Sinogram-Origin: " + sinogram.getOrigin()[0] + " " + sinogram.getOrigin()[1]);
		System.out.println("Sinogram-Spacing: " + sinogram.getSpacing()[0] + " " + sinogram.getSpacing()[1]);
		
		Grid2D backproj = backProj(sinogram, phantom);
		backproj.show("Back Projection");
		System.out.println("Backprojection-Origin: " + backproj.getOrigin()[0] + " " + backproj.getOrigin()[1]);
		System.out.println("Backprojection-Spacing: " + backproj.getSpacing()[0] + " " + backproj.getSpacing()[1]);
		
//		Grid2D filtSino1 = rampFilter(sinogram);
//		filtSino1.show("Filtered Sinogram Ramp");
//		Grid2D filtBackProj1 = backProj(filtSino1, phantom);
//		filtBackProj1.show("Filtered Back Projection Ramp");
		
		Grid2D filtSino2 = ramLakFilter(sinogram);
		filtSino2.show("Filtered Sinogram RamLak");
		Grid2D filtBackProj2 = backProj(filtSino2, phantom);
		filtBackProj2.show("Filtered Back Projection RamLak");
		
	}
}
