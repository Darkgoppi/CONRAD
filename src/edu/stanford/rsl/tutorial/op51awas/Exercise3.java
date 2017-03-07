package edu.stanford.rsl.tutorial.op51awas;

import java.util.ArrayList;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import edu.stanford.rsl.conrad.data.numeric.InterpolationOperators;
import edu.stanford.rsl.conrad.geometry.shapes.simple.Box;
import edu.stanford.rsl.conrad.geometry.shapes.simple.PointND;
import edu.stanford.rsl.conrad.geometry.shapes.simple.StraightLine;
import edu.stanford.rsl.conrad.geometry.transforms.Translation;
import edu.stanford.rsl.conrad.numerics.SimpleVector;
import ij.ImageJ;

public class Exercise3 {
	
	//-----------------------------------------------------------------------------------------
	public static Grid2D fanogram(Phantom phantom, double dSD, double dSI, int detSize, double detSpacing, int numOfProj, double fanAngle, double angularRange, double angularInc){
		
		double samplingRate= 1.0;
		
		Grid2D fano = new Grid2D(detSize, numOfProj);
		fano.setSpacing(detSpacing, angularInc);
		fano.setOrigin((-((double)detSize-1.0)*(detSpacing/2)), 0);
		
		Box b = new Box(((double)phantom.getSize()[0])*phantom.getSpacing()[0], ((double)phantom.getSize()[1])*phantom.getSpacing()[1], 2);
		Translation trans = new Translation(-(((double)phantom.getSize()[0]-1.0)*phantom.getSpacing()[0]/2), -(((double)phantom.getSize()[1]-1.0)*phantom.getSpacing()[1]/2), -1);
		b.applyTransform(trans);
		
		double actualAngle = 0;
		for (int i = 0; i < numOfProj; i++) {
			actualAngle = i * angularInc;
			
			// Source Point berechnen (f�r einen Winkel immer gleich), bereits in Weltkoords
			double x = Math.cos(actualAngle);
			double y = Math.sin(actualAngle);
			PointND source = new PointND(dSI*x, dSI*y, 0);
			
			// DetectorMiddle
			PointND detectorMiddle = new PointND((-dSD+dSI)*x, (-dSD+dSI)*y, 0);
			SimpleVector detMiddle = detectorMiddle.getAbstractVector().clone();
			
			double normDirX = y;
			double normDirY = -x;
			SimpleVector dir = new SimpleVector(normDirX, normDirY, 0);
			
			double detHalf = (((double)detSize-1.0)*detSpacing)/2; 		//!!!!!!hier vielleicht die +1.0??!!!!!!!!!!!!
			
			for (int j = 0; j < detSize; j++) {
				double step = 0.0;
				if(j < detSize/2){
					step = -detHalf + j*detSpacing;
				}else{
					step = j*detSpacing - detHalf;
				}
				SimpleVector curStep = dir.multipliedBy(step).clone();
				curStep.add(detMiddle);
				PointND curPoint = new PointND(curStep);
				StraightLine line = new StraightLine(source, curPoint);
				ArrayList<PointND> hits = b.intersect(line);
				
				float sum = 0f;
				if(hits.size() == 2){
					PointND end = hits.get(1);
					PointND start = hits.get(0);
					SimpleVector dir2 = end.getAbstractVector().clone();
					dir2.subtract(hits.get(0).getAbstractVector());
					double length = dir2.normL2();
					dir2.divideBy(length / samplingRate);
					
					for (double s = 0.0; s < length / samplingRate; s++) {
						PointND current = new PointND(start);
						current.getAbstractVector().add(dir2.multipliedBy(s));
						double realX = current.getCoordinates()[0] * phantom.getSpacing()[0];
						double realY = current.getCoordinates()[1] * phantom.getSpacing()[1];
						sum += InterpolationOperators.interpolateLinear(phantom, phantom.physicalToIndex(realX, realY)[0], phantom.physicalToIndex(realX, realY)[1]);
					}
				}
				sum /= samplingRate;
				fano.setAtIndex(j, i, sum);
			}
		}
		return fano;
	}
	
	
	//-----------------------------------------------------------------------------------------
	public static Grid2D rebinning(Grid2D fanogram, double dsi, double dsd){
		
		Grid2D sinogram = new Grid2D(fanogram.getSize()[0], (int)(Math.PI/fanogram.getSpacing()[1]));		//!!!!!!hier abh�ngige Gr��e!!!!!!!!!!!!
		sinogram.setSpacing(fanogram.getSpacing()[0], fanogram.getSpacing()[1]);
		sinogram.setOrigin(fanogram.getOrigin());
		
		//for-Schleifen gedreht
		for (int theta = 0; theta < sinogram.getHeight(); theta++) {
			for (int s = 0; s < sinogram.getWidth(); s++) {
				double[] physCoords = sinogram.indexToPhysical(s, theta);
				
				double gamma = Math.asin(physCoords[0]/dsi);
				double beta = physCoords[1] - gamma;
				if(beta < 0){
					gamma = - gamma;
					beta += Math.PI - (2 * gamma);
				}
				double t = Math.tan(gamma) * dsd;
				
				double[] fanoPhysCoords = fanogram.physicalToIndex(t, beta);
				float intval = InterpolationOperators.interpolateLinear(fanogram, fanoPhysCoords[0], fanoPhysCoords[1]);
				sinogram.setAtIndex(s, theta, intval);
			}
		}
		return sinogram;
	}
	
	//-----------------------------------------------------------------------------------------
	//-----------------------------------MAIN--------------------------------------------------
	//-----------------------------------------------------------------------------------------
	public static void main(String[] args){
		
		new ImageJ();
		
		//--------Variablen--------
		//Phantom
		int width = 200;
		int height = 200;
		double spacingX = 1.0;		//0.5
		double spacingY = 1.0;		//0.5
		
		//fanogram
		int detSize = 400;				//gr��ere Det-Size ben�tigt auch gr��ere numOfProj (z.B. detSize = 400 --> numOfProj = 250)
		double detSpacing = 0.75;		//0.9; bei detSize = 400 und numOfProj = 250 geht auch 0.75
		int numOfProj = 250;
		//int numOfDetEl = 200;
		double dSD = 1200.0;
		double dSI = 750.0;
		double fanAngle = 2.0*(Math.atan((double)(detSize/2/dSD)));
		double angularRange = Math.PI + fanAngle;
		double angularInc = angularRange / numOfProj;
		
		Phantom phantom = new Phantom(width, height, spacingX, spacingY);
		phantom.show("Phantom");
		System.out.println("Phantom-Origin: " + phantom.getOrigin()[0] + " " + phantom.getOrigin()[1]);
		System.out.println("Phantom-Spacing: " + phantom.getSpacing()[0] + " " + phantom.getSpacing()[1]);
		
		Grid2D fanogram = fanogram(phantom, dSD, dSI, detSize, detSpacing, numOfProj, fanAngle, angularRange, angularInc);
		fanogram.show("Fanogram");
		System.out.println("Fanogram-Origin: " + fanogram.getOrigin()[0] + " " + fanogram.getOrigin()[1]);
		System.out.println("Fanogram-Spacing: " + fanogram.getSpacing()[0] + " " + fanogram.getSpacing()[1]);
		
		Grid2D sinogram = rebinning(fanogram, dSI, dSD);
		sinogram.show("Sinogram");
		System.out.println("Sinogram-Origin: " + sinogram.getOrigin()[0] + " " + sinogram.getOrigin()[1]);
		System.out.println("Sinogram-Spacing: " + sinogram.getSpacing()[0] + " " + sinogram.getSpacing()[1]);
		
		long start = System.currentTimeMillis();
		Grid2D backproj = Exercise2.backProj(sinogram);
		long end = System.currentTimeMillis();
		System.out.println("Zeit f�r BackProj: " + (end - start) + "ms");
		backproj.show("Backprojection");
		
		Grid2D filtered = Exercise2.ramLakFilter(sinogram);
		Grid2D backproj2 = Exercise2.backProj(filtered);
		backproj2.show("Filtered Backprojection");
	}
}
