package edu.stanford.rsl.tutorial.op51awas;

import edu.stanford.rsl.conrad.data.numeric.Grid2D;
import ij.ImageJ;

public class Phantom extends Grid2D{

	public Phantom(int width, int height, double spacingX, double spacingY) {
		super(width, height);
		double[] spacing = new double[2];
		spacing[0] = spacingX;
		spacing[1] = spacingY;
		this.setSpacing(spacing);
		double[] origin = new double[2];
		origin[0] =  (double)(-(width-1.0) * (spacing[0]/2.0));
		origin[1] =  (double)(-(height-1.0) * (spacing[1]/2.0));	
		this.setOrigin(origin);
		Quadrat q1 = this.new Quadrat((int)(width*0.12), (int)(height*0.06), (int)(width*0.48), (int)(height*0.36));	//Distanz zu origin
		Quadrat q2 = this.new Quadrat((int)(width*0.06), (int)(height*0.1), (int)(width*0.75), (int)(height*0.56));
		q1.zeichnenQ(this);
		q2.zeichnenQ(this);
		Dreieck d1 = this.new Dreieck((int)(width*0.32), (int)(height*0.46), height);
		d1.zeichnenD(this);
		Dreieck d2 = this.new Dreieck((int)(width*0.26), (int)(height*0.52), height);
		d2.zeichnenD(this);
		Dreieck d3 = this.new Dreieck((int)(width*0.38), (int)(height*0.52), height);
		d3.zeichnenD(this);
		Ellipse e1 = this.new Ellipse((int)(width*0.02), (int)(height*0.14), (int)(width*0.5), (int)(height*0.56), Math.PI/4, 0.75f);
		e1.zeichnenE(this);
		Ellipse k1 = this.new Ellipse((int)(width*0.3), (int)(height*0.3), (int)(width*0.5), (int)(height*0.5), 0, 0.5f);
		k1.zeichnenE(this);
		Ellipse e2 = this.new Ellipse((int)(width*0.34), (int)(height*0.4), (int)(width*0.5), (int)(height*0.5), 0, 0.75f);
		e2.zeichnenE(this);
	}
	
	//---------Quadrat-Klasse--------------------------
	
	private class Quadrat{
		
		private int width;
		private int height;
		private int PosX;
		private int PosY;
		
		public Quadrat(int width, int height, int PosX, int PosY){
			this.width = width;
			this.height = height;
			this.PosX = PosX;
			this.PosY = PosY;
		}
		
		public void zeichnenQ(Phantom p){
			for (int i = 0; i < p.getWidth()-1 ; i++) {
				for (int j = 0; j < p.getHeight()-1; j++) {
					
					//Quadrat
					if(Math.pow((double)i - PosX, 2.0) <= Math.pow((double)width/2, 2.0) 
							&& Math.pow((double)j - PosY, 2.0) <= Math.pow((double)height/2, 2.0)){
						if(p.getAtIndex(i, j) != 0){
							//hier ist schon was
						}else{
							p.setAtIndex(i, j, 1);
						}
					}
				}
			}
		}
	}

	//---------Ellipsen-Klasse--------------------------
	
		private class Ellipse{
			
			private int dirx;
			private int diry;
			private int PosX;
			private int PosY;
			private double angle;
			private float color;
			
			public Ellipse(int dirx, int diry, int PosX, int PosY, double angle, float color){ 	
				this.dirx = dirx;
				this.diry = diry;
				this.PosX = PosX;
				this.PosY = PosY;
				this.angle = angle;
				this.color = color;
			}
			
			public void zeichnenE(Phantom p){
				for (int i = 0; i < p.getWidth()-1 ; i++) {
					for (int j = 0; j < p.getHeight()-1; j++) {
						
						//Ellipse
						if(Math.pow((double)(i - PosX), 2.0) / Math.pow((double)dirx, 2.0) + 
								Math.pow((double)(j - PosY), 2.0) / Math.pow((double)diry, 2.0) <= 1){
							if(p.getAtIndex(i, j) != 0){
								//hier ist schon was
							}else{
								int newx = (int)( (i - PosX) * Math.cos(angle) + (j - PosY) * (-Math.sin(angle))) + PosX;
								int newy = (int)( (i - PosX) * Math.sin(angle) + (j - PosY) * Math.cos(angle)) + PosY;
								p.setAtIndex(newx, newy, color);
							}
						}
					}
				}
			}
		}
		
	//---------Kreis-Klasse--------------------------
		
		@SuppressWarnings("unused")
		private class Kreis{
				
			private int radius;
			private int PosX;
			private int PosY;
				
			public Kreis(int radius, int PosX, int PosY){ 	
				this.radius = radius;
				this.PosX = PosX;
				this.PosY = PosY;
			}
			
			public void zeichnenK(Phantom p){
				for (int i = 0; i < p.getWidth()-1 ; i++) {
					for (int j = 0; j < p.getHeight()-1; j++) {
							
						//Kreis
						if(Math.pow((double)radius, 2.0) >= 
								Math.pow((double)(i - PosX), 2.0) + Math.pow((double)(j - PosY), 2.0)){
							if(p.getAtIndex(i, j) != 0){
								//hier ist schon was
							}else{
								p.setAtIndex(i, j, 0.5f);
							}
						}
					}
				}
			}
		}
	
	//---------Dreiecks-Klasse--------------------------
		
			private class Dreieck{	//feste Gr��e und rechtwinklig
						
				private int PosX;
				private int PosY;
				private int height;
						
				public Dreieck(int PosX, int PosY, int height){ 	
					this.PosX = PosX;
					this.PosY = PosY;
					this.height = height;
				}
					
				public void zeichnenD(Phantom p){
					for(int i = 0; i <= (int)(height * 0.04); i++){
						for (int j = 0; j <= i; j++) {
							p.setAtIndex(PosX, PosY, 1);
							p.setAtIndex(PosX - j, PosY, 1);
							p.setAtIndex(PosX + j, PosY, 1);
						}
						PosY = PosY + 1;
					}
				}
			}
		
			
	//------------------------------------------
	//-----------------MAIN---------------------
	//------------------------------------------
	
	public static void main(String args[]){
		
		new ImageJ();
		
		//Bildgr��e
		int width = 250;
		int height = 250;
		
		//Spacing
		double spacingX = 0.5;
		double spacingY = 1.5;

		//Phantom erstellen
		Phantom p1 = new Phantom(width, height, spacingX, spacingY);
		p1.show("Phantom p1");
		System.out.println("Spacing: " + p1.getSpacing()[0] + ", " + p1.getSpacing()[1]);
		System.out.println("Origin: " + p1.getOrigin()[0] + ", " + p1.getOrigin()[1]);
	}
}