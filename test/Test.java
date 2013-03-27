package test;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.sgl.*;

public class Test {
	public static void main(String[] args){
		int nx = 2142;
		int ny = 100;
		int nt = 1500;
		float dt = (float) .004;
		float maxf = 85;
		float minf = 5;
		float tstart = 0;
		float[][][] testdata = new float[nx][ny][nt];
		
		//sweep build
		testdata = BuildTestData.build3DSweepData(nt, nx, ny, dt, maxf, minf, tstart);
		
		BuildTestData.write3DData(testdata, "sweep2142x100x15003D.binary");
		
		
		//impulse response 
		//testdata = BuildTestData.impulseData(nt, nx, dt, 150);
		//BuildTestData.write2DData(testdata, "sweep1.binary");
		//SimplePlot.asPixels(testdata);
		
		
		
	
		
		System.out.println("test");
	
	}
}
