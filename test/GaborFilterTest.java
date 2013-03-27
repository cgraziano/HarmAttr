package test;

import morletwavelettransform.MorletTransform;
import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.ArrayMath.*;

public class GaborFilterTest {
	public static void main(String[] args){
		double fhz = 5; //(hz)
		double dt = .004; //(seconds/sample)
		double f = fhz*dt; //(cycles/sample)
		int n = 1500;
		
		float[] one = new float[n];//constant ones
		/*for (int i = 0; i < n; ++i){
			one[i] = (float) 1.0;
		}*/
		one[1300] = (float)1.0;
		float[] filtered = new float[2*n];
		MorletTransform gf = new MorletTransform(f, dt, n);
		
		gf.apply(one, filtered);
		double[] dfiltered = new double[2*n];
		for (int i = 0; i < 2*n; ++i){
			dfiltered[i] = (double) filtered[i];
		}
		SimplePlot.asSequence(dfiltered);
		SimplePlot.asSequence(carg(dfiltered));
		
		
		
	}
}
