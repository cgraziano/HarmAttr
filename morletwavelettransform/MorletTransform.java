package morletwavelettransform;

import static java.lang.Math.*;
import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.ArrayMath.copy;
import static edu.mines.jtk.util.ArrayMath.cabs;
import static edu.mines.jtk.util.ArrayMath.carg;
import static edu.mines.jtk.util.ArrayMath.log10;
import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.mul;
import static edu.mines.jtk.util.ArrayMath.sqrt;



public class MorletTransform {
	
	/***
	 * Creates a Morlet Transform with a time sampling, the number of frequencies the signal will be broken up
	 * into(has to be atleast 1 frequency), and the frequency band that is to be analyzed. Note: if the frequency band 
	 * is 0-80 Hz and the number of frequencies is two, the frequencies that will be analyzed will be 0 and 80 Hz.
	 * @param st The time sampling
	 * @param ni The number of filters
	 * @param fmin The minimum frequency
	 * @param fmax The maximum frequency
	 */
	public MorletTransform(Sampling st, int ni, double fmin, double fmax){
		centerfreqs = CenterFreqs.calcDoubleCenterFreqs(fmin, fmax, st.getDelta(), ni);
		this.ni = ni;
		
for (int i=0; i<ni; ++i){
	System.out.println(centerfreqs[i]/st.getDelta() + " Hz ----- "+centerfreqs[i]+ " cycles/sample");
}
		this.n = st.getCount();

		c = new float[ni][n];
		s = new float[ni][n];
		
		sigma = new double[ni];
		rgf = new RecursiveGaussianFilter[ni];
		double w0 = 5.336;
		double sigmaScale = w0/(2.0*PI);
		for (int i = 0; i < ni; ++i){
			sigma[i] = sigmaScale/centerfreqs[i];
			rgf[i] = new RecursiveGaussianFilter(sigma[i]);
		}
		
		
		
		
		double inside = 0;
		for (int i = 0; i < ni; ++i){
			
			inside = 2*PI*centerfreqs[i];
			
			for (int j = 0; j < n; ++j){
				c[i][j] = (float) cos(inside*j);
				s[i][j] = (float) sin(inside*j);	
			}
		}
		
	}
	
	public void apply(float[] x, int i, float[][] subbands){
		float[] xreal = new float[n];
		float[] ximag = new float[n];
		
		float[] yreal = new float[n];
		float[] yimag = new float[n];
		
		float[] scale = new float[n];
		//create the real and imaginary parts of the input trace (step 1)
		for (int j = 0; j < n; ++j){
			scale[j] = (float) sqrt(sqrt(PI)*2.0*sigma[i]);
			xreal[j] = c[i][j]*x[j]*scale[i];
			ximag[j] = -s[i][j]*x[j]*scale[i];
		}
		//convolve gaussian to x real and imaginary to create y real and imaginary (step 2),
		//which gives the Gabor-Morlet wavelet its gaussian envelope.
		rgf[i].apply0(xreal, yreal);
		rgf[i].apply0(ximag, yimag);
		
		
		
		
		//Apply complex exponential that was taken out before convolution(step 3)
		//This is done seperately to increase speed by taking this out of the convolution, step2.
		for (int j = 0; j < n; ++j){
			subbands[i][2*j] = (yreal[j] * c[i][j] - yimag[j] * s[i][j]);
			subbands[i][2*j + 1] = (yimag[j] * c[i][j] + yreal[j] * s[i][j]);
			
			/*subbands[i][2*j] = (yreal[j] * c[i][j] - yimag[j] * s[i][j]);
			subbands[i][2*j + 1] = (yimag[j] * c[i][j] + yreal[j] * s[i][j]);*/
			
		
		}
		
	}
	
	/**
	 * Applies the Morlet transform to the input x array of floats to 
	 * create a complex array of floats.
	 * NOTE: the output y array needs to be double the input x array.
	 * @param x Input array
	 * @param y Output array
	 */
	public void apply(float[] x, float[][] subbands){
	
		for (int i = 0; i < subbands.length; ++i){
			apply(x, i, subbands);
		}
		
	}
	
	/**
	 * Applies the Morlet transform to each array in a 2D array of floats  to 
	 * create a complex array of floats..
	 * NOTE: Each output y array needs to be double the input x array.
	 * @param x Input array
	 * @param y Output array
	 */
	public void apply(float[][] x, float[][][] subbands){
		int n2 = x.length;
		for (int tr = 0; tr < n2; ++tr){
			apply(x[tr], subbands[tr]);
		}
		//reorg(subbands);
		
		
	}
	
	/**
	 * Applies the Morlet transform to each array in a 3D array of floats  to 
	 * create a complex array of floats..
	 * NOTE: Each output y array needs to be double the input x array.
	 * @param x Input array
	 * @param y Output array
	 */
	public void apply(float[][][] x, float[][][][] subbands){
		int n3 = x.length;
		for (int tr = 0; tr < n3; ++tr){
			apply(x[tr], subbands[tr]);
		}
	}
	
	/***
	 * Reorganizes the output of the morlet data so that 
	 * the first index is the frequency, followed by the lengths of
	 * the data. For Example for 2D data, the output would go from
	 * [totaltraces][freq.getCount()][2*time.getCount()] to 
	 * [freq.getCount()][totaltraces][2*time.getCount()]
	 * @return
	 */
	private void reorg(float[][][] orig){
		int n3 = orig.length;
		int n2 = orig[0].length;
		int n1 = orig[0][0].length;
		
		
		float[][][] reO = new float[n3][n2][n1];
		
		for (int i3 = 0; i3<n3; ++i3){
			for (int i2=0; i2<n2; ++i2){
				for (int i1=0; i1<n1; ++i1){
					reO[i2][i3][i1] = orig[i3][i2][i1]; 
				}
			}
		}
		
		orig = copy(reO);
		
	}
	

	
	
	
	/***
	 * Constructs a new sampling of the log of the center frequencies.
	 * The frequencies should be in cycles/sample. 
	 * @param centerfreqs 
	 * @return The uniform sampling of the log frequencies.
	 */
	public static Sampling getLogFrequencySampling(double dt){
		double[] cfHz = div(centerfreqs, dt);
		//get log of centerfreqs
		int nf = cfHz.length;
		double ff = log10(cfHz[0]);
		double lf = log10(cfHz[nf - 1]);
		double df = (lf - ff)/(nf - 1);
		
		
		/*double[] logcf = log10(centerfreqs);
		Sampling cf = new Sampling(logcf);
		int ncf = cf.getCount();
		double df = cf.getDelta();
		
		double ff = cf.getFirst();//1-cf.getFirst()+cf.getFirst();//to start sampling at 1
		*/
		Sampling slogcf = new Sampling(nf, df, ff);
		return slogcf;
		
	}
	

	private double[] sigma, scalefact;
	private int n, ni;
	private double a;
	private float[][] c, s;
	private RecursiveGaussianFilter[] rgf;
	private static double[] centerfreqs;
	
}
