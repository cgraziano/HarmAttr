package morletwavelettransform;

import static java.lang.Math.*;


public class CenterFreqs {
	/**
	 * 
	 * @param freqmin minimum bandwidth set by user (Hertz)
	 * @param freqmax maximum bandwidth set by user (Hertz)
	 * @param dt The sampling interval of the input data (seconds/sample)
	 * @param numfilters the number of filters set by user (default is 21). Keep in 
	 * mind that this number will set the number of filters above the minimum frequency
	 * Example: Filter number (j=0(minimum frequency), 1, 2, 3... 20, 21)
	 * @return a array of floats that has a length equal to the number of filters
	 * and contains the center frequencies that will be used for sub-band calculation.
	 * THESE FREQUENCIES ARE IN CYCLES/SAMPLE, NOT HERTZ ANYMORE!
	 */
	public static float[] calcFloatCenterFreqs(double freqmin, double freqmax, double dt, int numfilters){
		
		double fmin = freqmin*dt;//fmin in (cycles/sample)
		double fmax = freqmax*dt;//fmax in (cycles/sample)
		
		float[] cfreq = new float[numfilters]; //array with center frequencies
												   //the length is the numfilters + 1 because
											       //the first value is the minimum frequency, which
												   //is not a filter.
		
		//calculate mew value for logarithmic sampling (mew = (fmax/fmin)^(1/numfilters))
		//Used in the equation Fcenterj = (mew^j)*Fstart
		double mew = pow(fmax/fmin,(1.0/(double) (numfilters-1.0)));
		for (int i = 0; i < numfilters; ++i){
			cfreq[i] = (float) (pow(mew, i)*fmin);
		}
		
		return cfreq;
	
	}
	
	/**
	 * 
	 * @param freqmin minimum bandwidth set by user (Hertz)
	 * @param freqmax maximum bandwidth set by user (Hertz)
	 * @param dt The sampling interval of the input data (seconds/sample)
	 * @param numfilters the number of filters set by user (default is 21). Keep in 
	 * mind that this number will set the number of filters above the minimum frequency
	 * Example: Filter number (j=0(minimum frequency), 1, 2, 3... 20, 21)
	 * @return a array of floats that has a length equal to the number of filters
	 * and contains the center frequencies that will be used for sub-band calculation.
	 * THESE FREQUENCIES ARE IN CYCLES/SAMPLE, NOT HERTZ ANYMORE!
	 */
	public static double[] calcDoubleCenterFreqs(double freqmin, double freqmax, double dt, int numfilters){
		double fmin = freqmin*dt;//fmin in (cycles/sample)
		double fmax = freqmax*dt;//fmax in (cycles/sample)
		
		double[] cfreq = new double[numfilters]; //array with center frequencies
												   //the length is the numfilters + 1 because
											       //the first value is the minimum frequency, which
												   //is not a filter.
		
		//calculate mew value for logarithmic sampling (mew = (fmax/fmin)^(1/numfilters))
		//Used in the equation Fcenterj = (mew^j)*Fstart
		double mew = pow(fmax/fmin,(1.0/(double) (numfilters-1.0)));
		for (int i = 0; i < numfilters; ++i){
			cfreq[i] = (float) (pow(mew, i)*fmin);
		}
		
		return cfreq;
	
	}
	
	
	
}
