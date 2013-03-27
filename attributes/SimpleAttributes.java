package attributes;

import static edu.mines.jtk.util.ArrayMath.cabs;
import static edu.mines.jtk.util.ArrayMath.carg;
import static edu.mines.jtk.util.ArrayMath.mul;
import static edu.mines.jtk.util.ArrayMath.add;
import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;

import static java.lang.Math.PI;


/**
 * 
 * @author Chris
 * Calculates the amplitude spectrum or phase spectrum by using sub-bands(the imaginary and 
 * real parts of a certain number of frequencies that were created by the Morlet transform 
 * of a single trace). This class also has the capability to calculate the the overall amplitude 
 * spectrum from the sub-bands of multiple traces.
 */
public class SimpleAttributes {
	/**
	 * Finds the amplitude spectrum of the sub-bands of a single trace.
	 * @param subbands - First dimension: number of sub-bands. Second dimension: the number of samples in 
	 * each sub-band(twice the number of samples in the orginal trace to account for 
	 * real and imaginary parts).
	 * @param amplitude - The magnitude of the real and imaginary parts of each sub-band. 
	 * amplitude[f][t] = SQRT(subbands[f][2*t]*subbands[f][2*t] + subbands[f][2*t + 1]*subbands[f][2*t +1 ]),
	 * First dimension: number of sub-bands. Second dimension: number of samples in the original
	 * trace. 
	 */
	public static void findAmplitude(float[][] subbands, float[][] amplitude){
		for (int cf = 0; cf < subbands.length; ++cf){
			amplitude[cf] = cabs(subbands[cf]);
		}
	}
	
	/**
	 * Finds the amplitude spectrum of the sub-bands of a multiple trace.
	 * @param subbands - First dimension: number of traces. Second dimension: number of sub-bands.
	 * Third dimension: the number of samples in 
	 * each sub-band(twice the number of samples in the orginal trace to account for 
	 * real and imaginary parts).
	 * @param amplitude - The magnitude of the real and imaginary parts of each sub-band. 
	 * amplitude[f][t] = SQRT(subbands[f][2*t]*subbands[f][2*t] + subbands[f][2*t + 1]*subbands[f][2*t +1 ]),
	 * First dimension: number of sub-bands. Second dimension: number of samples in the original
	 * traces. 
	 */
	public static void findAmplitudeStack(float[][][] subbands, float[][] amplitude){
		int numTraces = subbands.length;
		int numFreqs = subbands[0].length;
		
		//If there is only one trace
		if (subbands.length == 1){
			for (int cf = 0; cf < numFreqs; ++cf){
				amplitude[cf] = cabs(subbands[0][cf]);
			}
		}
		
		else {
			//adds up all the traces' amplitudes together into one 2D amplitude array.
			for (int tr = 0; tr < numTraces; ++tr){	
				for (int cf = 0; cf < numFreqs; ++cf){
					amplitude[cf] = add(cabs(subbands[tr][cf]), amplitude[cf]);
				}
			}
			//Divide summed amplitude array by the number of traces.
			//float maxAll = max(amplitude);
			for (int cf = 0; cf < numFreqs; ++cf){
				amplitude[cf] = div(amplitude[cf], numTraces);
			}
		}
		
		
	}
	
	public static void findAmplitude(float[][][] subbands, float[][][] amplitude){
		int numTraces = subbands.length;
		int numFreqs = subbands[0].length;
		
		//If there is only one trace
		if (subbands.length == 1){
			for (int cf = 0; cf < numFreqs; ++cf){
				amplitude[0][cf] = cabs(subbands[0][cf]);
			}
		}
		
		else {
			//adds up all the traces' amplitudes together into one 2D amplitude array.
			for (int tr = 0; tr < numTraces; ++tr){	
				for (int cf = 0; cf < numFreqs; ++cf){
					amplitude[tr][cf] = cabs(subbands[tr][cf]);
				}
			}
			
		}
		
		
	}
	
	public static void findAmplitude(float[][][][] subbands, float[][][][] amplitude){
		int n3 = subbands.length;
		int n2 = subbands[0].length;
		int n1 = subbands[0][0].length;
		int n = subbands[0][0][0].length;
		
		int an3 = amplitude.length;
		int an2 = amplitude[0].length;
		int an1 = amplitude[0][0].length;
		int an = amplitude[0][0][0].length;
		
		float[][][][] convertAmp = new float[n3][n2][n][n1];
		//If there is only one trace
		if (n3 == 1 && n2 == 1){
			for (int i1 = 0; i1 < n1; ++i1){
				amplitude[0][0][i1] = cabs(subbands[0][0][i1]);
			}
			
		}
		
		else {
			for (int i3 = 0; i3 < n3; ++i3){	
				for (int i2 = 0; i2 <n2; ++i2){
					for (int i1 = 0; i1 <n1; ++i1){
						amplitude[i3][i2][i1] = cabs(subbands[i3][i2][i1]);
					}
				}
			}
			
		}
		
		
		
		
		
		
	}
	
	/**
	 * Finds the phase spectrum of the sub-bands of a single trace.
	 * @param subbands - First dimension: number of sub-bands. Second dimension: the number of samples in 
	 * each sub-band or twice the number of samples in the orginal trace to account for 
	 * real and imaginary parts.
	 * @param phase - The phase of the real and imaginary parts of each sub-band. 
	 * phase[f][t] = arcTan(subbands[f][2*t + 1]/subbands[f][2*t]),
	 * First dimension: number of sub-bands. Second dimension: number of samples in the original
	 * trace. 
	 */
	public static void findPhase(float[][][] subbands, float[][] phase){
		int numTraces = subbands.length;
		int numFreqs = subbands[0].length;
		
		//If there is only one trace
		if (subbands.length == 1){
			for (int cf = 0; cf < numFreqs; ++cf){
				phase[cf] = mul(carg(subbands[0][cf]), (float)(180/PI));
			}
		}
		
		else {
			//adds up all the traces' phases together into one 2D amplitude array.
			for (int tr = 0; tr < numTraces; ++tr){	
				for (int cf = 0; cf < numFreqs; ++cf){
					phase[cf] = mul(carg(subbands[tr][cf]), (float)(180/PI));
				}
			}
			
			//Divide summed amplitude array by the number of traces.
			for (int cf = 0; cf < numFreqs; ++cf){
				phase[cf] = div(phase[cf], numTraces);
			}
		}
		
	}
}
