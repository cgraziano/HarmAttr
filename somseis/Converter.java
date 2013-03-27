package somseis;

import static edu.mines.jtk.util.ArrayMath.max;


public class Converter {
	/**
	 * 
	 * @param minValue1
	 * @param maxValue1
	 * @param minValue2
	 * @param maxValue2
	 */
	public Converter(float minValue1, float maxValue1, float minValue2, float maxValue2){
		this.minValue1 = minValue1;
		this.minValue2 = minValue2;
		this.maxValue1 = maxValue1;
		this.maxValue2 = maxValue2;
	}
	
	/**
	 * Will take value1, which has a range between minValue1 and maxValue1, and
	 * scale it appropriately to be between maxValue2 and minValue2.
	 * @param value1 the value to be scaled.
	 * @param maxValue1
	 * @param minValue1
	 * @param maxValue2
	 * @param minValue2 
	 * @param value2 the scaled value.
	 */
	public float convert(float value1){
		float v1scale = (value1-minValue1)/(maxValue1-minValue1);
		float value2 = v1scale*(maxValue2-minValue2)+minValue2;
		return value2;
	}
	
	public void convert(float[] values1, float[] values2){
		
		for (int i=0; i<values1.length; ++i){
			values2[i] = convert(values1[i]);
		}
	}
	
	public void convert(float[][] values1, float[][] values2){
		
		for (int i=0; i<values1.length; ++i){
			convert(values1[i], values2[i]);
		}
	}
	
	public void convert(float[][][] values1, float[][][] values2){
		
		for (int i=0; i<values1.length; ++i){
			convert(values1[i], values2[i]);
		}
	}
	
	/**
	 * Switch 1st and 2nd indices.
	 * @param original
	 * @param switched
	 */
	public static void switchIndices(float[][][] original, float[][][] switched){
		for (int i=0; i<original.length; ++i){
			for (int j=0; j<original[0].length; ++j){
				for (int k=0; k<original[0][0].length; ++k){
					switched[j][i][k] = original[i][j][k];
				}
			}
		}
	}
	
	private float maxValue1, minValue1, maxValue2, minValue2;//max amp of all input amplitudes
}
