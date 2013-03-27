package test;

import morletwavelettransform.CenterFreqs;

public class CenterFrequencyTest {
	public static void main(String[] args){
		float[] cfreq = CenterFreqs.calcCenterFreqs(5, 120, 21);
		
		for (int i = 0; i < cfreq.length; ++i){
			System.out.println(cfreq[i]);
		}
	}
}
