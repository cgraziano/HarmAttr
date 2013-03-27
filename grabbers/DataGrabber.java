package grabbers;

import java.io.DataInputStream;

import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.sgl.*;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;


public class DataGrabber {
	public static float[][] grab2DDataFromFile(String fileName, int x, int t){
			int ni = x;
	    	int nj = t;
	    	
	    	float[][] data2 = new float[ni][nj];
	    	
		    try {

		    	ArrayFile af = new ArrayFile(fileName, "r",ByteOrder.LITTLE_ENDIAN, ByteOrder.BIG_ENDIAN);

		    	af.readFloats(data2);
		 
		    	
		    	System.out.println("Data extracted from file.");

		    	af.close();
		    	
		    } catch (IOException e) {
		    	System.out.println("File Cannot be found!");
		    	throw new RuntimeException(e);
		    }	
		  return data2;
	}
	
	public static float[][][] grab3DDataFromFile(String fileName, int n3, int n2, int n1){
		
    	
    	float[][][] data = new float[n3][n2][n1];
    	
	    try {

	    	ArrayFile af = new ArrayFile(fileName, "r",ByteOrder.BIG_ENDIAN, ByteOrder.LITTLE_ENDIAN);
	    	af.readFloats(data);
	 
	    	
	    	System.out.println("Data extracted from file.");

	    	af.close();
	    	
	    } catch (IOException e) {
	    	System.out.println("File Cannot be found!");
	    	throw new RuntimeException(e);
	    }	
	  return data;
}
	
	public static float[][][] convert2DTo3D(float[][] data2){
		int ni = data2.length;
		int nj = data2[0].length;
		
		float[][][] data3 = new float[ni][nj][1];
		
		for (int i=0; i<ni; ++i){
    		for (int j=0; j<nj; ++j){
    			data3[i][j][0] = data2[i][j];
    		}
    	}
		return data3;
	}
	
	
	/**
	 * Shortens the 2D array in the 1st and 2nd dimension
	 * @param originalarray
	 * @param start1D
	 * @param end1D
	 * @param start2D
	 * @param end2D
	 * @return
	 */
	public static float[][] shorten2DData(float[][] originalarray, int start1D, int end1D, int start2D, int end2D){
		
		float[][] shortarray = new float[end1D-start1D + 1][end2D-start2D + 1];
		
		for (int i = start1D, ii = 0; i <= end1D; ++i, ++ii){
			for (int j = start2D, jj = 0; j < end2D; ++j, ++jj){
				
				shortarray[ii][jj] = originalarray[i][j];
			}
		}
		return shortarray;
	}
	
	public static float[][][] shorten3DData(float[][][] originalarray, int sN3, int eN3, int sN2, int eN2, int sN1, int eN1){
		
		float[][][] shortarray = new float[eN3-sN3 + 1][eN2-sN2 + 1][eN1-sN1+1];
		System.out.println(eN3);
		for (int i = sN3, ii = 0; i < eN3; ++i, ++ii){
			for (int j = sN2, jj = 0; j < eN2; ++j, ++jj){
				for (int k=sN1, kk=0; k<eN1;++k,++kk){
					shortarray[ii][jj][kk] = originalarray[i][j][k];
					
				}
			}
		}
		
		return shortarray;
	}
	
	
}
