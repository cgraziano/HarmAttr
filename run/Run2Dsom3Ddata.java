package run;

import static edu.mines.jtk.util.ArrayMath.div;
import edu.mines.jtk.mosaic.*;

import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import somseis.Converter;
import somseis.SOM2;
import attributes.SimpleAttributes;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import grabbers.DataGrabber;

public class Run2Dsom3Ddata {
	public static void main(String[] args){
		//File file = new File("C://Users/Chris/Documents/!gb.seismic.binary");
		//File file = new File("C://Users/Chris/workspace/Harmonic_Attributes/tpst.binary");
		File file = new File("C:/Users/Chris/401_4_600TPData/tpsz.binary");

		//File file = new File("C:/Users/Chris/401_4_600TPData/tpgg.dat");

		
		//File file = new File("C://Users/Chris/workspace/Harmonic_Attributes/sweep1.binary");
		//File file = new File("C://Users/Chris/workspace/Harmonic_Attributes/impulseat150.binary");
		String sFile = file.getPath();
		
	
		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile,161,357,401);
		

		//File Information
		sN3 = 0;
		eN3 = 160;
		
		sN2 = 0;
		eN2 = 356;
		
		sN1 = 0;
		eN1 = 400;
		dt = .004;
		
		n1 = eN1 - sN1 + 1;
		n2 = eN2 - sN2 +1;
		n3 = eN3 - sN3 +1;
		totaltraces = n3*n2;
		
		Sampling sampN3 = new Sampling(n3, 1.0f, sN3);
		Sampling sampN2 = new Sampling(n2, 1.0f, sN2);
		
		//Frequencies to analyze and number of filters
		fmin = 8;//Hz
		fmax = 30.0;//Hz
		numfilters = 2;
		
		//SOM
		int numIter = 100000;
		int numAttributes = numfilters;
		int n2SOM = 4;
		int n1SOM = 4;
		SOM2 som = new SOM2(numIter, n2SOM, n1SOM, numAttributes);
		
		//final float[][][] shortData = DataGrabber.shorten3DData(rawData, sN3, eN3, sN2, eN2, sN1, eN1);
		System.out.println("Trimmed Data Built");
		
		
	
		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);
		
		float[][][][] subbands = new float[n3][n2][freq.getCount()][2*time.getCount()];
		float[][][][] amp = new float[n3][n2][time.getCount()][freq.getCount()];

		gf.apply(rawData, subbands);
		System.out.println("Morlet Complete");
		
		SimpleAttributes.findAmplitude(subbands, amp);
		System.out.println("Amplitude Complete");

		
		
		
		som.train3D(amp);
		 
		
		for (int i=0; i<n2SOM; ++i){
			for (int j=0; j<n1SOM; ++j){
				System.out.println("nodes = "+som.node(i,j));
			}
		}
		
		float[][][] classData = som.classify3DData(amp);
		ColorMap cm = som.getColorMap(0, n1*n2-1);//new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		ImagePanelGroup img1 = new ImagePanelGroup(classData);
  		img1.setColorModel(cm.getColorModel());
  		
  		ImagePanelGroup img2 = new ImagePanelGroup(rawData);
  		//img2.setColorModel(cm.getColorModel());
		SimpleFrame sf1 = new SimpleFrame();
		sf1.addImagePanels(img1);
		sf1.addImagePanels(img2);
		
		PlotPanelPixels3 ppp31 = new PlotPanelPixels3(PlotPanelPixels3.Orientation.X1RIGHT_X2UP, PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, time,sampN2,sampN3,classData);
		ppp31.setLineColor(Color.GREEN);
		ppp31.addColorBar();
		int slipor1 = (int) ((eN3-sN3)/2.0f);
		int sliro1 = (int) ((eN2-sN2)/2.0f);
		int slivp1 = (int) ((eN1-sN1)/2.0f);
		System.out.println("Por = "+slipor1+"Ro = "+sliro1+"Vp = "+slivp1);
		ppp31.setSlices(slivp1, sliro1, slipor1);
		
		ppp31.setColorModel(cm.getColorModel());
		//ppp31.addColorBar("Normalized Count");
		PlotFrame pf1 = new PlotFrame(ppp31);
		pf1.setVisible(true);
		
		
		//plot
		SwingUtilities.invokeLater(new Runnable() {
		      public void run() {
		    	//plotStacked(shortData, amplitudeStack, phase);
		        //plotFrequencies(shortData, frtrtiamp);
		      }
		    });
	}
	public static void plotFrequencies(float[][] orig, float[][][] amp){
		float fontsize=14;
		//cdp = new Sampling(totaltraces, 1.0, startCDP);
		//time = new Sampling(totalsamples, dt, startSample);
		amp[0][amp[0].length-1][amp[0][0].length-1] = max(amp);
		amp[1][amp[0].length-1][amp[0][0].length-1] = max(amp);
		//amp[2][amp[0].length-1][amp[0][0].length-1] = max(amp);

		ColorMap cm = new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		cm.setValueRange(min(amp),max(amp));

		float[] freq = div(CenterFreqs.calcFloatCenterFreqs(fmin, fmax, dt, numfilters),(float) dt);
		
		//MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp, amp[1], cdp, amp[2]);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp, amp[1]);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Morlet Amplitudes Separated by Frequency");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0,  "CDP Number for Original Data");
		mp1.setXLabel(1,  "CDP Number for "+freq[0]+" Hz");
		mp1.setXLabel(2,  "CDP Number for "+freq[1]+" Hz");
		//mp1.setXLabel(3,  "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(fontsize);
		mp1.setColorModel(1, cm.getColorModel());
		mp1.setColorModel(2, cm.getColorModel());
		//mp1.setColorModel(3, cm.getColorModel());
		mp1.setColorBar("Amplitude");
		
		

		
	}
	
	public static void plotStacked(float[][] origdata, float[][] amplitude, float[][] phase){
		
		
		//cdp = new Sampling(totaltraces, 1.0, startCDP);
		//time = new Sampling(totalsamples, dt, startSample);
		freq = MorletTransform.getLogFrequencySampling(dt);
		System.out.println("*****************");
		float fontsize = 14;
		
		IndexColorModel cmPhase = ColorMap.getHue(0.0,1.0);
		IndexColorModel cmAmplitude = ColorMap.getHueBlueToRed();
		
		SinglePlot sp = new SinglePlot(origdata);
		sp.setTitle("Hi");
		
		
		MultiplePlot mp1 = new MultiplePlot(time, cdp, origdata, freq, amplitude);
		mp1.setColorModel(1, cmAmplitude);
		mp1.setNearestInterpolation(1);
		mp1.setColorBar("Amplitude");
		mp1.setTitle("Gabor-Morlet Amplitude Scalogram");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0, "Amplitude");
		mp1.setXLabel(1, "log10[Frequency(Hz)]");
		mp1.createPNG(1000, 5, "Amplitude, Viking-Graben, 5-85 Hz, 1500 samples.png");
		mp1.setFont(fontsize);
		
		MultiplePlot mp2 = new MultiplePlot(time, cdp, origdata, freq, phase);
		mp2.setColorModel(1, cmPhase);
		mp2.setNearestInterpolation(1);
		mp2.setColorBar("Phase (degrees)");
		mp2.setTitle("Gabor-Morlet Phase Scalogram");
		mp2.setYLabel("Time (s)");
		mp2.setXLabel(0, "Amplitude");
		mp2.setXLabel(1, "log10[Frequency(Hz)]");
		mp2.createPNG(1000, 5, "Phase, Viking-Graben, 5-85 Hz, 1500 samples.png");
		mp2.setFont(fontsize);
		
	}
	
	private static int numfilters, totalsamples, totaltraces, n1Trace, n2Trace;
	private static int sN1, eN1,n1, sN2, eN2, n2, sN3, eN3, n3, startSample, endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq; 
	private static float[][][][] conAmp, frtrtiamp;
}
