package run;

import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;
import grabbers.DataGrabber;
import gui.UIOccurenceCounter;

import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import somseis.Converter;
import somseis.SOM2;
import somseis.SOM3;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import attributes.SimpleAttributes;

public class Run3Dsom2Ddata {
	public static void main(String[] args){
		File file = new File("C:/inversion_ws/Senior/gb.seismic.binary");
		//File file = new File("C://Users/Chris/workspace/Harmonic_Attributes/sweep1.binary");
		//File file = new File("C://Users/Chris/workspace/Harmonic_Attributes/impulseat150.binary");
		String sFile = file.getPath();


		float[][] rawData = DataGrabber.grab2DDataFromFile(sFile,2142,1500);
		//File Information
		startCDP = 0;
		endCDP = 2141;
		startSample = 0;
		endSample = 1499;
		dt = .004;

		//Frequencies to analyze and number of filters
		fmin = 8;//Hz
		fmax = 30.0;//Hz
		numfilters = 2;

		totalsamples = endSample - startSample + 1;
		totaltraces = endCDP - startCDP + 1;

		final float[][] shortData = DataGrabber.shorten2DData(rawData, startCDP, endCDP, startSample, endSample);


		Sampling time = new Sampling(totalsamples, .004, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);



		float[][] subbands = new float[freq.getCount()][2*time.getCount()];
		final float[][] amplitudeStack = new float[freq.getCount()][time.getCount()];
		final float[][][] amplitude = new float[totaltraces][freq.getCount()][time.getCount()];

		final float[][] phase = new float[freq.getCount()][time.getCount()];
		float[][][] bigsubbands = new float[totaltraces][freq.getCount()][2*time.getCount()];

		gf.apply(shortData, bigsubbands);
		System.out.println("Morlet Complete");


		SimpleAttributes.findAmplitudeStack(bigsubbands, amplitudeStack);
		SimpleAttributes.findAmplitude(bigsubbands, amplitude);

		SimpleAttributes.findPhase(bigsubbands, phase);
		frtrtiamp = new float[freq.getCount()][totaltraces][time.getCount()];

		float[][][] somAmp = new float[totaltraces][time.getCount()][numfilters];
		for (int i=0; i<totaltraces; ++i){
			for (int j=0; j<time.getCount(); ++j){
				for (int k=0; k<numfilters; ++k){
					somAmp[i][j][k] = amplitude[i][k][j];
				}
			}
		}
		//SOM
		int numIter = 100000;
		int numAttributes = numfilters;
		int n1 = 4;
		int n2 = 4;
		int n3 = 4;
		SOM3 som = new SOM3(numIter, n1, n2, n3, numAttributes);


		som.train(somAmp);
		float[][][][] nodes = som.getNodes();
		/*float[][] weightr1 = new float[n1][n2];
		float[][] weightr2 = new float[n1][n2];
		float[][] weightc1 = new float[n2][n1];
		float[][] weightc2 = new float[n2][n1];
		for (int i=0; i<n2; ++i){
			for (int j=0; j<n1; ++j){
				weightc1[i][j] = nodes[i][j][0];
				weightc2[i][j] = nodes[i][j][1];
				weightr1[j][i] = nodes[i][j][0];
				weightr2[j][i] = nodes[i][j][1];
			}
		}

		String iterations = Integer.toString(numIter);
		PointsView pointsr = new PointsView(weightr1,weightr2);
		PointsView pointsc = new PointsView(weightc1,weightc2);

		pointsr.setMarkColor(java.awt.Color.RED);
		pointsr.setMarkStyle(PointsView.Mark.CROSS);
		pointsr.setMarkSize(4);
		pointsc.setMarkColor(java.awt.Color.RED);
		pointsc.setMarkStyle(PointsView.Mark.CROSS);
		pointsc.setMarkSize(4);
		PlotPanel pp = new PlotPanel();
		pp.addTiledView(pointsr);
		pp.addTiledView(pointsc);

		pp.addTitle(iterations);
		pp.setLimits(0, 0, 7, 7);
		PlotFrame pf = new PlotFrame(pp);
		pf.paintToPng(1000, 5, "100000 iterations, 8x8 som, 2 frequencies");
		pf.setVisible(true);*/

		/*for (int i=0; i<n1; ++i){
			for (int j=0; j<n2; ++j){
				for (int k=0; j<n3; ++k){
					System.out.println("nodes = "+som3.node(k, j ,i));//nodes[i][j][2]);
				}
			}
		}*/
		float[][] classData = som.classifyAllData(somAmp);
		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, shortData,cdp, classData);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Morlet Amplitudes Seperated by Frequency");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0,  "CDP Number for Original Data");
		ColorMap cm = new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		cm.setValueRange(min(classData),max(classData));
		//mp1.setXLabel(3,  "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(14);
		mp1.setColorModel(1, cm.getColorModel());
		//mp1.setColorModel(2, cm.getColorModel());
		//mp1.setColorModel(3, cm.getColorModel());
		mp1.setColorBar("Amplitude");

		ColorMap cm = som.getColorMap(0, n1);//new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		//mp1.setXLabel(3,  "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(14);
		mp1.setColorModel(1, cm.getColorModel());
		//mp1.setColorModel(2, cm.getColorModel());
		//mp1.setColorModel(3, cm.getColorModel());
		mp1.setColorBar("Category Number");
		
		MultiplePlot mpT = new MultiplePlot(time, cdp, shortData,cdp, classData,cm, .5f);
		mpT.setNearestInterpolation(1);
		mpT.setTitle("Categorized Data");
		mpT.setYLabel("Time (s)");
		mpT.setXLabel(0,  "CDP Number");

		//mp1.setXLabel(3,  "CDP Number for "+freq[2]+" Hz");
		mpT.setFont(14);
		//mpT.setColorModel(1, cm.getColorModel());
		//mp1.setColorModel(2, cm.getColorModel());
		//mp1.setColorModel(3, cm.getColorModel());
		mpT.setColorBar("Category Number");
		
		Converter.switchIndices(amplitude, frtrtiamp);




		//plot
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				//plotStacked(shortData, amplitudeStack, phase);
				plotFrequencies(shortData, frtrtiamp);
			}
		});
	}
	public static void plotFrequencies(float[][] orig, float[][][] amp){
		float fontsize=14;
		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		amp[0][amp[0].length-1][amp[0][0].length-1] = max(amp);
		amp[1][amp[0].length-1][amp[0][0].length-1] = max(amp);
		//amp[2][amp[0].length-1][amp[0][0].length-1] = max(amp);

		ColorMap cm = new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		cm.setValueRange(min(amp),max(amp));

		float[] freq = div(CenterFreqs.calcFloatCenterFreqs(fmin, fmax, dt, numfilters),(float) dt);

		//MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp, amp[1], cdp, amp[2]);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp, amp[1]);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Morlet Amplitudes Seperated by Frequency");
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


		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
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

	private static int numfilters, totalsamples, totaltraces, startCDP, endCDP, startSample, endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq; 
	private static float[][][] conAmp, frtrtiamp;
}
