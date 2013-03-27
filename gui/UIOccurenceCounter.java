/**
 * The MIT License (MIT)
 
Copyright © 2012 Chris Graziano
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the “Software”), to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and 
to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package gui;

import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.awt.ColorMapListener;
import edu.mines.jtk.awt.ColorMapped;


import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.log10;
import static edu.mines.jtk.util.ArrayMath.min;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.div;




import grabbers.DataGrabber;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.IndexColorModel;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;

import attributes.SimpleAttributes;

import somseis.Converter;

public class UIOccurenceCounter {
	private JButton open, save, run;
	private JLabel instructions, inlineNum, crosslineNum, ixsettings;
	private JTextField input, output; 
	private JFrame frame;
	private String file;
	private static String inputfilepath, outputfilepath;

	/**
	 * Constructs the user interface for the geobody counter.
	 */
	public UIOccurenceCounter(String file) {
		this.file = file;
		
		frame = new JFrame("gui");

		Container c = frame.getContentPane();// sets this container as the
												// container for this frame
		GridBagLayout gbl = new GridBagLayout();

		c.setLayout(gbl);// sets the layout manager from this container to be
							// the GridBagLayout

		GridBagConstraints gbc = new GridBagConstraints();// constraints put to
															// the layout
															// manager

		
		open = new JButton("Select Input File");
		gbc.weightx = .5;
		gbc.fill = GridBagConstraints.HORIZONTAL;
		gbc.gridx = 0;
		gbc.gridy = 1;
		gbl.setConstraints(open, gbc);
		c.add(open);
		

		save = new JButton("Select Output File Location");
		gbc.gridx = 0;
		gbc.gridy = 4;
		gbl.setConstraints(save, gbc);
		c.add(save);

		output = new JTextField();
		gbc.gridx = 1;
		gbc.gridy = 4;
		output.setPreferredSize(new Dimension(300, 30));
		gbl.setConstraints(output, gbc);
		c.add(output);
		
		
		
	
		run = new JButton("Run Harmonic Attributes");
		gbc.weightx = 0.0;
		gbc.ipady = 40;
		gbc.gridwidth = 2;
		gbc.gridheight = 2;
		gbc.gridx = 0;
		gbc.gridy = 10;
		gbl.setConstraints(run, gbc);
		c.add(run);

		// buttons, action listerners, labels added to JPanel
		open.addActionListener(new Open());
		save.addActionListener(new Save());
		run.addActionListener(new Run());
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setSize(1200, 800);
		frame.setVisible(true);
	}

	/**
	 * The outcome of pushing on the "open" button.
	 * @author Chris
	 *
	 */
	class Open implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			JFileChooser c = new JFileChooser();// needed to create ease with
												// selecting files
			int optionselected = c.showOpenDialog(frame);// return result of
															// what user did
															// with JFileChooser

			if (optionselected == JFileChooser.APPROVE_OPTION) {// if file is
																// opened, file
																// paths and
																// names stored
																// in variables
																// below
				
				inputfilepath = c.getSelectedFile().getAbsolutePath();// sent to
																		// arrayinputter

				if (inputfilepath.endsWith(".segy")) {// if selected file does
														// end with .xml, save
														// button enabled, if
														// not already enabled
					save.setEnabled(true);

				}

				else {// if selected file does not end with .xml, save button is
						// disabled to get user to change file selected
					input.setText("Need a .segy file format. Please select another file");
					save.setEnabled(false);
				}

			}
			if (optionselected == JFileChooser.CANCEL_OPTION) {//
				input.setText("You pressed cancel");
				output.setText("");
			}
		}
	}
	

	/**
	 * The outcome of pushing on the "Save" button
	 * @author Chris
	 *
	 */
	class Save implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			JFileChooser c = new JFileChooser();// needed to create ease with
												// selecting files
			int optionselected = c.showSaveDialog(frame);// return result of
															// what user did
															// with JFileChooser

			if (optionselected == JFileChooser.APPROVE_OPTION) {// if file is
																// saved, file
																// paths and
																// names stored
																// in variables
																// below

				output.setText(c.getSelectedFile().getName());

				outputfilepath = c.getSelectedFile().getAbsolutePath();// sent
																		// to
																		// arrayinputter

			}
			if (optionselected == JFileChooser.CANCEL_OPTION) {
				output.setText("You pressed cancel");
				input.setText("");
			}
		}
	}
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * The outcome of pushing on the "run" button
	 * @author Chris
	 *
	 */
	class Run implements ActionListener {
		
		public void actionPerformed(ActionEvent e) {
			
			
			
			float[][] originaldata2 = DataGrabber.grab2DDataFromFile(file);
			startCDP = 0;
			endCDP = 2141;
			startSample = 0;
			endSample = 1499;
			final float[][] origdata = DataGrabber.shorten2DData(originaldata2, startCDP, endCDP, startSample, endSample);
			fmin = 8;//Hz
			fmax = 30.0;//Hz
			dt = .004;
			totalsamples = endSample - startSample + 1;
			totaltraces = endCDP - startCDP + 1;
			numfilters = 3;
			Sampling time = new Sampling(totalsamples, .004, 0);
			
			//3D image like Transforms
			/*float[][][] originaldata3 = DataGrabber.convert2DTo3D(originaldata2);
			ImagePanelGroup panel = new ImagePanelGroup(originaldata3);
			panel.setPercentiles(0, 95);
			SimpleFrame frame = new SimpleFrame();
			frame.addImagePanels(panel);
			*/
			//Gabor Transform
			
			MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
			Sampling freq = MorletTransform.getLogFrequencySampling(dt);
			
			float[][] subbands = new float[freq.getCount()][2*time.getCount()];
			final float[][] amplitudeStack = new float[freq.getCount()][time.getCount()];
			final float[][][] amplitude = new float[totaltraces][freq.getCount()][time.getCount()];

			final float[][] phase = new float[freq.getCount()][time.getCount()];
			float[][][] bigsubbands = new float[totaltraces][freq.getCount()][2*time.getCount()];
			
			gf.apply(origdata, bigsubbands);
			System.out.println("Morlet Complete");
			
			
			SimpleAttributes.findAmplitudeStack(bigsubbands, amplitudeStack);
			SimpleAttributes.findAmplitude(bigsubbands, amplitude);

			SimpleAttributes.findPhase(bigsubbands, phase);
			
			frtrtiamp = new float[freq.getCount()][totaltraces][time.getCount()];
			/*Converter.switchIndices(amplitude, frtrtiamp);
			

			float min = min(frtrtiamp);
			float max = max(frtrtiamp);
			
			conAmp = new float[freq.getCount()][totaltraces][time.getCount()];
			Converter con = new Converter(min, max, (float) 1.0, (float) 256.0);
			con.convert(frtrtiamp, conAmp);*/

		
			
			
				
		
		
			 SwingUtilities.invokeLater(new Runnable() {
			      public void run() {
			    	//plotStacked(origdata, amplitudeStack, phase);
			        plotFrequencies(origdata, frtrtiamp);
			      }
			    });
		
			
		}
	}
	
	public static void plotFrequencies(float[][] orig, float[][][] amp){
		float fontsize=14;
		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		amp[0][amp[0].length-1][amp[0][0].length-1] = max(amp);
		amp[1][amp[0].length-1][amp[0][0].length-1] = max(amp);
		amp[2][amp[0].length-1][amp[0][0].length-1] = max(amp);

		ColorMap cm = new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		cm.setValueRange(min(amp),max(amp));
		

		
		
		
		
		
		float[] freq = div(CenterFreqs.calcFloatCenterFreqs(fmin, fmax, dt, numfilters),(float) dt);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp, amp[1], cdp, amp[2]);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Morlet Amplitudes Seperated by Frequency");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0,  "CDP Number for Original Data");
		mp1.setXLabel(1,  "CDP Number for "+freq[0]+" Hz");
		mp1.setXLabel(2,  "CDP Number for "+freq[1]+" Hz");
		mp1.setXLabel(3,  "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(fontsize);
		mp1.setColorModel(1, cm.getColorModel());
		mp1.setColorModel(2, cm.getColorModel());
		mp1.setColorModel(3, cm.getColorModel());
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
	
	public float[][][] getAmp(){
		return conAmp;
	}
	
	public float[][][] getFreqs(){
		return frtrtiamp;
	}
	

	/*public static void main(String args[]) {
		UIOccurenceCounter ui = new UIOccurenceCounter(args[0]);
		
	}*/
	
	private static int numfilters, totalsamples, totaltraces, startCDP, endCDP, startSample, endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq; 
	private static float[][][] conAmp, frtrtiamp;
	
}
