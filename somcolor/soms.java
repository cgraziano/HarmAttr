package somcolor;
/*<applet
  archive = soms.jar
  code = soms.class
  width=260
  height=400>
  </applet>
*/
  

/*
  soms.java

  Demonstration of Kohonen Self-Organizing Maps

  Copyright (c) 1999, Tom Germano
  
  Worcester Polytechnic Institute course CS563
  http://www.cs.wpi.edu/~matt/cs563
  3/20/99
*/

/*

  The class soms basically sets up and controls the user interface.

  The class SOM_Thread is a thread that updates the view constantly. I
    put in four modes of playing the SOM, so the thread will update 
    accordingly. The only interesting class here is:
      public void run()

*/

import java.net.*;
import java.applet.*;
import java.util.*;
import java.lang.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.swing.*;

import somcolor.fpoint.ipoint;

public class soms extends JApplet {

  //class of weights and samples which performs operations on the SOM
  static Screen screen = new Screen();

  //thread responsible for updating Screen
  static public SOM_thread thread = new SOM_thread();    

  //Java specific user interface material
  static JApplet japplet = new JApplet();

  static JLabel jbar_label = new JLabel("Percentage Complete = 0%",JLabel.LEFT);
  static JProgressBar jbar = new JProgressBar(0,100);

  static JLabel slabel = new JLabel("Number of Iterations = 100",JLabel.LEFT);
  static JScrollBar scrollbar = new JScrollBar(JScrollBar.HORIZONTAL,100,100,0,1100);
  
  //VCR buttons
  static JButton pause;
  static JButton go1;
  static JButton play;
  static JButton ff;

  //Initialization control
  static JButton reset = new JButton("Reset");
  static JComboBox combo = new JComboBox();

  //Radio buttons for controlling the view
  static JRadioButton use_rgb = new JRadioButton("Colored SOM");
  static JRadioButton use_bw  = new JRadioButton("Similarity SOM");
  static ButtonGroup rgroup = new ButtonGroup();

  //Used for laying out the user interface
  static Container c;
  static JPanel p1 = new JPanel();
  static JPanel p2 = new JPanel();
  static JPanel p3 = new JPanel();
  static JPanel p4 = new JPanel();
  static JPanel p5 = new JPanel();
  static JPanel p6 = new JPanel();


    ImageProducer MyImgProd;
    ImageProducer imgprod;
    URL MyURL;
    Toolkit tool;
/*
  Sets up the radio buttons
*/
  public void setup_radios() {

    class RadioListener implements ActionListener {
      public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand()=="use_rgb")
          screen.set_draw_type(1);
        else
          screen.set_draw_type(0);
      }
    }

    RadioListener rlistener = new RadioListener();

    use_rgb.setActionCommand("use_rgb");
    use_bw.setActionCommand("use_bw");
    use_rgb.addActionListener(rlistener);
    use_bw.addActionListener(rlistener);
    use_rgb.setSelected(true);

    rgroup.add(use_rgb);
    rgroup.add(use_bw);
  }


/*
  More user interface setup
*/
  public void setup_applet() {
    japplet.setVisible(true);

    //Load files from the Jar
    
    ImageIcon ipause;
    ImageIcon igo1;
    ImageIcon iplay;
    ImageIcon iff;
        
    MyURL = getClass().getResource("pause.jpg");
    try {
      imgprod = (ImageProducer) MyURL.getContent();
    } catch(Exception ex) {
      System.out.println(ex.getMessage());
    }
    ipause = new ImageIcon(tool.createImage(imgprod));
            
    MyURL = getClass().getResource("go1.jpg");
    try {
      MyImgProd = (ImageProducer) MyURL.getContent();
    } catch(Exception ex) {
      System.out.println(ex.getMessage());
    }
    igo1 = new ImageIcon(tool.createImage(MyImgProd));

    MyURL = getClass().getResource("play.jpg");
    try {
      MyImgProd = (ImageProducer) MyURL.getContent();
    } catch(Exception ex) {
      System.out.println(ex.getMessage());
    }
    iplay = new ImageIcon(tool.createImage(MyImgProd));

    MyURL = getClass().getResource("ff.jpg");
    try {
      MyImgProd = (ImageProducer) MyURL.getContent();
    } catch(Exception ex) {
      System.out.println(ex.getMessage());
    }
    iff = new ImageIcon(tool.createImage(MyImgProd));

    pause = new JButton(ipause);
    go1   = new JButton(igo1);
    play  = new JButton(iplay);
    ff    = new JButton(iff);

    pause.addActionListener(new ActionListener() {public void actionPerformed(ActionEvent e) {thread.set_state(0); }});
    go1.addActionListener(new ActionListener() {public void actionPerformed(ActionEvent e) {thread.set_state(1); }});
    play.addActionListener(new ActionListener() {public void actionPerformed(ActionEvent e) {thread.set_state(2); }});
    ff.addActionListener(new ActionListener() {public void actionPerformed(ActionEvent e) {thread.set_state(3); }});

    pause.setToolTipText("Pause updating");
    go1.setToolTipText("Advance one iteration");
    play.setToolTipText("Continue updating");
    ff.setToolTipText("Continue updating w/o display");

    pause.setBackground(Color.white);
    go1.setBackground(Color.white);
    play.setBackground(Color.white);
    ff.setBackground(Color.white);
    
    reset.addActionListener( new ActionListener() {public void actionPerformed(ActionEvent e) {thread.reset();} });
    reset.setToolTipText("Initialize the weights");

    combo.getAccessibleContext().setAccessibleName("Initialize Style");
    combo.addItem("Random Points");
    combo.addItem("Four Corners");
    combo.addItem("Three Color Fade");

    scrollbar.addAdjustmentListener( new AdjustmentListener() {
      public void adjustmentValueChanged(AdjustmentEvent e) {
        slabel.setText("Number of Iterations = "+scrollbar.getValue());
        thread.set_MAX_ITER(scrollbar.getValue());
      }
    });    

    jbar_label.setForeground(Color.black);
    slabel.setForeground(Color.black);
  }


/*
  Lays out all the user interface in a neat manner
*/
  public void setup_layout()
    {

    //Setup and layout the screen 
    screen.setLayout(new GridLayout(1,1));
    screen.setVisible(true); 
    screen.setAlignmentX(JPanel.CENTER_ALIGNMENT);
    screen.setAlignmentY(JPanel.CENTER_ALIGNMENT);
    screen.setPreferredSize(new Dimension(200+screen.getX(),200+screen.getY()));
    
    //Layout buttons layout
    p1.setLayout(new FlowLayout());
    p1.add(pause); p1.add(go1); p1.add(play); p1.add(ff);

    //Layout the initialization interface
    p2.setLayout(new FlowLayout());
    p2.add(combo); p2.add(reset);

    //Layout the radio buttons
    p3.setLayout(new FlowLayout());
    p3.add(use_rgb); p3.add(use_bw);

    //Layout the two bars
    p5.setLayout(new GridLayout(4,1));
    p5.add(jbar_label); p5.add(jbar); p5.add(slabel); p5.add(scrollbar);

    p6.setLayout(new BoxLayout(p6,BoxLayout.Y_AXIS));
    p6.add(Box.createRigidArea(new Dimension(5,5)));
    p6.add(screen); p6.add(p5); p6.add(p1); p6.add(p2); p6.add(p3);
    p6.add(Box.createRigidArea(new Dimension(5,5)));


    //Add everything to the main panel
    c = getContentPane();
    c.setLayout(new FlowLayout());
    c.add(Box.createRigidArea(new Dimension(5,5)));
    c.add(p6);
    c.add(Box.createRigidArea(new Dimension(5,5)));
    
    
    //Cleanup
    resize(c.getPreferredSize());

    invalidate();
    validate();
    repaint();
    }

  public void init() {
    tool = Toolkit.getDefaultToolkit();

    setup_applet();
    setup_radios();
    setup_layout();
    
    screen.init_Screen();
     
    thread.init();
    thread.start();
  }

}


class SOM_thread extends Thread {

  //0 - pause
  //1 - go1
  //2 - play
  //3 - ff
  int state=-1;
  int MAX_ITER=0;
  float t;
  float T_INC;

  //Random sample retrieved from class Screen
  fpoint r_sample = new fpoint();

  //Location of the best matching unit of r_sample retrieved from
  //  Screen
  ipoint bmu_loc = r_sample.new ipoint();

  
  public void set_MAX_ITER(int M)  {MAX_ITER=M; T_INC=1.0f/(float)(MAX_ITER);}    
  public void set_state(int f)     {state=f;}
  public void reset() {
       
    soms.screen.set_INIT_STYLE(soms.combo.getSelectedIndex());
    soms.screen.init_v_weights();

    //Get rid of this line to keep the same colors when resetting
    soms.screen.init_v_samples();

    t=0.0f;
    T_INC=1.0f/(float)(MAX_ITER);
    soms.jbar.setValue(0);
    soms.screen.paint(soms.screen.getGraphics());
  }


/*
  Updates the class Screen according the to values retrieved from 
    the user interface
*/
  public void run() {
    while(true) {
      if (t<1.0f && state!=0) {
  
        //Get a random sample
        r_sample = soms.screen.get_random();

        //Find its best matching unit
        bmu_loc  = soms.screen.get_bmu(r_sample);

        //Scale the neighbors according to t
        soms.screen.scale_neighbors(bmu_loc,r_sample,t);

        //Increase t to decrease the number of neighbros and
        //  the amount each weight can learn
        t+=T_INC;
      }
    
      //To prevent it from displaying 101%
      if (t>=1.0f) {
        t=1.0f;  
        state=0;
      }
    
      //Stop the process if it's completed or just advanced one frame
      if (t==1.0f || state==1 || state==2) {
        soms.screen.paint(soms.screen.getGraphics());
        if (state==1)
          state=0;
      }
      
      //Update progress bar
      soms.jbar.setValue(Math.round(t*100.0f));
      soms.jbar_label.setText("Percentage Complete = "+Math.round(t*100.0f)+"%");
      soms.japplet.getContentPane().invalidate();

      //Allow other operations
      yield();
    } //while(true) 
  }


/* 
  Performs some simple initializations
*/
  public void init() {
    if (MAX_ITER==0)
      MAX_ITER=100;
    if (state==-1)
      state=0;

    reset();
  }
}