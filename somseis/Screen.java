package somseis;
/*
  Screen.java
  
  Demonstration of Kohonen Self-Organizing Maps

  Copyright (c) 1999, Tom Germano
  
  Worcester Polytechnic Institute course CS563
  http://www.cs.wpi.edu/~matt/cs563
  3/20/99
*/

/*  
  This class does the actual operations on the SOM. All these classes
  are called by the SOM_Thread class in soms.java. The interesting
  procedures here are: 
    public void init_v_map() 
    public ipoint get_bmu(fpoint input) 
    public void scale_neighbors(ipoint loc, fpoint factual, float t2) 
    public void update_bw()
*/


import gui.UIOccurenceCounter;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import somseis.fpoint.ipoint;

import java.applet.*;
import java.util.*;
import java.io.File;
import java.lang.*;
import java.net.URI;

 
public class Screen extends JPanel {
	
	 private float[][][] amp;
	
  //Weight values of the SOM
  fpoint v_weights[][]   = new fpoint [100][100];

  //Randomly selected colors
  fpoint v_samples[] = new fpoint [100];

  //Dimensions of the Map (multiplied by 4 in paint(Graphics)) 
  int HEIGHT     = 100;
  int WIDTH      = 100;

  //Number of samples generated
  int MAX_PTS    = 15;
  int MAX_ACTUAL = 100;

  //Initial radius of influence (should start large)
  int RADIUS     = 120;

  //How many neighbors are used in the calculation of the
  //  similarity map, the more the higher the quality
  int SIMWEIGHT  = 10;

  //Whether to use color or black and white palette
  int RGB        = 1;

  //Way to initialize the weight values
  int INIT_STYLE = -1;

  java.util.Random r = new java.util.Random();

  //Temporary buffer for the weights
  public Color vrscr[][] = new Color [100][100];

  //Red, green, blue palette
  Color rgb_pal[] = new Color [256];

  //Black and white palette
  Color bw_pal[]  = new Color [256];

  //Is filled with either rgb or bw depending on the variable RGB
  Color cur_pal[] = new Color [256];

  //Lookup table for finding the right Color number for the given
  //  [r][g][b] components
  int rgb_table[][][] = new int[6][6][6];


  public void set_MAX_PTS(int m)    {MAX_PTS=m;}
  public void set_RADIUS(int m)     {RADIUS=m;}
  public void set_INIT_STYLE(int m) {INIT_STYLE=m;}


/*
  Calculates the Euclidean distance
*/
  public float get_dist(fpoint imap, fpoint iactual) {
    fpoint d= new fpoint();
    d = imap.sub(iactual);
    d.set(d.x*d.x,d.y*d.y,d.z*d.z);
    
    return ((float)Math.sqrt(d.x+d.y+d.z));
  }


/*
  Calculates the Euclidean distance without the square root, saves time
*/
  public float fget_dist(fpoint imap, fpoint iactual) {
    fpoint d= new fpoint();
    d = imap.sub(iactual);
    d.set(d.x*d.x,d.y*d.y,d.z*d.z);
    
    return (d.x+d.y+d.z);
  }


/*
  Sets the current palette to rgb
*/
  public void set_rgb() {
    for (int loop=0; loop<256; loop++) 
      cur_pal[loop] = rgb_pal[loop];
  }

  
/*
  Sets the current palette to black and white
*/  
  public void set_bw() {
    for (int loop=0; loop<256; loop++) 
      cur_pal[loop] = bw_pal[loop];
  }

  
/*
  Initializes variables, is only called at the beginning of the applet
*/  
  public void init_Screen() {
	  
	 
	
    for (int loop=0; loop<256; loop++) {
      bw_pal[loop] = new Color(loop,loop,loop);
      cur_pal[loop] = new Color(0,0,0);
    }
    
    int col_loop=0;

    //Generates the rgb lookup table and palette
    for (int r=0; r<6; r++)
      for (int g=0; g<6; g++)
        for (int b=0; b<6; b++) {
        rgb_table[r][g][b]= col_loop;
        rgb_pal[r*36+g*6+b] = new Color (r*(256/5),g*(256/5),b*(256/5));
        col_loop++;
        }

    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
        v_weights[loop2][loop] = new fpoint();
        vrscr[loop2][loop] = new Color(0);
        }

    for (int loop=0; loop<MAX_ACTUAL; loop++)
      v_samples[loop] = new fpoint();

    set_rgb();
    INIT_STYLE=0;

    init_v_samples();
    init_v_weights();
    amp = ui.getAmp();
   }


/*
  Changes the palette
*/
  public void set_draw_type(int r) {
    if (r!=RGB) {
      RGB = r;
      if (RGB==1)
        set_rgb();
      else
        set_bw();
      paint(getGraphics());
    }
  }


/*
  Initializes the weight vectors in three different ways defined by
    INIT_STYLE.  
    INIT_STYLE==0 The wieghts all have random values for r,g,b, 
    INIT_STYLE==1 Red, green, blue, and black are all at corners of the screen and 
      their intensity falls off as the pixel moves away from the corner  
    INIT_STYLE==2 Red, green, and blue are all equidistant from the center
      and their intensities fall off with the distance. Probably the best map 
      to start with
*/
  public void init_v_weights() {

    //Random initiation
    if (INIT_STYLE==0)
    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
       v_weights[loop2][loop].x = (float)(r.nextInt(500))/100.0f;
       v_weights[loop2][loop].y = (float)(r.nextInt(500))/100.0f;
       v_weights[loop2][loop].z = (float)(r.nextInt(500))/100.0f;
       }

    //Have red,green,blue,and black radiate from the corners
    if (INIT_STYLE==1) {
      float Hmul,Wmul;
      for (int loop=0; loop<HEIGHT; loop++) {
        Hmul = (float)(loop)/(float)(HEIGHT)*5.0f; //[0..6)
        for (int loop2=0; loop2<WIDTH; loop2++) {
          Wmul = (float)(loop2)/(float)(WIDTH);   //[0..1]
          v_weights[loop2][loop].x = (1.0f-Wmul)*Hmul;
          v_weights[loop2][loop].y = Wmul*Hmul;
          v_weights[loop2][loop].z = Math.abs(Wmul)*(5.0f-Hmul);
        }//for HEIGHT
      }//for WIDTH
    }//if init_style

    //Have red,green,and blue form a circle around the center
    if (INIT_STYLE==2) {
      fpoint center = new fpoint((float)(WIDTH),(float)(HEIGHT),0.0f);
      fpoint outer = new fpoint(0.0f,0.0f,0.0f);

      float max_dist=get_dist(center,outer)/5.0f;
      float theta1=90.0f*((float)Math.PI/180.0f);
      float theta2=210.0f*((float)Math.PI/180.0f);
      float theta3=330.0f*((float)Math.PI/180.0f);
      float H2=(float)(HEIGHT/2),H4=(float)(HEIGHT/4),W2=(float)(WIDTH/2);

      fpoint rcenter = new fpoint((float)((Math.cos(theta1))*H4),(float)((Math.sin(theta1))*H4),0.0f);
      fpoint gcenter = new fpoint((float)((Math.cos(theta2))*H4),(float)((Math.sin(theta2))*H4),0.0f);
      fpoint bcenter = new fpoint((float)((Math.cos(theta3))*H4),(float)((Math.sin(theta3))*H4),0.0f);

      rcenter.set(rcenter.x+W2,rcenter.y+H2,0.0f);
      gcenter.set(gcenter.x+W2,gcenter.y+H2,0.0f);
      bcenter.set(bcenter.x+W2,bcenter.y+H2,0.0f);

      for (int loop=0; loop<HEIGHT; loop++) {
        for (int loop2=0; loop2<WIDTH; loop2++) {
          outer.set((float)(loop2),(float)(loop),0.0f);

          v_weights[loop2][loop].x = (get_dist(outer,rcenter)/max_dist);
          v_weights[loop2][loop].y = (get_dist(outer,gcenter)/max_dist);
          v_weights[loop2][loop].z = (get_dist(outer,bcenter)/max_dist);
        }//for WIDTH
      }//for HEIGHT
    }//if init_style
  }// init_map(int)


/*
  Initialize the colors to be organized
*/
  public void init_v_samples() {
    for (int loop=0; loop<MAX_PTS; loop++)
      {
      v_samples[loop].x = (float)(r.nextInt(6));
      v_samples[loop].y = (float)(r.nextInt(6));
      v_samples[loop].z = (float)(r.nextInt(6));
      }
  }


/*
  Returns a random sample
*/
  fpoint get_random() {
    int index = r.nextInt(MAX_PTS);

    return v_samples[index];
  }

  
/*
  Returns the best matching unit (bmu). In order to determine the bmu
    a distance is calculated from each weight to the randomly chosen
    sample. The weight with the smallest distance is closest to the 
    sample color. If there are more than one winners with the same 
    distance, then one is selected randomly
*/
  ipoint get_bmu(fpoint input) {  
    int match_amt=0;
    int WH=WIDTH*HEIGHT;
    float max_dist=10000.0f;
    ipoint match_list[] = new ipoint[WIDTH*HEIGHT];


    for (int loop=0; loop<WIDTH*HEIGHT; loop++)
      match_list[loop] = input.new ipoint();

    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
        float t_dist=fget_dist(input,v_weights[loop2][loop]);

        if (t_dist<max_dist) {
          max_dist = t_dist;
          match_list[0].set(loop2,loop,0);
          match_amt=1;
          }
        else if (t_dist==max_dist && match_amt<(WH))
          match_list[match_amt++].set(loop2,loop,0);
    }

    return match_list[r.nextInt(2147483647)%match_amt];
    }//get_bmu


/*
  Scales the neighboring weights. There are two parts to this operation:
    determining the neighbors and determining how much the neighbors will 
    learn. There are many ways to go about doing this, but I chose to use 
    gaussian function. The amount of neighbors and amount each weight can 
    learn all fall off with time
*/
  public void scale_neighbors(ipoint loc, fpoint factual, float t2) {

    int R2 = Math.round(((float)(RADIUS)*(1.0f-t2))/2.0f);
    fpoint outer  = new fpoint((float)(R2),(float)(R2),0.0f);
    fpoint center = new fpoint(0.0f,0.0f,0.0f);
    float d_normalize = get_dist(center,outer);

    for (int loop=-R2; loop<R2; loop++)
      for (int loop2=-R2; loop2<R2; loop2++)
        if ((loop+loc.y)>=0 && (loop+loc.y)<HEIGHT && (loop2+loc.x)>=0 && (loop2+loc.x)<WIDTH) {
        
          //Get distance from center point and normalize it
          outer.set((float)(loop2),(float)(loop),0.0f);
          float distance = get_dist(outer,center);
          distance/= d_normalize;

          //Get how much to scale it by
          float t=(float)(Math.exp(-1.0f*(Math.pow(distance,2.0f))/0.15f));

          //Amount a neuron can learn decreases with time
          //The 4 is chosen and the +1 is to avoid divide by 0's
          t/=(t2*4.0f+1.0f);

          //Scale it with the parametric equation
          fpoint temp = (factual.mult(t)).add(v_weights[loc.x+loop2][loc.y+loop].mult(1.0f-t));
          v_weights[loc.x+loop2][loc.y+loop] = temp;
          }
  }//scale neighbors


/*
  Updates the screen if in color mode by copying in the proper
    color into the viral screen
*/  
  public void update_rgb() {
    ipoint t = get_random().new ipoint();

    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
    	  
        t.x = v_weights[loop2][loop].fti().x;
        t.y = v_weights[loop2][loop].fti().y;
        t.z = v_weights[loop2][loop].fti().z;

        vrscr[loop2][loop] = cur_pal[rgb_table[t.x][t.y][t.z]];
      }
  }

  
/*
  Updates the screen if in black and white mode by computing the sum
    of the distances that the neighboring weights make from each weight
*/
  public float update_bw() {
    float t_map[][] = new float[100][100];
    float max_dist=0.0f;

    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
        fpoint center = new fpoint();
        center = v_weights[loop2][loop];

        int numinave=0;
        float total=0.0f;

        //Total up the distances
        for (int loopi=-SIMWEIGHT; loopi<=SIMWEIGHT; loopi++)
          for (int loopi2=-SIMWEIGHT; loopi2<=SIMWEIGHT; loopi2++)
            if ((loopi +loop)>=0   && (loopi +loop)<HEIGHT &&
               (loopi2+loop2)>=0  && (loopi2+loop2)<WIDTH){
             total+=get_dist(v_weights[loop2+loopi2][loop+loopi],center);
             numinave++;
             }

        //-1 is for the center, no cuenta
        total/=(float)(numinave-1);
        if (total>max_dist)
          max_dist=total;

        //Put all the totals into a buffer for later scaling
        t_map[loop2][loop]=total;
      }

    //The colors here are scaled in order to use all the colors
    float total=0.0f;
    for (int loop=0; loop<HEIGHT; loop++)
      for (int loop2=0; loop2<WIDTH; loop2++) {
        vrscr[loop2][loop]=cur_pal[255-Math.round((t_map[loop2][loop]/max_dist)*255.0f)];
        total+=t_map[loop2][loop];
        }

    //total is now the sum of the average distances which is a value
    //  describing how good the map is at showing similarities
    return total;
  }


/* 
  First copies all the correct colors into vrscr, then draws the
    colors onto the screen
*/
  public void paint(Graphics g) {
    int x = 15;

    if (INIT_STYLE!=-1) {
      if (RGB==1)
        update_rgb();
      else
        update_bw();
      
      for (int loop=0; loop<HEIGHT; loop++)
        for (int loop2=0; loop2<WIDTH; loop2++) {
          g.setColor(vrscr[loop2][loop]);
          g.fillRect((loop2<<2)+x,(loop<<2),4,4);
        }
    }//if INIT_STYLE    
    /*
    g.setColor(Color.black);
    g.fillRect(x,y,200,200);
    g.setColor(Color.white);
    g.fillRect(x,y,32,32);
    g.fillRect(x+200-32,y,32,32);
    g.fillRect(x,y+200-32,32,32);
    g.fillRect(x+200-32,y+200-32,32,32);
    */
  }//paint(Graphics)
}