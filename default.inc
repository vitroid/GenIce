// -*- c++ -*-
#include "glass.inc"


      /*
  global_settings {
    radiosity {
      pretrace_start 0.08
      pretrace_end   0.04
      count 35

      nearest_count 5
      error_bound 1.8
      recursion_limit 3

      low_error_factor 0.5
      gray_threshold 0.0
      minimum_reuse 0.015
      brightness 1

      adc_bailout 0.01/2
    }
  }
*/
global_settings { max_trace_level 100 } 

light_source { <0,0,0> 1 } //center


#declare enable_focus=0;
#declare draw_detail = 0;
#declare draw_media=1;
#declare enable_radiosity=1;
#declare enable_focus=1;


#if (enable_radiosity=1)
global_settings {
  radiosity {   }
}
#end

camera {
  location <0, 3, 0>
    look_at <0,0,0>
    right x
#if (enable_focus=1)
    focal_point < 0,0,0>
    aperture 0.05    
    blur_samples 30
#end
    angle 90  //2.7
    }

#declare INTOH = interior { ior 1.0 }
#declare INTHB = interior { 
ior 1.1  //water
}


#declare FA = finish {
    diffuse 1
    ambient 0
    specular 1
    roughness 0.005
}

#declare TEXO = texture {
  pigment {
    rgb <1,0,0.05>
  }
  finish {FA}
}
#declare TEXOH = texture {
  pigment {
    rgb <1,1,1>
  }
  finish {FA}
}
#declare TEXH = texture {
  pigment {
    rgb <0,0.7,1>
  }
  finish {FA}
}

//#declare TEXH = TEXO;
#declare TEXMe = TEXO;
//#declare TEXOH = TEXO;


#declare TEXHB = texture {
pigment { rgb <0.9,1,0.2> transmit 0 }
finish {FA}
/* for clear finish    ambient 0
    diffuse 1
    specular 1
    roughness 0.02
    }*/
}

#declare MATHB = material {
texture { TEXHB }
interior { INTHB }
/*
texture {
  pigment { rgbt 1 }

  finish {F_Glass1}
  }*/
//interior {I_Glass1}
}
  
#declare MATOH = material {
texture { TEXOH }
}

#declare MATO = material {
texture { TEXO }
}
#declare MATH = material {
texture { TEXH }
}

#declare MATMe = MATO;


//background { rgb <0.01,0.01,0.01> }


#declare RO = 0.05;
#declare RMe = 0.05;
#declare RH = 0.03;
#declare ROH = 0.02;
#declare RHB = 0.015;


/*
 sphere // transparent sphere containing media
 { 0,1.5 pigment { rgbt 1 } hollow
   interior
   { media
     { emission 1
       density
       { spherical density_map
         { [0 rgb 0]
	     //[0.3 rgb <0,0,1>]
           [0.6 rgb <0.5,1,1>]
           [1 rgb 1]
         }
       }
     }
  }
 }
*/
 // Scattering media box:
/*
 box
 { <-3, -3, -3>, <3, 3, 3>
   pigment { rgbt 1 } hollow
   interior
   { media
     { scattering { 1, 0.07 extinction 0.01 }
       samples 30
     }
   }
 }
*/


#if (draw_media=1)
sphere{ 0, 3.5 
pigment { rgbt 1 } hollow
  interior {
  media {
       scattering { 1, 0.01 extinction 0.01 } //normal
    #if (draw_detail=1)
      intervals 120
      samples 1, 10
      confidence 0.9999
      variance 1/1000
    #else
      intervals 10
      samples 1, 1
      confidence 0.9
      variance 1/100
    #end
    //ratio 0.9
  }
}}
#end
