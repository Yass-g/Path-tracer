/**
    main.cpp
    outputs a .ppm file, representing a scene by solving the rendering equation (Monte Carlo)
    the number of samples is set as an argument (default : 4*smp = 4)
    @author Yass-g
    @version 1.2 10/04/21 
*/
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <random>
#include <vector>
#include "geometry.h"

using namespace std;
default_random_engine generator;
uniform_real_distribution<double> distr(0.0,1.0);
#define M_PI 3.1415926

Vect tc(0.0588, 0.361, 0.0941);
double R=60;
Vect Cen(50,40.8,-860);
double T=30*M_PI/180.;
double D=R/cos(T);
double Z=60;

vector<Object*> Scene = {
    new  Sphere(1600, Vect(1,0,2)*3000, Vect(.8,.8,.8)*1e100,Vect(), DIFF), // Light
    //new  Sphere(1560, Vect(1,0,2)*3500,Vect(.4,.4,.4)*2e100, Vect(),  DIFF),
    new  Sphere(110000, Vect(50, -110048.5, 0),  Vect(.5,.5,.5)*4,Vect(),DIFF),
    new  Sphere(10000,Vect(50,40,-1000), Vect(0.00063842, 0.02001478, 0.28923243)*5e-1, Vect(.3,.3,1)*.50,  DIFF), // Backgrnd
    new  Sphere(100000, Vect(50, -100000, 0),  Vect(),Vect(.3,.3,.3),SPEC), // floor
    new  Sphere(26.5,Vect(60,26.5,22),   Vect(),Vect(1,1,1)*0.6, DIFF),
    //new  Sphere(13,Vect(22,13,62),   Vect(),Vect(.96,.96,.96), IMP_SPEC),
    //new  Sphere(8,Vect(80,8,80),   Vect(),Vect(.66,1,.96), REFR),
    //new  Sphere(10,Vect(60,10,80),   Vect(),Vect(.66,.66,.96), DIFF),
};

bool intersect(const Ray &r, double &t, int &id){ // closest intersection with object in scene
  double n= Scene.size(), d, inf=t=1e20;
  for(int i=int(n);i--;)
  {
      d=Scene[i]->intersect(r);
      if(d&&d<t)
        t=d;id=i;
  }
  return t<inf;
}

Vect rendering(const Ray &r_, int depth_){ //solving the rendering equation (Monte Carlo+Russian roulette)
  double t;
  int id=0;
  Ray r=r_;
  int depth=depth_;
  Vect cl(0,0,0);   // color vector
  Vect cf(1,1,1);  //reflectance vector
  while (1){
    if (!intersect(r, t, id)) return cl; // if miss, return black
    const auto obj = Scene[id];        // the hit object
    Vect x=r.o+r.d*t;//x : intersection point,
    Vect n=obj->normal(x);// n : sphere normal,
    Vect nl=n.dot(r.d)<0?n:n*-1; //n1 : corrected normal,
    Vect f=obj->c; //  f : object color

    Ray ref_Ray(x, r.d-n*2*n.dot(r.d));//reflected ray
    double p = (f.x>f.y && f.x>f.z) ? f.x : (f.y>f.z) ? f.y : f.z; // max reflection
    cl = cl + cf.mult(obj->e);
    if (++depth>5)//Russian roulette after (at least) a depth of 5
        if (erand48(1)<p)
            f=f*(1/p);
        else return cl;
    cf = cf.mult(f);
    if (obj->refl == DIFF){ // DIFFUSE reflection
        double r1=2*M_PI*erand48(1);
        double r2=erand48(1);
        double r2s=sqrt(r2);
        Vect w=nl, u=((fabs(w.x)>.1?Vect(0,1):Vect(1))%w).norm(), v=w%u;
        Vect d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
        r = Ray(x,d);
        continue;
    }
    else if (obj->refl == SPEC){           // SPECULAR reflection
        r = ref_Ray;
        continue;
    }
     else if (obj->refl == IMP_SPEC){           // SPECULAR reflection
        double fi=2*M_PI*erand48(1);
        double r2=erand48(1);
        double expo = 4;
        double costhet = pow(r2, 1/expo);
        double sinthet = sqrt(1-costhet*costhet);
        Vect w=r.d-n*2*n.dot(r.d), u=((fabs(w.x)>.1?Vect(0,1):Vect(1))%w).norm(), v=w%u;
        Vect d = (u*cos(fi)*sinthet + v*sin(fi)*sinthet + w*costhet).norm();
        r = Ray(x,d);
        continue;
     }
                                    // Ideal dielectric REFRACTION
        bool into = n.dot(nl)>0;                // normal orientation
        double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
        if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0){    // Total internal reflection
              r = ref_Ray;
              continue;
            }
        Vect tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        if (erand48(1)<P){
                cf = cf*RP;
                r = ref_Ray;
        } else {
                cf = cf*TP;
                r = Ray(x,tdir);
        }

  }
}


int main(int argc, char *argv[]){
  int w=1024, h=768, smp = 2; // resolution and samples
  Ray cam(Vect(50,25,295), Vect(0,0,-1).norm()); // camera settings
  Vect cx=Vect(w*0.5135/h), cy=(cx%cam.d).norm()*0.5135, r, *out=new Vect[w*h];
  int N = 150;
    // generate scene
  for(int i=0; i<N; i++)
  {
      double r = 2*erand48(1) + 2.0;
      double x = erand48(1)*170.0 -50;
      double y = -erand48(1)*500.0+100.0;
      double red = erand48(1);
      double  green = erand48(1);
      double blue = erand48(1);
      int m = rand()%4;
      cout<<r<<"  "<<x<<"  "<<y<<"  "<<m<<"  "<<endl;
      Scene.push_back(new  Sphere(r,Vect(x,r,y),   Vect(),Vect(i%2,(i+1)%2,rand()%2), (m==0)?DIFF:(m==1)?SPEC:(m==1)?IMP_SPEC:REFR));
  }
   // generating outage
  for (int y=0; y<h; y++){
    fprintf(stderr,"\rDone at (%d spp) %5.2f%%",smp*4,100.*y/(h-1));
    for (unsigned short x=0; x<w; x++)
        //loop over 2x2 subpixels (antialiasing)
      for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)
        for (int sx=0; sx<2; sx++, r=Vect()){
          for (int s=0; s<smp*2; s++){
            double r1=2*erand48(1), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
            double r2=2*erand48(1), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
            Vect d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) + cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
            r = r + rendering(Ray(cam.o+d*140,d.norm()),0)*(1./smp);
          }
          out[i] = out[i] + Vect(bound(r.x),bound(r.y),bound(r.z))*.25;
        }
  }
  FILE *f = fopen("metaltest8smp+.ppm", "w");
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(f,"%d %d %d ", toInt(out[i].x), toInt(out[i].y), toInt(out[i].z));
}
