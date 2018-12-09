// random_dist - Program to compute and plot a gaussian distribution of 
//               random numbers on the interval [-6*sigma,6*sigma]
#include <iostream>
#include "math.h"
#include "TStopwatch.h"
#include "TRandom.h"

// Random number generator; Uniform dist. in [0,1]
double rand_ejb( long& seed ) {
//                                                                          
// Input:                                                                           
//   	seed    Integer seed (DO NOT USE A SEED OF ZERO)
// Output:
//   	rand    Random number uniformly distributed in [0,1)
//
// Numerical Recipes
  const double a = 1664525.0;
  const double m = 2147483648.0*2.0;
  const double c = 1013904223.0;

  double temp = a * seed + c;
  seed = (long)(fmod(temp,m));
  double rand = seed/m;
  return(rand);
}

void graph_random_dist_gaussian_compare() {
//
//Author : Edward J. Brash
//	       

int npoints, i;
double r1,r2,x,x0,pgauss,sigma;
long seed = 1;
    
cout << "Enter the number of data points to generate:" << endl;
cin >> npoints;

cout << "Enter the centroid of the gaussian:" << endl;
cin >> x0;

cout << "Enter the sigma of the gaussian:" << endl;
cin >> sigma;

const int nbins=1000;
double width = 12*sigma;
double dx = width/nbins;
int n_randoms=0;

TH1F *h1 = new TH1F("h1","Random Number Distribution (EJB)",nbins,x0-width/2.0,x0+width/2.0);
TH1F *h2 = new TH1F("h2","Random Number Distribution (ROOT)",nbins,x0-width/2.0,x0+width/2.0);

TStopwatch timer;
TStopwatch timer3;
timer.Start();

for (i=0;i<npoints; i++){
	bool check=true;
	while(check){
		r1 = rand_ejb(seed);
		n_randoms++;
		x=(x0-width/2.0+width*r1);
		pgauss=exp(-1.0*(x-x0)*(x-x0)/(2.0*sigma*sigma));
		r2 = rand_ejb(seed);
		n_randoms++;
		//cout << "Random number: " << i << " = " << seed << " = " << pgauss << " " << r2 << endl;
		if (r2<pgauss) check=false;
	}	
	h1->Fill( x , 1.0/width/sqrt(2.0*3.14159265*sigma*sigma) );
	//h1->Fill(x,1.0);
}
   
timer.Stop();
TRandom3 *rnum = new TRandom3();
timer3.Start();

for (i=0;i<npoints; i++){
	x=rnum->Gaus(x0,sigma);
	h2->Fill( x , 1.0/width/sqrt(2.0*3.14159265*sigma*sigma) );
	//h1->Fill(x,1.0);
}
timer3.Stop();

TH1F *h3 = (TH1F*)h1->Clone("h3");
h3->Divide(h2);

// Plotting stuff below this point ...

   const Int_t n = npoints;
   cout << "Total number of randoms generated = " << n_randoms << endl;
   cout << endl;
   cout << "EJB Random Generation" << endl;
   cout << "Total real time = " << timer.RealTime() << endl;
   cout << "Total CPU time = " << timer.CpuTime() << endl;
   cout << "Real time / random = " << timer.RealTime()/npoints << endl;
   cout << "CPU time /random = " << timer.CpuTime()/npoints << endl;
   cout << endl;
   cout << "ROOT Random Generation" << endl;
   cout << "Total real time = " << timer3.RealTime() << endl;
   cout << "Total CPU time = " << timer3.CpuTime() << endl;
   cout << "Real time / random = " << timer3.RealTime()/npoints << endl;
   cout << "CPU time /random = " << timer3.CpuTime()/npoints << endl;
   cout << endl;

   TCanvas *c2 = new TCanvas("c2","Cubic Data",200,10,700,800);
   TCanvas *c1 = new TCanvas("c1","Cubic Data",200,10,700,800);
   c1->Divide(1,2);
   c1->SetFillColor(42);
   c1->SetGrid();
   c2->SetFillColor(42);
   c2->SetGrid();
   //gStyle->SetOptFit(kFALSE);	
   gStyle->SetOptFit(1);	
   
   c2->cd();
   h3->Draw();
   c1->cd(1);
   h1->Fit("gaus");
   h1->Draw();
   c1->cd(2);
   h2->Fit("gaus");
   h2->Draw();
   
   // TCanvas::Update() draws the frame, after which one can change it
   c2->Update();
   c2->GetFrame()->SetFillColor(21);
   c2->GetFrame()->SetBorderSize(12);
   c2->Modified();
   c1->Update();
   c1->GetFrame()->SetFillColor(21);
   c1->GetFrame()->SetBorderSize(12);
   c1->Modified();
}
