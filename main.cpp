#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TRandom.h>
#include "TMinuit.h"
#include "TMultiGraph.h"

using namespace std;

void simulate(int numOfElect, double z0);

int main()
{
	int n = 200000;
	double z = 0.0;
	simulate(n, z);
}

void simulate(int numOfElect, double z0)
{
	// Exit coordinates
	TH1D *exit_x = new TH1D("exit_coordinate_x", "exit_coordinate_x", 500, -10, 10);
	TH1D *exit_y = new TH1D("exit_coordinate_y", "exit_coordinate_y", 500, -10, 10);
	TH1D *exit_r = new TH1D("exit_coordinate_r", "exit_coordinate_r", 500, 0, 10);

	// Constants
	double T = 163.15; // temperature in kelvin (-110 Celsius)
	double k = 8.617333e-5; // Boltzmann constant in eV K^-1
	// Boundary Conditions
	int upperBoundary = 4000;
	int lowerBoundary = -4000;

	map<double, vector<double>> mfp_table; // Table with stored mfp
	string line;
	ifstream file ("mfp.txt");
	if (!file.is_open())
	{
		cout <<"Error loading file." << endl;
		return;
	}
	// Read in mean free pass table line by line
	while (getline(file, line))
	{
		vector<double> l;
		double a, b, c, d, e;
		istringstream iss(line);
		iss >> a >> b >> c >> d >> e;
		l.push_back(b);
		l.push_back(c);
		l.push_back(d);
		l.push_back(e);
		mfp_table.insert(pair<double, vector<double>>(a, l));
	}
	// This is where simulation of electron diffusion starts
	for (int i = 0; i < numOfElect; i++)
	{
		vector<double> x;
		vector<double> y;
		vector<double> z;
		x.push_back(0);
		y.push_back(0);
		z.push_back(z0);
		int shift = 0;
		int j = 1;
		while (true) // Keep looping until electron reaches exit
		{
			int p;
			double mfp;
			double prob = ((double) rand()/ (RAND_MAX)); // Uniform distribution from 0 to 1
		  double boltzmann = gRandom->Exp(k*T); // Boltzmann distribution
			// decide to round either by the nearest tenth or hundreth
			if (boltzmann < 0.1)
				p = 2;
			else
				p = 1;
			// Doing the rounding
			double value = (int)(boltzmann * pow(10, p) + 0.5);
			double energy = (double)value / pow(10, p);
			if (energy == 0) // No electrons with zero energy
				energy = 0.01;
			if (prob <= mfp_table[energy][1])
				mfp = mfp_table[energy][0];
		  else
				mfp = mfp_table[energy][2];
			double r = gRandom->Exp(mfp); // Exponential distribution with the average equal to mean free path
			// Choose a random direction for electron to travel in cylindrical coordinates
			double theta = gRandom->Uniform(TMath::Pi());
			double phi = gRandom->Uniform(2*TMath::Pi());
			double z_check = r * TMath::Cos(theta) + z[j+shift-1];
			// if true then calculate where electron hits the upper boundary
			if (z_check > upperBoundary)
			{
				double r2 = (upperBoundary - z[j+shift-1])/TMath::Cos(theta);
				x.push_back(r2 * TMath::Sin(theta) * TMath::Cos(phi) + x[j+shift-1]);
				y.push_back(r2 * TMath::Sin(theta) * TMath::Sin(phi) + y[j+shift-1]);
				z.push_back(r2 * TMath::Cos(theta) + z[j+shift-1]);
				shift++;
				r = r - r2;
				theta = TMath::Pi() - theta;
				z_check = r * TMath::Cos(theta) + z[j+shift-1];
			}
			if (z_check >= lowerBoundary)
			{
				x.push_back(r * TMath::Sin(theta) * TMath::Cos(phi) + x[j+shift-1]);
				y.push_back(r * TMath::Sin(theta) * TMath::Sin(phi) + y[j+shift-1]);
				z.push_back(z_check);
			}
			else // Electron has reached lower boundary
			{
				double r2 = (lowerBoundary - z[j+shift-1])/TMath::Cos(theta);
				x.push_back(r2 * TMath::Sin(theta) * TMath::Cos(phi) + x[j+shift-1]);
				y.push_back(r2 * TMath::Sin(theta) * TMath::Sin(phi) + y[j+shift-1]);
				z.push_back(r2 * TMath::Cos(theta) + z[j+shift-1]);
				double e_x = x[j+shift]/8000;
				double e_y = y[j+shift]/8000;
				exit_x->Fill(e_x);
				exit_y->Fill(e_y);
				exit_r->Fill(TMath::Sqrt(TMath::Power(e_x, 2) + TMath::Power(e_y, 2)));
				break;
			}
			j++;
		}
	}

	TCanvas *r_plot = new TCanvas("Exit r coordinates", "Exit r coordinates", 30, 10, 900, 800);
	r_plot->SetBorderMode(0);
	r_plot->SetFrameFillColor(0);
	r_plot->Divide(1,1);
	r_plot->SetFillStyle(4100);
	r_plot->SetFillColor(0);

	r_plot->cd(1);
	gPad->SetLogy();
	exit_r->GetXaxis()->SetTitle("r/dff");
	exit_r->GetYaxis()->SetTitle("Entries");
	exit_r->Draw();

	r_plot->Update();
	r_plot->Print("Exit Coordinates r");

	TCanvas *x_plot = new TCanvas("Exit x coordinates", "Exit x coordinates", 30, 10, 900, 800);
	x_plot->SetBorderMode(0);
	x_plot->SetFrameFillColor(0);
	x_plot->Divide(1,1);
	x_plot->SetFillStyle(4100);
	x_plot->SetFillColor(0);

	x_plot->cd(1);
	gPad->SetLogy();
	exit_x->GetXaxis()->SetTitle("x/dff");
	exit_x->GetYaxis()->SetTitle("Entries");
	exit_x->Draw();

	x_plot->Update();
	x_plot->Print("Exit Coordinates x");

	TCanvas *y_plot = new TCanvas("Exit y coordinates", "Exit y coordinates", 30, 10, 900, 800);
	y_plot->SetBorderMode(0);
	y_plot->SetFrameFillColor(0);
	y_plot->Divide(1,1);
	y_plot->SetFillStyle(4100);
	y_plot->SetFillColor(0);

	y_plot->cd(1);
	gPad->SetLogy();
	exit_y->GetXaxis()->SetTitle("y/dff");
	exit_y->GetYaxis()->SetTitle("Entries");
	exit_y->Draw();

	y_plot->Update();
	y_plot->Print("Exit Coordinates y");

}
