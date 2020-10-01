/* 
Copyright 2020 ContacTUM
Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
and associated documentation files (the "Software"), to deal in the Software without 
restriction, including without limitation the rights to use, copy, modify, merge, publish, 
distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom 
the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or 
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.
*/

/*
 * This code reads in the files generated from CVMC and calculates summary outcomes.
 * Results are written into a new tree.
 * @author Tina Pollmann
 * @date 2020
 */


//c++
#include <list>
#include <vector>
#include <numeric>
#include <iostream>
#include <fstream>
#include <ctime>

//root
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TPaveText.h>
#include <TStopwatch.h>

using namespace std;
bool fDebug = 0; //! Enable or disable debug output
bool fWriteTextfile = 0; //!Write results to a text file instead of to the default .root output file.

//! Histograms to hold information of key outcomes over each set of MC runs
int nbins = 400;
TH1F *hMaxDailySick  = new TH1F("hMaxDailySick",";max sick per day",nbins,0,1);
TH1F *hMaxQuarantine = new TH1F("hMaxQuarantine",";max quarantined per day", nbins,0,1);
TH1F *hMaxTotalSick = new TH1F("hMaxTotalSick",";max total sick", nbins,0,1);
TH1F *hSick1yr = new TH1F("hSick1yr","; sick after 1 year", nbins,0,1);
TH1F *hQuarantineSum = new TH1F("hQuarantineSum","; sick after 1 year", nbins,0,1);
TH1F *hDoublingTime = new TH1F("hDoublingTime","; doubling time [days]", 240,-60,60);
TH1F *hRe = new TH1F("hRe","; hRe", nbins,0,10);
TH1F *hRePre = new TH1F("hRePre","; hRe", nbins,0,10);


TGraph *gNSusceptible[100]; //! The array length should be equal or bigger than the number of simulation run.
TGraph *gNExposed[100];
TGraph *gNReported[100];
TGraph *gNQuarantineToday[100];
TGraph *gNQuarantineNotInfectious[100];


void makeLegend(TLegend *ll, TVirtualPad *pad, int location);


TF1 *fexp = new TF1("fexp","[0]*TMath::Exp(x*TMath::Log(2)/[1])",0,500); //! Fit function to get the doubling time
int cc2 = 0;

//! This fits an exponential to the curve of 'newly infected each day' to extract the doubling time
//! It tries different starting values for the doubling time from -32 to 80 days
//! and returns as soon as the reduced chi2 is better than 3.
float CalculateDoublingTime(TH1F *hinfected, int tstart, int tend) {
  if (!hinfected || hinfected->GetEntries() < 10) return -99999;
    int total = 0;
  for (int ibin = hinfected->FindBin(tstart); ibin < hinfected->FindBin(tend); ibin++) {
    total = total + hinfected->GetBinContent(ibin);
  }
  if (total < 50) return -99999;
	int maxbin = hinfected->GetMaximumBin();
	int maxday = hinfected->GetXaxis()->GetBinLowEdge(maxbin);
	int maxnumber = hinfected->GetBinContent(maxbin);
	int fitToDay = 0.0;
	int fitFromDay = -10;
	float fitToFraction = 0.6;
	int killswitch = -16;
	float redchi = 10.;
	while (redchi > 3 && killswitch < 40) {
    float height = (hinfected->GetBinContent(hinfected->FindBin(tstart)));	
    fexp->SetParameters(1,killswitch*2.0);
    float A0 = height/fexp->Eval(tstart);
    fexp->SetParameters(A0,killswitch*2.0);
    int status = hinfected->Fit("fexp","q","",tstart, tend); // q: quiet (do not print fit status)
    redchi = fexp->GetChisquare()/fexp->GetNDF();    
    killswitch++;
	}
	cc2++;
	if (killswitch < 39) { return fexp->GetParameter(1);}
	return  -99999;
}

TCanvas *tt = new TCanvas("fit","fit");

//! All the work happens here
//! fin: The output file from the CVMC
//! runner: index used to append to debug output so that each output file a unique name
//! fname: The name for the default output .root file
void postprocess(TFile *fin, int runner, string fname) {
	if (!fin->IsOpen()) return; //! Can't do anything if the file isn't open.
	gStyle->SetOptTitle(0);
	
	//! Reset all the histograms that contain the outcomes
	hMaxDailySick->Reset();
	hMaxQuarantine->Reset();
	hMaxTotalSick->Reset();
	hSick1yr->Reset();
	hQuarantineSum->Reset();
	hDoublingTime->Reset();
	hRe->Reset();
	hRePre->Reset();
	
	//! Get from the input file the tree with the daily counts ...
	TTree* PopulationLevelInformation = (TTree*)fin->Get("PopulationLevelInformation");
	//! .. and the full list of settings this simulation was run with ... 
	TTree* settings = (TTree*)fin->Get("settings")	;
	if (PopulationLevelInformation == NULL || settings == NULL) return;

  //! ... and make all the variables available for use
	int fFileIndex, fDayForROOTTree, fNSusceptible, fNExposed, fNInfectious, fNRecovered, fNTraced, fNReported,fNQuarantineToday;
	float feffectiveR, feffectiveRUncertainty;
	PopulationLevelInformation->SetBranchAddress("fileIndex",&fFileIndex);
	PopulationLevelInformation->SetBranchAddress("day",&fDayForROOTTree);
	PopulationLevelInformation->SetBranchAddress("fNSusceptible",&fNSusceptible);
	PopulationLevelInformation->SetBranchAddress("fNExposed",&fNExposed);
	PopulationLevelInformation->SetBranchAddress("fNInfectious",&fNInfectious);
	PopulationLevelInformation->SetBranchAddress("fNRecovered",&fNRecovered);
	PopulationLevelInformation->SetBranchAddress("fNTraced",&fNTraced);
	PopulationLevelInformation->SetBranchAddress("fNReported",&fNReported);
	PopulationLevelInformation->SetBranchAddress("feffectiveR",&feffectiveR);
	PopulationLevelInformation->SetBranchAddress("feffectiveRUncertainty",&feffectiveRUncertainty);
	PopulationLevelInformation->SetBranchAddress("fNQuarantineToday",&fNQuarantineToday);	

	int fNPersons, fRunId, fNDays, fDaysInQuarantine, fTracingOrder;
	int fStartTestingOnDay, fDaysBackwardTrace, fdTTest, fSocialDistancingFrom,fSocialDistancingTo ;
	int fDaysToTestResult, fTracingDelay, fLastDayWithPatients;
	float fStartTracingTestingInfectedFraction, TestThreshold, fAppProbability, symptomProbability, testPositiveProbability, fReportingProbability;
	float falsePositiveRate, transmissionProbability, fTracingEfficiency, incubationGamma, incubationMu, incubationBeta, infectionGamma, infectionMu, infectionBeta;
	float fPeopleMetPerDay, R00, R0e, AsymptomaticTransmissionScaling, fSocialDistancingFactor, fSocialDistancingMaxPeople, RandomTestingRate;
	bool fBackwardTracing, fTraceUninfected, RandomTesting;
  settings->SetBranchAddress("index",&fRunId);
  settings->SetBranchAddress("nPersons",&fNPersons);
  settings->SetBranchAddress("nDays",   &fNDays);
  settings->SetBranchAddress("LastDayWithPatients", &fLastDayWithPatients);  
  settings->SetBranchAddress("peopleMetPerDay", &fPeopleMetPerDay);
  settings->SetBranchAddress("daysInQuarantine",&fDaysInQuarantine);
  settings->SetBranchAddress("tracingOrder",    &fTracingOrder);
  settings->SetBranchAddress("startTracingTestingInfectedFraction", &fStartTracingTestingInfectedFraction);
  settings->SetBranchAddress("startTracingOnDay", &fStartTestingOnDay);
  settings->SetBranchAddress("daysBackwardTrace", &fDaysBackwardTrace);
  settings->SetBranchAddress("backwardTracing", &fBackwardTracing);
  settings->SetBranchAddress("traceUninfected", &fTraceUninfected);
  settings->SetBranchAddress("dTTest", &fdTTest);
  settings->SetBranchAddress("testThreshold", &TestThreshold);
  settings->SetBranchAddress("AppProbability", &fAppProbability);
	settings->SetBranchAddress("tracingEfficiency", &fTracingEfficiency);  
  settings->SetBranchAddress("daysToTestResult",  &fDaysToTestResult);
  settings->SetBranchAddress("symptomProbability",     &symptomProbability);   
  settings->SetBranchAddress("ReportingProbability", &fReportingProbability);   
  settings->SetBranchAddress("AsymptomaticTransmissionScaling",&AsymptomaticTransmissionScaling); 
  settings->SetBranchAddress("TracingDelay", &fTracingDelay);
  settings->SetBranchAddress("randomTestingRate", &RandomTestingRate);
  settings->SetBranchAddress("randomTesting", &RandomTesting);  
  settings->SetBranchAddress("testPositiveProbability",&testPositiveProbability);
  settings->SetBranchAddress("falsePositiveRate",      &falsePositiveRate);
  settings->SetBranchAddress("transmissionProbability",&transmissionProbability);  
  settings->SetBranchAddress("incubationGamma",&incubationGamma); 
  settings->SetBranchAddress("incubationMu",   &incubationMu); 
  settings->SetBranchAddress("incubationBeta", &incubationBeta); 
  settings->SetBranchAddress("infectionGamma", &infectionGamma); 
  settings->SetBranchAddress("infectionMu",    &infectionMu); 
  settings->SetBranchAddress("infectionBeta",  &infectionBeta); 
  settings->SetBranchAddress("R00",  &R00);
  settings->SetBranchAddress("R0e",  &R0e);  
  settings->SetBranchAddress("SocialDistancingMaxPeople", &fSocialDistancingMaxPeople);
  settings->SetBranchAddress("SocialDistancingFrom", &fSocialDistancingFrom);
  settings->SetBranchAddress("SocialDistancingTo ", &fSocialDistancingTo);
  settings->SetBranchAddress("SocialDistancingFactor", &fSocialDistancingFactor);  
  settings->GetEntry(0);
  int NPersons = fNPersons;
	
	
  //! Define all the variables used to gather the statistics on the outcomes 	
	int lastindex = -1;
	int maxdailysick = 0, maxq=0, maxtotalsick = 0, sicktoday, qtoday, LastNSusceptible=NPersons;
	float fractionNoOutbreak = 0., sickByDay365, QuarantineManDays = 0.0, doublingTime= 0.0, avgRe =0.0 , avgRePre = 0.0, fractionNoOutbreak50 = 0.0;
	int filecounter = 0, recounter=0, reprecounter=0;
	//! Settings
	int calculateReFromInterventionStartToDay = 28;
	int numberSickToCountAsOutbreak = int(fStartTracingTestingInfectedFraction*float(NPersons));
  ofstream out, out2;
  if (fWriteTextfile) { out.open("Sum_CovidSimM5_3330.csv"); out2.open("Summary_CovidSimM5_3330.csv"); }
	if (fWriteTextfile) {
	  out << "day, fNSusceptible, fNExposed, fNInfectious, fNRecovered, fNTraced, fNReported, feffectiveR, feffectiveRUncertainty, fNQuarantineToday" << endl;
  	out2 << " fraction,  hMaxDailySick, hMaxQuarantine, hMaxTotalSick, hSick1yr, hQuarantineSum, hRe" << endl;
  }
	//! Loop on the PopulationLevelInformation tree to gather information on statistics in each 'box' of the disease.
	TH1F *hNewlyInfected = new TH1F("hNewlyInfected","hNewlyInfected", 250,0,500);
	for (int ientry = 0; ientry < PopulationLevelInformation->GetEntries(); ientry++) {
		PopulationLevelInformation->GetEntry(ientry);
		//! The daily numbers are in the same tree for all simulations run. This splits them up so that we gather statistics for each run.
		if ((lastindex != fFileIndex && ientry > 0) || ientry == PopulationLevelInformation->GetEntries()-1) {
			filecounter++;
			if (fWriteTextfile)out << endl;
		 	if (maxtotalsick < numberSickToCountAsOutbreak) { fractionNoOutbreak = fractionNoOutbreak + 1.0; }
		 	if (maxtotalsick < 50) { fractionNoOutbreak50 = fractionNoOutbreak50 + 1.0; }
		 	if ((maxtotalsick > 50 && fStartTracingTestingInfectedFraction < 0.0001) ||  maxtotalsick > numberSickToCountAsOutbreak) {
				hMaxDailySick->Fill(float(maxdailysick)/float(NPersons));
				hMaxQuarantine->Fill(float(maxq)/float(NPersons));
				hMaxTotalSick->Fill(float(maxtotalsick)/float(NPersons));
				hSick1yr->Fill(float(sickByDay365)/float(NPersons));
				hQuarantineSum->Fill(QuarantineManDays/float(NPersons)/365.);
				//! Calculation of doublingTime currently disabled for speed reasons.
				//doublingTime = CalculateDoublingTime(hNewlyInfected, fStartTestingOnDay, TMath::Min(fStartTestingOnDay+42, fLastDayWithPatients));
				doublingTime = -9999.999				;
				if (fDebug) hNewlyInfected->DrawCopy();
				if (fDebug) tt->SaveAs(Form("fit%d_%d.png",runner,filecounter));
				hDoublingTime->Fill(doublingTime); 
				hRe->Fill(avgRe/float(recounter));
				hRePre->Fill(avgRePre/float(reprecounter));
				if (fDebug) cout << "Re " << avgRePre/float(reprecounter) << ", " << avgRe/float(recounter) << endl;
			}
			settings->GetEntry(filecounter);
			maxdailysick = 0; maxq=0; maxtotalsick = 0; QuarantineManDays=0., doublingTime=0.0, avgRe =0.0, avgRePre = 0.0;
			recounter=0; reprecounter=0;
			hNewlyInfected->Reset(); LastNSusceptible=NPersons;			
		}		
	  if (fWriteTextfile) out << fDayForROOTTree << ", " << fNSusceptible<< ", " <<fNExposed<< ", " <<fNInfectious<< ", " <<fNRecovered<< ", " <<fNTraced<< ", " <<fNReported<< ", " <<feffectiveR<< ", " <<feffectiveRUncertainty<< ", " <<fNQuarantineToday << endl;		
		maxtotalsick = NPersons-fNSusceptible;
		sicktoday = fNExposed + fNInfectious;		
		if (sicktoday > maxdailysick) maxdailysick = sicktoday;
		if (fNQuarantineToday > maxq) maxq = fNQuarantineToday;		
		if (fDayForROOTTree <= 365+fStartTestingOnDay) sickByDay365 = maxtotalsick;
		if (fDayForROOTTree > fStartTestingOnDay+10 //! the first 10 days are still influenced by the people who got exposed before interventions started
		              && fDayForROOTTree < fStartTestingOnDay + calculateReFromInterventionStartToDay  //! we can't go too long past start of interventions lest we get into non-linear regime
		              && feffectiveR > -1 //! -1 means no one recovered today so we can't calculate Re
		              && maxtotalsick < NPersons/2.0)  //! do not go into non-linear regime. The good thing here is that if this many people are exposed, we have enough statistics to get a good Re even if the number of days we average over is smaller
		              { avgRe = avgRe+feffectiveR; recounter++;}
		if (fDayForROOTTree > fStartTestingOnDay-calculateReFromInterventionStartToDay/2 && fDayForROOTTree < fStartTestingOnDay - 1 && feffectiveR > -1)  { avgRePre = avgRePre+feffectiveR; reprecounter++;}
		if (fDayForROOTTree > fStartTestingOnDay && fDayForROOTTree < 365+fStartTestingOnDay) QuarantineManDays = QuarantineManDays + fNQuarantineToday;
		hNewlyInfected->SetBinContent(hNewlyInfected->FindBin(fDayForROOTTree), LastNSusceptible - fNSusceptible ); // this will quietly ignore out-of-range bins
		LastNSusceptible = fNSusceptible;
		lastindex = fFileIndex;
	}	
	fractionNoOutbreak = fractionNoOutbreak/float(filecounter);
	fractionNoOutbreak50 = fractionNoOutbreak50/float(filecounter);
	if (fWriteTextfile) {
    for (int ibin = 0; ibin < nbins; ibin++) {
      out2 << hMaxDailySick->GetBinCenter(ibin) << ", " << hMaxDailySick->GetBinContent(ibin) << ", " << hMaxQuarantine->GetBinContent(ibin) << ", " << hMaxTotalSick->GetBinContent(ibin) << ", " << hSick1yr->GetBinContent(ibin) << ", " << hQuarantineSum->GetBinContent(ibin) << endl;
     }
     out2 << endl;
    for (int ibin = 0; ibin < nbins; ibin++) { 
     out2 << hRe->GetBinCenter(ibin) << ", " << hRe->GetBinContent(ibin) << endl;
    }	
  }
	float averagemaxdailysick = hMaxDailySick->GetMean();
	float averagemaxdailysickWidth = hMaxDailySick->GetRMS();
	float averagemaxq = hMaxQuarantine->GetMean();
	float averagemaxqWidth = hMaxQuarantine->GetRMS();
	float averagemaxtotalsick = hMaxTotalSick->GetMean();
	float averagesickByDay365 = hSick1yr->GetMean();
	float averageQuarantineManDays = hQuarantineSum->GetMean();
	float averageQuarantineManDaysWidth = hQuarantineSum->GetRMS();
	float averagesickByDay365Width = hSick1yr->GetRMS();
	float averageDoublingTime = hDoublingTime->GetMean();
	float averageReAfterIntervention = hRe->GetMean();
	float averageReAfterInterventionWidth = hRe->GetRMS();
	float averageReBeforeIntervention = hRePre->GetMean();
	if (fDebug) cout << "T2 = " << averageDoublingTime << ", Re = " << averageReAfterIntervention << ", No luck " << fractionNoOutbreak << endl;
	if (fDebug) cout << " Re before = " << averageReBeforeIntervention << endl;
	if (fDebug) cout << "Drawing from PopulationLevelInformation" << endl;

if (fDebug) {
    for (int irun = 0; irun < filecounter; irun++) {			
      PopulationLevelInformation->Draw(Form("day:fNSusceptible/(%d):(fNExposed+fNInfectious)/(%d):fNInfectious/(%d)",fNPersons,fNPersons,fNPersons) , Form("fileIndex==%d",irun) ,"goff");
      gNSusceptible[irun] = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV2());
      gNExposed[irun] = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV3());
      //TGraph *gNInfectious = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV4());	
  
      PopulationLevelInformation->Draw(Form("day:fNReported/(%d):fNQuarantineToday/(%d):((fNQuarantineToday)/%d - (%d - fNExposed-fNInfectious)/%d)",fNPersons,fNPersons,fNPersons,fNPersons) , Form("fileIndex==%d",irun),"goff");
      gNQuarantineToday[irun] = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV3());	
      gNSusceptible[irun]->SetLineColor(1); gNSusceptible[irun]->SetMarkerColor(1);
      gNExposed[irun]->SetLineColor(kOrange-3); gNExposed[irun]->SetMarkerColor(kOrange-3);
      gNQuarantineToday[irun]->SetLineColor(kGreen+1); gNQuarantineToday[irun]->SetMarkerColor(kGreen+1);gNQuarantineToday[irun]->SetLineWidth(2);
       
      gNReported[irun] = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV2());	
      gNQuarantineNotInfectious[irun] = new TGraph(PopulationLevelInformation->GetSelectedRows(), PopulationLevelInformation->GetV1(), PopulationLevelInformation->GetV4());	
      gNQuarantineNotInfectious[irun]->SetLineColor(kCyan+4); gNQuarantineNotInfectious[irun]->SetMarkerColor(kCyan+4);
      gNReported[irun]->SetLineColor(kBlue+4); gNReported[irun]->SetMarkerColor(kBlue+4);
    }

    float yoffset = 2.0, xoffset = 2.3;
    TCanvas *tc1 = new TCanvas(Form("tc%d",runner),Form("tc%d",runner),1600,1000);
    tc1->Divide(2,3);
    tc1->cd(1);
    if (fDebug) cout << "Drawing Graphs" << endl;
    TH1F *hdummy = new TH1F("hdummy","; day; Fraction susceptible",10,0,365);
    hdummy->Draw();
    hdummy->GetYaxis()->SetRangeUser(0,1.01);
    hdummy->GetYaxis()->SetTitleOffset(yoffset);
    hdummy->GetYaxis()->SetTitle("Fraction susceptible");
    hdummy->GetXaxis()->SetTitleOffset(xoffset);
    //gNSusceptible[0]->Draw("samePL");
    for (int irun = 1; irun < filecounter; irun++) {	 gNSusceptible[irun]->Draw("samePL");	 }
  
    tc1->cd(2);
    float maxy = 0.1;
    gNExposed[0]->GetYaxis()->SetRangeUser(0,maxy);
    gNExposed[0]->GetYaxis()->SetTitleOffset(yoffset);
    gNExposed[0]->GetYaxis()->SetTitle("Fraction of population");
    gNExposed[0]->GetXaxis()->SetTitleOffset(xoffset);
    gNExposed[0]->Draw("APL");
    for (int irun = 1; irun < filecounter; irun++) {	 gNExposed[irun]->Draw("samePL");	 }
  
    /*gNInfectious->GetYaxis()->SetTitleOffset(yoffset);
    gNInfectious->GetYaxis()->SetTitle("Fraction infectious");	
    gNInfectious->GetXaxis()->SetTitleOffset(xoffset);
    gNInfectious->SetLineColor(kRed+1);
    gNInfectious->SetMarkerColor(kRed+2);
    gNInfectious->Draw("samePL");	*/ 
  
    gNReported[0]->GetYaxis()->SetTitleOffset(yoffset);
    gNReported[0]->GetYaxis()->SetTitle("Fraction reported");	
    gNReported[0]->GetXaxis()->SetTitleOffset(xoffset);	
    for (int irun = 0; irun < filecounter; irun++) {	 gNReported[irun]->Draw("samePL");	 }


    TLegend *fl = new TLegend();
    fl->AddEntry(gNExposed[0], "Sick","PL");
    fl->AddEntry(gNReported[0], "Reported","PL");
    makeLegend(fl, gPad, 0)	;
    
    tc1->cd(3);	
    maxy = 0.5;
    gNQuarantineToday[0]->GetYaxis()->SetTitleOffset(yoffset);
    gNQuarantineToday[0]->GetYaxis()->SetTitle("Fraction quarantined");	
    gNQuarantineToday[0]->GetXaxis()->SetTitleOffset(xoffset);	
    gNQuarantineToday[0]->Draw("APL");
    for (int irun = 1; irun < filecounter; irun++) {	 gNQuarantineToday[irun]->Draw("samePL");	 }
    for (int irun = 0; irun < filecounter; irun++) {	 gNQuarantineNotInfectious[irun]->Draw("samePL");	 }
  
    TLegend *fl3 = new TLegend();
    fl3->AddEntry(gNQuarantineToday[0], "Quarantined","PL");				
    fl3->AddEntry(gNQuarantineNotInfectious[0], "Quarantined but not sick","PL");					
    makeLegend(fl3, gPad, 0)	;
  
    tc1->cd(4);
    TH1F *fhDaysToQuarantinePostIntervention = (TH1F*)fin->Get("fhDaysToQuarantinePostIntervention;1");
    TH1F *fhDaysToReportedPostIntervention = (TH1F*)fin->Get("fhDaysToReportedPostIntervention;1");
    TH1F *fhDaysToTestedPostIntervention = (TH1F*)fin->Get("fhDaysToTestedPostIntervention;1");
  
    if (fhDaysToQuarantinePostIntervention != NULL && fhDaysToReportedPostIntervention != NULL && fhDaysToTestedPostIntervention != NULL) {
      fhDaysToQuarantinePostIntervention->SetLineColor(kGreen+1);
      fhDaysToReportedPostIntervention->SetLineColor(kBlue+4);	
      fhDaysToTestedPostIntervention->SetLineColor(kOrange-3);	

      fhDaysToQuarantinePostIntervention->GetXaxis()->SetTitleOffset(2.3);
      fhDaysToQuarantinePostIntervention->GetYaxis()->SetTitleOffset(yoffset);
      fhDaysToQuarantinePostIntervention->GetYaxis()->SetTitle("Fraction of population");
      fhDaysToQuarantinePostIntervention->GetXaxis()->SetTitle("time difference [day]")	;
      fhDaysToQuarantinePostIntervention->Draw("");
      fhDaysToReportedPostIntervention->Draw("same");	
      fhDaysToTestedPostIntervention->Draw("same");	
      TLegend *fl2 = new TLegend();
  
      fl2->AddEntry(fhDaysToQuarantinePostIntervention, "(Q-I)","FL");
      fl2->AddEntry(fhDaysToReportedPostIntervention, "(R-E)","FL");	
      fl2->AddEntry(fhDaysToTestedPostIntervention, "(T-E)","FL");	
      makeLegend(fl2, gPad, 3)	;
      if (fDebug) cout << "Done histograms" << endl;
    }

    tc1->cd(5);

    hMaxDailySick->SetLineColor(kBlue+2);
    hMaxQuarantine->SetLineColor(kGreen+1);
    hMaxTotalSick->SetLineColor(kRed+1);
    hSick1yr->SetLineColor(kOrange+1);
    hQuarantineSum->SetLineColor(kGreen-2);
    hMaxDailySick->GetXaxis()->SetTitle("Fraction of population");
    hMaxDailySick->GetXaxis()->SetTitleOffset(2.3);	
    hMaxDailySick->DrawCopy();
    hMaxQuarantine->DrawCopy("same");
    hMaxTotalSick->DrawCopy("same");
    hSick1yr->DrawCopy("Same");
    hQuarantineSum->DrawCopy("same");
    TLegend *fl4 = new TLegend();
  
    fl4->AddEntry(hMaxDailySick, "Peak sick same day","FL");
    fl4->AddEntry(hMaxQuarantine, "Peak quarantine","FL");	
    //fl4->AddEntry(hMaxTotalSick, "Total exposed after 500 days","FL");	
    fl4->AddEntry(hSick1yr, "Total exposed after one year","FL");	
    fl4->AddEntry(hQuarantineSum, "Average fraction in quarantine","FL");	
    makeLegend(fl4, gPad, 0)	;
  
    tc1->cd(6);

    TPaveText *pt = new TPaveText(.05,.01,.5,.9);
    pt->SetTextAlign(31);
    pt->SetTextFont(42);
    pt->AddText(Form("nPersons = %d",fNPersons));
    pt->AddText(Form("peopleMetPerDay = %d",int(fPeopleMetPerDay)));
    pt->AddText(Form("transmissionProbability = %5.4f",transmissionProbability));   
    TText *t1 = 	pt->AddText(Form("tracingOrder = %d",fTracingOrder));                t1->SetTextFont(62);
    TText *t2 = 	pt->AddText(Form("AppProbability = %3.2f",fAppProbability));   		   t2->SetTextFont(62);
    TText *t2a = 	pt->AddText(Form("fTracingEfficiency = %3.2f",fTracingEfficiency)); t2a->SetTextFont(62);	
    TText *t3 = 	pt->AddText(Form("daysToTestResult = %d",fDaysToTestResult)); 				 t3->SetTextFont(62);
    TText *t4 = 	pt->AddText(Form("symptomProbability = %3.2f",symptomProbability));  t4->SetTextFont(62);
    pt->AddText(Form("dTTest = %d",fdTTest));		
    pt->AddText(Form("testPositiveProbability = %4.3f",testPositiveProbability));
    pt->AddText(Form("falsePositiveRate = %4.3f",falsePositiveRate));
    pt->Draw();	 
  
    TPaveText *pt2 = new TPaveText(.5,.1,.95,.8);
    pt2->SetTextAlign(11);
    pt2->AddText(Form("backwardTracing = %d",fBackwardTracing));
    pt2->AddText(Form("traceUninfected = %d",fTraceUninfected));
  
    pt2->AddText(Form("daysInQuarantine = %d",fDaysInQuarantine));
    pt2->AddText(Form("startTTIFraction = %6.5f",fStartTracingTestingInfectedFraction));
    pt2->AddText(Form("startTracingOnDay = %d",fStartTestingOnDay));
    pt2->AddText(Form("daysBackwardTrace = %d",fDaysBackwardTrace));
    pt2->AddText(Form("testThreshold = %4.3f",TestThreshold));
    pt2->AddText(Form("incubationGamma = %4.3f",incubationGamma)); 
    pt2->AddText(Form("incubationMu = %4.3f",incubationMu)); 
    pt2->AddText(Form("incubationBeta = %4.3f",incubationBeta)); 
    pt2->AddText(Form("infectionGamma = %4.3f",infectionGamma)); 
    pt2->AddText(Form("infectionMu = %4.3f",infectionMu)); 
    pt2->AddText(Form("infectionBeta = %4.3f",infectionBeta));    
    pt2->Draw();
    
    
    if (fDebug) {
      cout << "Saving canvas" << endl;
      tc1->SaveAs(Form("overview_%d.png",runner));
    }
 

   }
 if (fDebug) cout << "Opening output" << endl;
 if (!fWriteTextfile) {
 TFile *fout = new TFile(Form("Summary_%s",fname.c_str()),"RECREATE");
 TTree *outtree = new TTree("MCSummary","MCSummary tree");
 

  outtree->Branch("nPersons",&fNPersons);
  outtree->Branch("nDays",   &fNDays);
  outtree->Branch("fractionNoOutbreak", &fractionNoOutbreak); 
  outtree->Branch("fractionNoOutbreak50", &fractionNoOutbreak50);
  outtree->Branch("averagemaxdailysick", &averagemaxdailysick); 
  outtree->Branch("averagemaxdailysickWidth", &averagemaxdailysickWidth);    
  outtree->Branch("averagemaxq", &averagemaxq);  
  outtree->Branch("averagemaxqWidth", &averagemaxqWidth);    
  outtree->Branch("averagemaxtotalsick", &averagemaxtotalsick); 
  outtree->Branch("averagesickByDay365", &averagesickByDay365);
  outtree->Branch("averagesickByDay365Width", &averagesickByDay365Width); 
  outtree->Branch("averageQuarantineManDays", &averageQuarantineManDays); 
  outtree->Branch("averageQuarantineManDaysWidth", &averageQuarantineManDaysWidth); 
  outtree->Branch("averageDoublingTime", &averageDoublingTime);   
  outtree->Branch("averageReAfterIntervention", &averageReAfterIntervention);
  outtree->Branch("averageReAfterInterventionWidth", &averageReAfterInterventionWidth);
  outtree->Branch("averageReBeforeIntervention", &averageReBeforeIntervention);
  outtree->Branch("calculateReFromInterventionStartToDay",&calculateReFromInterventionStartToDay);
  outtree->Branch("numberSickToCountAsOutbreak",&numberSickToCountAsOutbreak);
  outtree->Branch("peopleMetPerDay", &fPeopleMetPerDay);
  outtree->Branch("daysInQuarantine",&fDaysInQuarantine);
  outtree->Branch("tracingOrder",    &fTracingOrder);
  outtree->Branch("startTracingTestingInfectedFraction", &fStartTracingTestingInfectedFraction);
  outtree->Branch("startTracingOnDay", &fStartTestingOnDay);
  outtree->Branch("daysBackwardTrace", &fDaysBackwardTrace);
  outtree->Branch("backwardTracing", &fBackwardTracing);
  outtree->Branch("traceUninfected", &fTraceUninfected);
  outtree->Branch("dTTest", &fdTTest);
  outtree->Branch("testThreshold", &TestThreshold);
  outtree->Branch("AppProbability", &fAppProbability);
	outtree->Branch("tracingEfficiency", &fTracingEfficiency);  
  outtree->Branch("daysToTestResult",  &fDaysToTestResult);
  outtree->Branch("randomTestingRate", &RandomTestingRate);
  outtree->Branch("randomTesting", &RandomTesting);
  outtree->Branch("symptomProbability",     &symptomProbability);   
  outtree->Branch("ReportingProbability", &fReportingProbability);    
  outtree->Branch("testPositiveProbability",&testPositiveProbability);
  outtree->Branch("falsePositiveRate",      &falsePositiveRate);
  outtree->Branch("transmissionProbability",&transmissionProbability);  
  outtree->Branch("R00",  &R00);
  outtree->Branch("R0e",  &R0e);  
  outtree->Branch("TracingDelay", &fTracingDelay);
  outtree->Branch("AsymptomaticTransmissionScaling",&AsymptomaticTransmissionScaling);
  outtree->Branch("incubationGamma",&incubationGamma); 
  outtree->Branch("incubationMu",   &incubationMu); 
  outtree->Branch("incubationBeta", &incubationBeta); 
  outtree->Branch("infectionGamma", &infectionGamma); 
  outtree->Branch("infectionMu",    &infectionMu); 
  outtree->Branch("infectionBeta",  &infectionBeta);
  outtree->Branch("SocialDistancingMaxPeople", &fSocialDistancingMaxPeople);
  outtree->Branch("SocialDistancingFrom", &fSocialDistancingFrom);
  outtree->Branch("SocialDistancingTo ", &fSocialDistancingTo);
  outtree->Branch("SocialDistancingFactor", &fSocialDistancingFactor);
    
  outtree->Fill();
  fout->cd();
  outtree->Write();
  if (fDebug) cout << "Closing output" << endl;
	fout->Close();
	}
	if (fWriteTextfile) out.close();
}

//! This is the main function. It finds all the root files in the input directory and runs the post processing on each one in turn.
void run() {
  string dirin = "results/";
	int runner = 0;
	TSystemDirectory dir("",dirin.c_str());

  TList *files = dir.GetListOfFiles();
  if (!files) { cout << "Fatal error: Cannot access input directory, or directory empty." << endl; return; }
	TSystemFile *file;
	TString fname;
	TIter next(files);
	TFile *ff = NULL;
	string fnameonly = "";
	TStopwatch stopwatch;
 	while ( (file=(TSystemFile*)next())) {
    fname =(file->GetName());
    string fnamestring = fname.Data(); 	
 		if (fnamestring.find(".root") != std::string::npos && runner < 5500) {
 		fnameonly = fnamestring;
 		cout << fnamestring.c_str() << endl;
 		fnamestring = dirin + fnamestring;
 				ff = new TFile(fnamestring.c_str());
 				postprocess(  ff , runner , fnameonly);
 				ff->Close();
 				runner++;
      }
    }
  cout << "Processed " << runner << " files" << endl;
  stopwatch.Print();
}


// because ROOT can't make a decent legend automatically
void makeLegend(TLegend *ll, TVirtualPad *pad, int location)
  {
    if (!pad || !ll) {
      cout << "Cannet make legend dimensions without pad or legend pointer. Please use makeDSLegend(legend, pad)." << endl;
      return;
    }
    pad->cd();

    // hardcoded constants that determine the size and spacing of legend entries
    int legendTextSize = gStyle->GetLabelSize() - 2; //px
    float legendSampleLength = 0.3; // length of sample lines as fraction of legend width
    int legendEntryPadding = 2; // vertical padding between entries
    int gapBetweenSampleAndLabel = 15; // px  horizontal distance between the sample line and the label

    // get some info about pad size and legend rows; this is used to figure out the dimensions of the legend box
    int padHeightinPx = gPad->VtoPixel(0);
    int padWidthinPx = gPad->UtoPixel(1);
    int legendRowsN = ll->GetNRows(); // how many entries are in this legend

    // utility variables to store legend data we need to determine legend dimensions
    TLegendEntry* entry;	// a pointer to each entry
    int longestLabel = 0;	// count how long the labels are
    string legendHeader = ""; // is there a header?
    gStyle->SetLegendFont(gStyle->GetLegendFont());
    ll->SetTextAlign(12);
    ll->SetTextSize(legendTextSize);
    ll->SetMargin(legendSampleLength); // length of sample lines as fraction of legend width

    // I guess legends are exempt from the general gStyle rules ...
    ll->SetBorderSize(0);
    ll->SetFillColor(kWhite);
    ll->SetFillStyle(1001);

    // determine height
    float legendHeight = float((legendRowsN * (legendTextSize + 2*legendEntryPadding)))/float(padHeightinPx);
    // determine width based on number of characters in labels

    if (ll->GetHeader()) legendHeader = ll->GetHeader();
    int iLegendRow = (legendHeader.size() < 2 ? 0 : 1) ;  // if there is a header, we skip it in the loop
    for (; iLegendRow < ll->GetNRows(); iLegendRow++) {  // i = 0:  header, if there is one
      entry = (TLegendEntry*)ll->GetListOfPrimitives()->At(iLegendRow);
      string labelText = entry->GetLabel();
      if ((int)labelText.size() > longestLabel) {
        longestLabel = labelText.size();
      }
    }

    // determine width

    // average width of text = (number of chars) * (size in px)/(font constant)
    // font constant = 2.24 for Helvetica
    float legendTextWidth = longestLabel * legendTextSize/2.24;
    float legendWidth = legendTextWidth * (1 + legendSampleLength) + gapBetweenSampleAndLabel;
    float titleWidth = legendHeader.size()* legendTextSize/2.24;
    if (legendWidth < titleWidth) { // we cut off the title if the box isn't big enough
      legendSampleLength = (titleWidth - legendTextWidth - gapBetweenSampleAndLabel)/titleWidth; // make samples longer so they aren't so far away from their labels
      ll->SetMargin(legendSampleLength);
      legendWidth = titleWidth;
    }

    legendWidth = float(legendWidth)/float(padWidthinPx);

    /*** legend location ***/
    switch (location) {
      case 0: 	// top right corner aligned with histogram top right corner
        ll->SetY2NDC(1. - gPad->GetTopMargin() - 0.02);
        ll->SetX2NDC(1 - gPad->GetRightMargin() - 0.02);
        ll->SetY1NDC(ll->GetY2NDC() - legendHeight);
        ll->SetX1NDC(ll->GetX2NDC() - legendWidth);
        break;
      case 1: // bottom right
        ll->SetX2NDC(1 - gPad->GetRightMargin() - 0.02);
        ll->SetY1NDC(gPad->GetBottomMargin() + 0.02);

        ll->SetY2NDC(ll->GetY1NDC() + legendHeight);
        ll->SetX1NDC(ll->GetX2NDC() - legendWidth);
        break;

      case 2: // top left
        ll->SetY2NDC(1. - gPad->GetTopMargin() - 0.02);
        ll->SetX1NDC(gPad->GetLeftMargin() + 0.02);
        ll->SetX2NDC(ll->GetX1NDC() + legendWidth);
        ll->SetY1NDC(ll->GetY2NDC() - legendHeight);
        break;
      case 3: // bottom left
        ll->SetY1NDC(gPad->GetBottomMargin() + 0.02);
        ll->SetX1NDC(gPad->GetLeftMargin() + 0.02);
        ll->SetX2NDC(ll->GetX1NDC() + legendWidth);
        ll->SetY2NDC(ll->GetY1NDC() + legendHeight);
        break;
    }
    ll->Draw();
  }	
