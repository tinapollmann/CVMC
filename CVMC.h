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
 * Simulate an outbreak of a virus.
 * 
 * @author Tina Pollmann, Christoph Wiesinger
 * @date 2020
 */
 

#ifndef CVMC_H
#define CVMC_H


//c++
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>

//root
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TPad.h>
#include "CVDisease.h"
#include "CVPerson.h"

using namespace std;

class CVMC 
{
  public:
    //! create MC outbreak with .. people for .. days
    CVMC(int nPersons,int nDays,float appProb=0.5, float repProb = 0.5, string outputPrefix="CovidMCResult") 
    : fNPersons(nPersons),fNDays(nDays),fAppProbability(appProb),fReportingProbability(repProb),fOutputPrefix(outputPrefix) {
      //! general
      fMaxPeopleInDotFile = 400;
      fDebug = false;      
      //! default parameters
      fPeopleMetPerDay           =  11;
      fSocialDistancingMaxPeople =  10;
      fSocialDistancingFrom      =  20;
      fSocialDistancingTo        =  50;
      fSocialDistancingFactor    = 0.5;
      fDaysInQuarantine          =  14;
      //! .. for tracing
      fTracingOrder      =     2;      
      fStartTracingOnDay = 99999; 
      fDaysBackwardTrace =    10; 
      fBackwardTracing   =  true; 
      fTraceUninfected   =  true;
      fTracingEfficiency =    1.;
      fTracingDelay      =     0;
      //! .. for testing
      fStartTestingOnDay = 99999;      
      fDaysToTestResult  =     0;    
      fdTTest            =     5; 
      fRandomTesting     = false;
      fRandomTestingRate =  0.01;
      //! This gamma distribution, given the 'mean people per day' as par0 from the paper https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074
      //! approximatly gives the distributions as observed in the paper.
      fPeopleMetFunction = new TF1("fPeopleMetFunction","x>[1]?TMath::GammaDist(x,2.,0.,[0]/2.):0.", 0., 60.);
      fPeopleMetFunction->SetParameter(0,fPeopleMetPerDay);
      fPeopleMetFunctionDistancing = (TF1*) fPeopleMetFunction->Clone("fPeopleMetFunctionDistancing");
      fPeopleMetFunctionDistancing->SetParameter(0,fPeopleMetPerDay*fSocialDistancingFactor);
      fPeopleMetFunctionDistancing->SetRange(0,fSocialDistancingMaxPeople);      
      //! make the people 
      for(int i=0;i<fNPersons;i++) 
        fPersons.push_back(new CVPerson(fNDays,&fRandom));
      if(fDebug) cout << "* made " << fPersons.size() << " people" << endl;  
      //! ... and the disease
      fDisease = new CVDisease(&fRandom); 
      //! vector to store qurantine stats   
      fNQuarantine = vector<int>(fNDays,0);
      //! create histograms
      fhIncubationPeriod = new TH1F("fhIncubationPeriod",";Incubation period [day]",30,0,30);
      fhInfectiousness   = new TH1F("fhInfectiousness",";Days since symptom onset [day]; Infectiousness ",24,-4,20);
      fhLatentPeriod     = new TH1F("fhLatentPeriod",";Latend period [day]",15,0,15);
      fhInfectiousnessAgeSymptom        = new TH1F("fhInfectiousnessAgeSymptom",";Age of infection [day]",40,0,40);
      fhInfectiousnessAgeNoSymptom      = new TH1F("fhInfectiousnessAgeNoSymptom",";Age of infection [day]",40,0,40);
      fhInfectiousnessDurationSymptom   = new TH1F("fhInfectiousnessDurationSymptom","; duration from first to last person infected [day]",30,0,30);
      fhInfectiousnessDurationNoSymptom = new TH1F("fhInfectiousnessDurationNoSymptom","; duration from first to last person infected [day]",30,0,30);
      fhDaysToQuarantinePostIntervention= new TH1F("fhDaysToQuarantinePostIntervention","; (day quarantined - day infectious) [day]",40,-20,20);
      fhDaysToQuarantinePreIntervention = new TH1F("fhDaysToQuarantinePreIntervention","; (day quarantined - day infectious) [day]",40,-20,20);        
      fhDaysToReportedPostIntervention  = new TH1F("fhDaysToReportedPostIntervention","; (day reported - day exposed) [day]",32,-2,30);
      fhDaysToReportedPreIntervention    = new TH1F("fhDaysToReportedPreIntervention","; (day reported - day exposed) [day]",32,-2,30);
      fhDaysToTestedPostIntervention    = new TH1F("fhDaysToTestedPostIntervention","; (day tested - day exposed) [day]",32,-2,30);
      fhDaysToTestedPreIntervention      = new TH1F("fhDaysToTestedPreIntervention","; (day tested - day exposed) [day]",32,-2,30);    
      fhNumberInfectedPreIntervention      = new TH1F("fhNumberInfectedPreIntervention","; number of infectees",60,0,60);    
          
      //! link output tree branches
      fPopulationLevelInformation = new TTree("PopulationLevelInformation","CV MC");
      fPopulationLevelInformation->Branch("fileIndex",&fRunId);            
      fPopulationLevelInformation->Branch("day",&fDayForROOTTree);
      fPopulationLevelInformation->Branch("fNSusceptible",&fNSusceptible);
      fPopulationLevelInformation->Branch("fNExposed",&fNExposed);
      fPopulationLevelInformation->Branch("fNInfectious",&fNInfectious);
      fPopulationLevelInformation->Branch("fNRecovered",&fNRecovered);
      fPopulationLevelInformation->Branch("fNTraced",&fNTraced);
      fPopulationLevelInformation->Branch("fNReported",&fNReported);
      fPopulationLevelInformation->Branch("feffectiveR",&feffectiveR);
      fPopulationLevelInformation->Branch("feffectiveRUncertainty",&feffectiveRUncertainty);
      fPopulationLevelInformation->Branch("fNQuarantineToday",&fNQuarantineToday);  
    }
    ~CVMC() {
      //! delete the people
      for(auto kv : fPersons) delete kv;
      //! ... and the disease
      delete fDisease;  
    }
      
    //! getters
    int    GetNPersons()               { return fNPersons;       }
    int    GetNDays()                  { return fNDays;          }
    float  GetAppProbability()         { return fAppProbability; }
    float  GetReportingProbability()   { return fReportingProbability; }     
    string GetOutputPrefix()           { return fOutputPrefix;   }   
    
    int    GetDaysInQuarantine()            { return fDaysInQuarantine;            }
    float  GetPeopleMetPerDay()             { return fPeopleMetPerDay;             }
    float  GetSocialDistancingMaxPeople()   { return fSocialDistancingMaxPeople;   }
    int    GetSocialDistancingFrom()        { return fSocialDistancingFrom;        }
    int    GetSocialDistancingTo()          { return fSocialDistancingTo;          }   
    float  GetSocialDistancingFactor()      { return fSocialDistancingFactor;      }
    TF1*   GetPeopleMetFunction()           { return fPeopleMetFunction;           }
    TF1*   GetPeopleMetFunctionDistancing() { return fPeopleMetFunctionDistancing; }    
    int    GetPeopleMetToday(int day=0)     {
      if(day >= fSocialDistancingFrom && day < fSocialDistancingTo) return (int) fPeopleMetFunctionDistancing->GetRandom(); 
      return (int) fPeopleMetFunction->GetRandom();
    }
    
    int    GetTracingOrder()           { return fTracingOrder;      }       
    float  GetStartTracingTestingInfectedFraction() { return fStartTracingTestingInfectedFraction;  } 
    int    GetDaysBackwardTrace()      { return fDaysBackwardTrace; } 
    bool   GetBackwardTracing()        { return fBackwardTracing;   } 
    bool   GetTraceUninfected()        { return fTraceUninfected;   } 
    float  GetTracingEfficiency()      { return fTracingEfficiency; } 
    int    GetTracingDelay()           { return fTracingDelay;      }

    int    GetStartTestingOnDay()      { return fStartTestingOnDay; }
    int    GetDaysToTestResult()       { return fDaysToTestResult;  } 
    bool   GetRandomTesting()          { return fRandomTesting;     }
    float  GetRandomTestingRate()      { return fRandomTestingRate; }
    
    CVDisease*        GetDisease()     { return fDisease; }
    vector<CVPerson*> GetPersons()     { return fPersons; }   
    
    //! setters
    void SetPeopleMetPerDay(int peopleMetPerDay) { 
      fPeopleMetPerDay = peopleMetPerDay; 
      fPeopleMetFunction->SetParameter(0,fPeopleMetPerDay); 
      fPeopleMetFunctionDistancing->SetParameter(0,fPeopleMetPerDay*fSocialDistancingFactor); 
    }
    void SetSocialDistancingMaxPeople(float socialDistancingMaxPeople) { 
      fSocialDistancingMaxPeople = socialDistancingMaxPeople;
      fPeopleMetFunctionDistancing->SetRange(0,fSocialDistancingMaxPeople);
    } 
    void SetSocialDistancingFactor(float socialDistancingFactor) { 
      fSocialDistancingFactor = socialDistancingFactor;  
      fPeopleMetFunctionDistancing->SetParameter(0,fPeopleMetPerDay*fSocialDistancingFactor);     
    }
    void SetSocialDistancingFrom(int socialDistancingFrom) { fSocialDistancingFrom=socialDistancingFrom; cout << "Error, fixed start of social distancing is currently not supported." << endl; }
    void SetSocialDistancingTo(int socialDistancingTo)     { fSocialDistancingTo=socialDistancingTo;     }  
    
    void SetDaysInQuarantine(int daysInQuarantine)     { fDaysInQuarantine = daysInQuarantine;     }    
    void SetTracingOrder(int tracingOrder)             { fTracingOrder = tracingOrder;             }       
    void SetStartTracingTestingInfectedFraction(float startTracingOnFrac) { fStartTracingTestingInfectedFraction = startTracingOnFrac;   } 
    void SetDaysBackwardTrace(int daysBackwardTrace)   { fDaysBackwardTrace = daysBackwardTrace;   } 
    void SetBackwardTracing(int backwardTracing)       { fBackwardTracing = backwardTracing;       }
    void SetTraceUninfected(int traceUninfected)       { fTraceUninfected = traceUninfected;       }
    void SetTracingEfficiency(float tracingEfficiency) { fTracingEfficiency = tracingEfficiency;   }
    void SetTracingDelay(int tracingDelay)             { fTracingDelay = tracingDelay; };
        
    void SetStartTestingOnDay(int startTestingOnDay)   { fStartTestingOnDay = startTestingOnDay; cout << "Error, fixed day for testing start is currently not supported." << endl;  }
    void SetDaysToTestResult(int daysToTestResult)     { fDaysToTestResult = daysToTestResult;     }
    void SetDTTest(int dTTest)                         { fdTTest = dTTest; };
    void SetRandomTesting(bool randomTesting)          { fRandomTesting = randomTesting;     }
    void SetRandomTestingRate(float randomTestingRate) { fRandomTestingRate = randomTestingRate; }
    
    void SetDebug(bool debug=true)                     { fDebug=debug;                           }
    void SetMaxPeopleInDotFile(int maxPeopleInDotFile) { fMaxPeopleInDotFile=maxPeopleInDotFile; }
    
    //! run mc with id .. and seed ..
    void Run(int runId=0,int seed=0) {
      fRunId = runId;
      fRandom.SetSeed(seed);
      fRandomSeed = fRandom.GetSeed();
      gRandom->SetSeed(fRandom.Integer(4294967295)); //! for TF1::GetRandom
      //! reset all counters etc.
      Reset();
      //! .. and run it
      DoMC();
    }      
  
  protected:
    void Reset() {
      //! reset the population
      for(auto kv : fPersons) {
        kv->Reset(); 
        //! ... and roll out the app
        if(fRandom.Uniform()<fAppProbability) kv->SetHasApp();
        if(fRandom.Uniform()<fReportingProbability) kv->SetDoesReport(); 
      }
      //! reset population statistics
      fNExposed    = 0;
      fNInfectious = 0;
      fNRecovered  = 0;
      fNTraced     = 0;
      fNReported   = 0;   
      fStartTracingOnDay = 99999;
      fStartTestingOnDay = 99999;
      fSocialDistancingFrom = 99999;
      fill(fNQuarantine.begin(), fNQuarantine.end(), 0);      
      fNumberRecoveredByDay.clear();
      fNumberInfectiousByDay.clear();
      numberNewlyInfectedByDay.clear();      
      //! reset output tree
      fPopulationLevelInformation->Reset();       
      //! reset histograms  
      fhInfectiousnessAgeSymptom->Reset();
      fhInfectiousnessAgeNoSymptom->Reset();
      fhInfectiousnessDurationNoSymptom->Reset();
      fhInfectiousnessDurationSymptom->Reset();
      fhLatentPeriod->Reset();
      fhInfectiousness->Reset();
      fhIncubationPeriod->Reset();
      fhDaysToQuarantinePostIntervention->Reset();
      fhDaysToQuarantinePreIntervention->Reset();
      fhDaysToReportedPostIntervention->Reset();
      fhDaysToReportedPreIntervention->Reset();
      fhDaysToTestedPostIntervention->Reset();
      fhDaysToTestedPreIntervention->Reset();
      fhNumberInfectedPreIntervention->Reset();
      //! 
      fLastDayWithPatients   = 0;
      fLastDayForDot          = 0;
      fLastDayForDotSimple   = 0;
      fPeopleInDotFile       = 0; 
      fDotStringPeople       = "";
      fDotStringPeopleSimple = "";
      fDotStringLegend       = "";
      //!     
      fDayForROOTTree        = 0;
      feffectiveR            = 0.0;
      feffectiveRUncertainty = 0.0;
      fNQuarantineToday      = 0;     
    }
    
    void DoMC();
    void Trace(CVPerson* aperson,int day,int fromPersonID, int direction);
    bool DoTrace(CVPerson* aperson,int day,int fromPersonID, int tlevel, int direction);
    void TraceUninfected(CVPerson* aperson,int day);
    void TraceBack(int id,int day,int fromPersonID);    
    void AddPersonToDot(CVPerson *kv);
    void AddConnectionToDot(int IDfrom, int IDto, int tlevel=0, int tday=-1);
    void AddDaysToDot() ;
    void AddLegendToDot();
    void AddHeaderToTSV();   
    void WriteGnuplotScript(); 
    void IncrementPopulationStatistics(CVInfectionStatus is, CVTracingStatus ts);
    void AddDayToTSVAndROOT(int day);
    void WriteROOTFile();
    void makeLegend(TLegend *ll, TPad *pad, int location=0); // I need this to stay sane
    void FillDiagnostics(CVPerson *kv, bool posttracing);
    
  private:
    TRandom3 fRandom = TRandom3(0);
    unsigned int fRandomSeed;
    bool  fDebug;
    int   fRunId = -1;
    
    int   fNPersons;
    int   fNDays;
    
    float fStartTracingTestingInfectedFraction; //! At what fraction of cummulative reported people tracing and testing starts
    float fPeopleMetPerDay;    //! how many people one meets per day    
    float fSocialDistancingMaxPeople;
    int   fSocialDistancingFrom;
    int   fSocialDistancingTo;   
    float fSocialDistancingFactor;
    int   fDaysInQuarantine;   //! when a person gets traced, how long are they in quarantine for?
    int   fdTTest;             //! How many days before people in quarantine get tested
    int   fTracingOrder;       //! do we trace contacts (1st order) or contacts of contacts (2nd order)    
    int   fStartTracingOnDay;  //! after how many days is contact tracing initiated
    int   fDaysBackwardTrace;  //! when we trace, how many days do we go back into the past of a person and find their contacts?
    bool  fTraceUninfected;    //! Do we trace uninfected contacts?
    bool  fBackwardTracing;    //! switch on/off backward tracing
    float fAppProbability;     //! probability that a person who got infected by a known carrier will be traced
    float fReportingProbability;//! probability that a person who has symptoms goes to the doctor and is reported
    float fTracingEfficiency;  //! efficiency that a contact is properly recorded
    int   fTracingDelay; 

    int   fStartTestingOnDay;   //! after how many days do we start testing traced individuals
    int   fDaysToTestResult;    //! how many days to get the result of the test (or how many days before people get tested)
    bool  fRandomTesting;
    float fRandomTestingRate;
    
    TF1*  fPeopleMetFunction;           //! Distribution of how many others one meets in a day   
    TF1*  fPeopleMetFunctionDistancing; //! .. under social distancing 
    CVDisease*        fDisease; //! this is the disease
    vector<CVPerson*> fPersons; //! this will hold every person in the population 
    
    //! Output and gathering of statistics 
    string fOutputPrefix;
    ofstream fOutputDot;       //! Graphviz dot format
    ofstream fOutputDotSimple; //! Graphviz dot format without numbers in the plot
    int fPeopleInDotFile;      //! How many people we have tracked in the dot file so far
    int fMaxPeopleInDotFile;   //! How many we track at most
    int fLastDayWithPatients;  //! The last day in the simulation where people are still infected    
    int fLastDayForDot;        //! where dot output stops or we get too much
    int fLastDayForDotSimple;
    string fDotStringPeople;
    string fDotStringPeopleSimple;    
    string fDotStringLegend;
    ofstream fOutputGnuplot;
    
    ofstream fOutputTSV;   //! Tab separated file
    TFile *fOutputRoot;    //! Root output file    
    //! diagnostic output histograms
    TH1F *fhIncubationPeriod;
    TH1F *fhInfectiousness;
    TH1F *fhInfectiousnessAgeNoSymptom;
    TH1F *fhInfectiousnessAgeSymptom;
    TH1F *fhInfectiousnessDurationSymptom;
    TH1F *fhInfectiousnessDurationNoSymptom;
    TH1F *fhLatentPeriod;
    TH1F *fhDaysToQuarantinePostIntervention;
    TH1F *fhDaysToQuarantinePreIntervention;
    TH1F *fhDaysToReportedPostIntervention;
    TH1F *fhDaysToReportedPreIntervention;    
    TH1F *fhDaysToTestedPostIntervention;
    TH1F *fhDaysToTestedPreIntervention;
    TH1F *fhNumberInfectedPreIntervention;
        
    TTree *fPopulationLevelInformation;
    // R_e and doubling time
    vector<int> fNumberRecoveredByDay;  //! number of people who recovered on each day
    vector<int> fNumberInfectiousByDay; //! number of people who are infectioues on each day

    vector<int> numberNewlyInfectedByDay; //! number of people who got newly infected on each day
    int fDayForROOTTree;
    int fNSusceptible;
    int fNExposed;
    int fNExposedTotal;
    int fNInfectious;
    int fNRecovered;
    int fNTraced;
    int fNReported;  
    float feffectiveR;
    float feffectiveRUncertainty;
    int fNQuarantineToday;
    vector<int> fNQuarantine;    
};

#endif
