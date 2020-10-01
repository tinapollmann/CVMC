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
 * Container for disease properties
 * 
 * @author Tina Pollmann, Christoph Wiesinger
 * @date 2020
 */


#ifndef CVDisease_H
#define CVDisease_H

//! c++
#include <vector>

//! root
#include <TRandom3.h>
#include <TRandom.h>
#include <TF1.h>
#include <TMath.h>
using namespace std;

//! the state of the infection a person is in
enum CVInfectionStatus : uint8_t {
  S_Susceptible,       ///< 0 - never been infected
  S_Exposed,           ///< 1 - exposed, not yet infectious or symptomatic
  S_ExposedInfectious, ///< 2 - infectious, but not symptomatic
  S_Infectious,        ///< 3 - infectious; symptomatic with probability p_symptom
  S_Recovered,         ///< 4
  S_Invalid,           ///< 5 
  Infection_COUNT      ///< Do not add anything after COUNT
};

//! base covid19
class CVDisease
{
  public:
    //! construct and pass global random number generator
    CVDisease(TRandom3* random) : fRandom(random) {
      //! keep track if intialized with internal generator
      fInternalRandom = false;
      //! create pdfs
      fIncubationPeriod = new TF1("fIncubationPeriod","x>[1]?TMath::GammaDist(x,[0],[1],[2]):0",     0.,31.);
      fInfectiousness   = new TF1("fInfectiousness",  "x>[1]?TMath::GammaDist(x,[0],[1],[2])*[3]:0",-5.,31.);
      //! set defaults
      SetSymptomProbability(0.14);
      SetTestPositiveProbability(0.95);
      SetFalsePositiveRate(0.01);
      SetIncubationParameters(3.93,-0.99,2.15);
      fTransmissionProbability = 0.064;
      SetAsymptomaticTransmissionScaling(0.1);
      SetInfectiousnessParameters(2.08,-2.42,1.56);   
      fTestThreshold=0.01;            
    }
    //! contruct with internal random number generator
      CVDisease() : CVDisease(new TRandom3(0)) {
      //! .. and keep track that intialized with internal generator
      fInternalRandom = true;
    }
    //! destruct
    ~CVDisease() {
      //! .. also the internal random number generator
      if(fInternalRandom) delete fRandom;
      delete fIncubationPeriod;
      delete fInfectiousness;
    }
    
    //! getters
    float GetSymptomProbability()      { return fSymptomProbability;      }
    float GetTestPositiveProbability() { return fTestPositiveProbability; }
    float GetFalsePositiveRate()       { return fFalsePositiveRate;       }
    float GetTransmissionProbability() { return fTransmissionProbability; }
    float GetTestThreshold()           { return fTestThreshold;           }
    float GetAsymptomaticTransmissionScaling() { return fAsymptomaticTransmissionScaling; }
    
    TF1*  GetIncubationPeriod() { return fIncubationPeriod; }
    TF1*  GetInfectiousness()   { return fInfectiousness;   }
    
    float GetAsymptomaticTime() { return fAsymptomaticTime; } //! days
    int   GetDaysForRecovery()  { return fDaysForRecovery;  } //! days
    
    //! special getters
    float GetInfectiousness(int day,bool symptomatic=true) { 
      //! evaluate infectiousness at day since symptom onset ...
      float infectiousness = 0.0;
      //! .. but the curve starts at -2.5 days and we do not simulate half days. However, ignoring the last 0.5 days of the curve 
      // chates the normailzation. So we sample randomly at this time.
      if (day < 0) infectiousness = fInfectiousness->Eval(fRandom->Uniform(day,day+1));
      //! .. at positive times there is no problem and we just find the infectivity for the day
      else infectiousness = fInfectiousness->Eval(day);
      //! .. and scale it for asymptomatic and presymptomatic cases
      if(!symptomatic || day <= 0) infectiousness *= fAsymptomaticTransmissionScaling;
      return infectiousness;
    }
    bool GetTestsPositive(int day) { //! day is day since symptom onset
      if(fRandom->Uniform()<GetTestPositiveProb(day)) return true;
      else return false;
    }
    float GetTestPositiveProb(int day) {
      //! test efficiency depends on infectiousness
      float infectiousness = GetInfectiousness(day);
      //! .. above threshold the test efficiency stays constant 
      if(infectiousness > fTestThreshold) return fTestPositiveProbability;
      //! .. before the onset it is zero, except where the infectiousness is still above threshold
      if(day<0) return 0;
      //! .. afterwards it decays with the slope of the infectiousness
      float prob = fTestPositiveProbability/fTestThreshold * infectiousness;
      return prob;
    }
      
    //! setters
    void SetSymptomProbability(float prob)      { fSymptomProbability      = prob; }
    void SetTestPositiveProbability(float prob) { fTestPositiveProbability = prob; }
    void SetFalsePositiveRate(float prob)       { fFalsePositiveRate       = prob; }
    void SetTestThreshold(float thresh)         { fTestThreshold = thresh; }
    void SetTransmissionProbability(float prob) { 
      fTransmissionProbability = prob; 
      //! reset infectiousness curve
      SetInfectiousnessParameters(
        fInfectiousness->GetParameter(0),
        fInfectiousness->GetParameter(1),
        fInfectiousness->GetParameter(2));
    }
    void SetAsymptomaticTransmissionScaling(float scale) { fAsymptomaticTransmissionScaling=scale; }
    void SetIncubationParameters(float gamma, float mu, float beta) {
      fIncubationPeriod->SetParameters(gamma,mu,beta);
    }
    void SetInfectiousnessParameters(float gamma, float mu, float beta) {
      if(mu<-5) cout << "warning: infection onset too early" << endl;
      //! start with normalized pdf
      fInfectiousness->SetParameters(gamma,mu,beta,1.0);
      //! infectioussness starts curve onset
      fAsymptomaticTime = (float) -mu;
      //! infectiousness ends at the 99% quantile 
      double quantile;
      double probSum = 0.99;
      fInfectiousness->GetQuantiles(1,&quantile,&probSum);
      fDaysForRecovery = TMath::Nint(quantile);
      //! update normalization
      fInfectiousness->SetParameter(3,(quantile-mu)*fTransmissionProbability);
    } 
    
    //! generate random course of infection
    vector<CVInfectionStatus> DrawCourse() {
      vector<CVInfectionStatus> course;
      //! draw day of symptom onset
      int dayOfSymptomOnset = (int) fIncubationPeriod->GetRandom();
      //! calculate backward from this day
      int latentDays = max(0,(int)(dayOfSymptomOnset-fAsymptomaticTime));  
      int asymptomaticInfectiousDays = dayOfSymptomOnset-latentDays;  
      //! fill course
      course.insert(course.end(),latentDays,S_Exposed);
      course.insert(course.end(),asymptomaticInfectiousDays,S_ExposedInfectious);
      course.insert(course.end(),fDaysForRecovery,S_Infectious);
      return course;       
    }
    //! draw if infection will develop symptomatic
    bool DrawHasSymptoms() { 
      if(fRandom->Uniform()<fSymptomProbability) return true;
      return false;
    }
   
  private:
    TRandom3* fRandom;
    bool      fInternalRandom; 
    
    //! probabilities/parameters
    float fSymptomProbability;      //! probability that an infected person will have symptoms and go to the doctor
    float fTestPositiveProbability; //! probability that an infected person when tested tests positive
    float fFalsePositiveRate;       //! probability to get a false positive when tested
    float fTransmissionProbability; //! probability to pass the infection on per contact per day
    float fTestThreshold;           //! Threshold infection probability per day when test starts to lose sensitivity to presence of viral RNA
    float fAsymptomaticTransmissionScaling;
    //! pdfs
    TF1* fIncubationPeriod;
    TF1* fInfectiousness;
    //! ... derived quantities
    float fAsymptomaticTime; //! days
    int   fDaysForRecovery;  //! days
};

#endif
