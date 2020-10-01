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
 * Defines a person in the population in which COVID-19 spread is simulated.
 * 
 * @author Tina Pollmann, Christoph Wiesinger
 * @date 2020
 */
 
#ifndef CVPerson_H
#define CVPerson_H

//! c++
#include <vector>
#include <algorithm>

#include "CVDisease.h"

using namespace std;

//! the state of tracing or reporting a person is in
enum CVTracingStatus : uint8_t {
  S_Unreported,     ///< 0 
  S_Reported,       ///< 1
  S_Traced,         ///< 2
  S_InvalidTracing, ///< 3
  Tracing_COUNT     ///< Do not add anything after COUNT
};

class CVPerson 
{
  public:
    //! construct and pass global random number generator
    CVPerson(int days,TRandom3* random) : fRandom(random) {
      //! keep track if intialized with internal generator
      fInternalRandom = false;
      //! automatic unique id
      fId = nextId++;
      //! max lifetime
      fNDays = days;
      //! status containers
      fInfectionStatus  = vector<CVInfectionStatus>(days,S_Susceptible);
      fTracingStatus    = vector<CVTracingStatus>(days,S_Unreported);
      fQuarantineStatus = vector<bool>(days,false);
      //! reset the persons memory
      Reset(); //! -> happens anyhow before MC run
    }
    //! construct with internal random number generator
    CVPerson(int days=100) 
      : CVPerson(days,new TRandom3(0)) {
      //! .. and keep track that intialized with internal generator
      fInternalRandom = true;
    }
    //! destruct
    ~CVPerson() {
      //! .. also the internal random number generator
      if(fInternalRandom) delete fRandom; 
    }
    
    void Reset() {
      //! reset memorized events  
      fExposedOn       = -1;  
      fInfectiousOn    = -1; 
      fSymptomOnset    = -1;
      fReportedOn      = -1;  
      fRecoveredOn     = -1; 
      fDayLastTestedOn = -1;
      fQuarantinedOn.clear();
      fTracedOn.clear();  
      //! reset status containers
      fill(fInfectionStatus.begin(), fInfectionStatus.end(), S_Susceptible);
      fill(fTracingStatus.begin(),   fTracingStatus.end(),   S_Unreported );
      fill(fQuarantineStatus.begin(),fQuarantineStatus.end(),false        ); 
      //! reset if symptomatic / has app
      fHasSymptoms = false; //! <- set at exposure
      fHasApp      = false; //! <- app will be given to the person
      fDoesReport  = false; //! <- given when object created
      //! reset contact memory
      fInfectedBy = -1;
      fExposed.clear();
      fDays.clear();
      fUninfectedContacts.clear();
    };
   
    //! getters   
    int        GetId()      { return fId;      }     
    CVDisease* GetDisease() { return fDisease; }
    
    int GetExposedOn()       { return fExposedOn;       }
    int GetInfectiousOn()    { return fInfectiousOn;    }
    int GetSymptomOnset()    { return fSymptomOnset;    }
    int GetReportedOn()      { return fReportedOn;      }
    int GetRecoveredOn()     { return fRecoveredOn;     }    
    int GetDayLastTestedOn() { return fDayLastTestedOn; }
    int GetQuarantinedOn(int tindex) { 
      try{ return fQuarantinedOn.at(tindex);} 
      catch(...) { return -1; } 
    }  
    int GetNQuarantinedOn()  { return fQuarantinedOn.size(); }
    int GetTracedOn(int tindex) {  
      try{ return fTracedOn.at(tindex); } 
      catch(...) { return -1; } 
    }     
    int GetNTracedOn()       { return fTracedOn.size(); }
    int GetNTracedBetween(int startday, int endday) { 
      int ntraced = 0; 
      for(auto kv : fTracedOn) { 
        if(kv >= startday && kv <= endday) ntraced++; 
      } 
      return ntraced; 
    }     
    
    vector<CVInfectionStatus> GetInfectionStatus()  { return fInfectionStatus;  }
    vector<CVTracingStatus>   GetTracingStatus()    { return fTracingStatus;    }
    vector<bool>              GetQuarantineStatus() { return fQuarantineStatus; }
    CVInfectionStatus GetInfectionStatus(int day) { 
      if(day < fNDays) { return fInfectionStatus.at(day); } 
      return S_Invalid ;  
    }
    CVTracingStatus GetTracingStatus(int day) { 
      if(day < fNDays) {return fTracingStatus.at(day); } 
      return S_InvalidTracing; 
    }
    bool GetQuarantineStatus(int day) { 
      if(day < fNDays) { return fQuarantineStatus.at(day); } 
      return S_InvalidTracing; 
    }    
    bool GetHasSymptoms()  { return fHasSymptoms; }  
    bool GetHasApp()       { return fHasApp;      } 
    bool GetDoesReport()   { return fDoesReport;      } 
    
    int          GetInfectedBy()                { return fInfectedBy;         }    
    vector<int>  GetExposed()                   { return fExposed;            }
    int          GetNExposed()                  { return fExposed.size();     }
    int          GetExposedPersonID(int pindex) { return fExposed.at(pindex); }    
    vector<int>  GetDays()                      { return fDays;               }
    int          GetNExposedDays()              { return fDays.size();        }
    int          GetExposedDays(int pindex)     { return fDays.at(pindex);    }  
    map<int,int> GetUninfectedContatcts()       { return fUninfectedContacts; }    
    int          GetNUninfectedContacts()       { return fUninfectedContacts.size(); }
    int          GetUninfectedContactsDay(int pindex) { 
      if(!fUninfectedContacts.count(pindex)) return -99;
      return fUninfectedContacts.at(pindex);
    }
    
    //! special getters
    bool GetTestsPositive(int day) { //! this is day into the outbreak
      //! record testing
      SetDayLastTestedOn(day);
      //! we have never been exposed
      if(GetExposedOn() < 0){
        //! .. and might get a false postive result
        if(fRandom->Uniform()<fDisease->GetFalsePositiveRate()) return true; 
        return false;
      }
      //! we recovered already
      if(GetRecoveredOn() < day && GetRecoveredOn() > 0)  {
        //! .. and if we know that we do not test positive
        if(GetReportedOn() < day) return false;
        //! .. but if we dont know that
        else {
          //! .. a false postive result might show up
          if(fRandom->Uniform()<fDisease->GetFalsePositiveRate()) return true; 
          else return false;
      } }
      //! we are still sick
      return fDisease->GetTestsPositive(day-GetSymptomOnset()); //! relative to symptom onset      
    }
    
    //! setters
    void SetExposedOn(int day)       { fExposedOn = day;         }    
    void SetTracedOn(int day)        { fTracedOn.push_back(day); }
    void SetDayLastTestedOn(int day) { fDayLastTestedOn = day;   }
    
    void AddUninfectedContacts(int day, int number) { 
      fUninfectedContacts[day] = number; 
    }
    void SetInfectionStatus(int day,CVInfectionStatus status) {
       //! do not re-infect if previously recovered
      if (day >=  fNDays) return;
      if(fRecoveredOn  > 0 && status != S_Recovered) return;
      if(fExposedOn    < 0 && status == S_Exposed  ) fExposedOn   = day;
      if(fRecoveredOn  < 0 && status == S_Recovered) fRecoveredOn = day;
      if(fInfectiousOn < 0 &&(status == S_ExposedInfectious || status == S_Infectious)) fInfectiousOn = day;
      fInfectionStatus.at(day) = status;        
    } 
    void SetQuarantineStatus(int day, bool status) { 
      if (day <  fNDays) fQuarantineStatus.at(day) = status; 
    }
    void SetInfectedBy(int id) { fInfectedBy = id; }
        
    //! give the person the app
    void SetHasApp(bool hasApp=true) { fHasApp = hasApp; }
    void SetDoesReport(bool doesReport=true) { fDoesReport = doesReport; }
    void ClearExposed()     { fExposed.clear();       }
    void AddExposed(int id) { fExposed.push_back(id); }
    void ClearDays()        { fDays.clear();          }
    void AddDay(int day)    { fDays.push_back(day);   }
    void AddExposedPerson(int id, int day) { 
      AddExposed(id);
      AddDay(day); 
    }
    
    //! methods
    bool Trace(int startday, int endday) {
      SetTracedOn(startday);
      bool traced = false;
      for (int iday = startday; iday <=endday && iday < fNDays; iday++) {
        if(fTracingStatus.at(iday) == S_Reported) break;
        fTracingStatus.at(iday) = S_Traced;
        traced = true;               
      }
      return traced;
    }
    void Expose(int day,CVDisease* disease) {
      fExposedOn = day;
      fDisease = disease;
      //! get course of disease
      vector<CVInfectionStatus> course = fDisease->DrawCourse();
      int availableDays = fInfectionStatus.size() - day;
      if (availableDays < 0) cout << "ERROR: fInfectionStatus out of range" << endl;
      int sicknessDuartion = course.size();
      int effectiveDays = TMath::Min(availableDays, sicknessDuartion);    
      copy(course.begin(),course.begin() + effectiveDays,fInfectionStatus.begin()+day); //FIXME    
      fill(fInfectionStatus.begin()+day+effectiveDays,fInfectionStatus.end(),S_Recovered);
      fRecoveredOn = day+course.size();
      fSymptomOnset = fRecoveredOn - fDisease->GetDaysForRecovery();
      //! evaluate if person gets sympomatic
      fHasSymptoms = fDisease->DrawHasSymptoms();     
    }
    void Expose(int day,CVPerson* infector) {
      fInfectedBy = infector->GetId();
      infector->AddExposedPerson(fId,day);
      Expose(day,infector->GetDisease());      
    }        
    int Quarantine(int day,int quarantine=40) {
      //! recovered persons that knew they were sick, do not go in quarantine
      if(fRecoveredOn > 0 && fRecoveredOn < day && fReportedOn > 0) return day;
      fQuarantinedOn.push_back(day);
      int endday = min(day+quarantine, fNDays);
      //! persons that recover in quarantine and knew they were sick, leave quarantine at recovery
      if(fRecoveredOn > 0 && fReportedOn > 0) endday = min(endday,fRecoveredOn); 
      int counter = 0;    
      for(int jday=day;jday<endday && jday < fNDays;jday++) { 
        if (GetQuarantineStatus(jday) == false) { 
          counter++; 
          SetQuarantineStatus(jday, true); 
        } 
      }
      return day + counter;
    }
    int ReleaseFromQuarantine(int day,int quarantine=40) {
      int daysremoved = 0;
      for(int iday=day; iday<day+quarantine && iday<fNDays; iday++) {
        if(GetQuarantineStatus(iday) == true) { 
          SetQuarantineStatus(iday, false);
          daysremoved++; 
        }
        else break;
      }
      return daysremoved;
    }    
    //! set reported and quarantine
    int Report(int day,int quarantine=40) {
      //bool notrecovered = true;
      //int iday = day;
      //while (notrecovered && iday < fInfectionStatus.size()) {
      //  if (fInfectionStatus.at(iday)!=S_Recovered) {
      //    if (GetReportedOn() < 0) fReportedOn = iday;
      //    try { fTracingStatus.at(iday) = S_Reported; }
      //    catch(...){ ; }
      //  }
      //  else notrecovered = false;
      //  iday++;
      //}
      for(int iday=day; iday < fNDays; iday++) {
        if(fInfectionStatus.at(iday)!=S_Recovered) {
          if(GetReportedOn()<0) fReportedOn = iday;
          fTracingStatus.at(iday) = S_Reported;
        }
        else break;
      }
      return Quarantine(day,quarantine);
    }

  protected:
    int        fId;      //! each person has a numeric ID
    CVDisease* fDisease; //! ... and might carry a disease
    int        fNDays;   //! Number of days in simulation
    
    int  fExposedOn;            //! day of exposure    
    int  fInfectiousOn;         //! the day this person becomes infectious    
    int  fSymptomOnset;         //! the day this person becomes symptomatic
    int  fReportedOn;           //! day when we were reported
    int  fRecoveredOn;          //! day on which we recovered
    int  fDayLastTestedOn;      //! day of most recent test   
    vector<int> fQuarantinedOn; //! days when person was put in quarantine 
    vector<int> fTracedOn;      //! days when we were traced
    
    bool fHasSymptoms;          //! does this person show symptoms if infected?
    bool fHasApp;               //! does this person use the tracing app?
    bool fDoesReport;           //! does this person go to the doctor if they have symptoms?
    vector<CVInfectionStatus> fInfectionStatus;  //! infection status for each day in the pandemic
    vector<CVTracingStatus>   fTracingStatus;    //! tracing or reporting status for each day in the pandemic
    vector<bool>              fQuarantineStatus; //! is this person in quarantine?
    
    int          fInfectedBy; //! ID of person who infected this one
    vector<int>  fExposed;    //! all the IDs of people this person infected
    vector<int>  fDays;       //! ... and the days it happened
    map<int,int> fUninfectedContacts; //! map of day and how many uninfected contacts were had
    
  private:
    static int nextId;
    TRandom3* fRandom;
    bool      fInternalRandom; 
};
//int CVPerson::nextId = 0;

#endif
