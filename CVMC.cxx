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
 
#include "CVMC.h"

using namespace std;

int CVPerson::nextId = 0; //don't like it here, but otherwise we get compiler issues

void CVMC::DoMC() { 
  if(fDebug) cout << "* starting run " << fRunId << endl;
  
  //! initialize stopwatch to take runtime
  TStopwatch stopwatch;
  Reset();
  //! people stacks
  list<CVPerson*> inRotation;   //! people who are exposed and not yet recovered
  list<CVPerson*> inQuarantine; //! people who are healthy but in quarantine
  
  list<CVPerson*> toErase; //! people who have recovered and should be removed before the next loop
  
  vector<int> timeOrderedListOfInfectedIDs; //! People are added here sorted by day they were infected; this is used in the dot output chart
  //! text ouput
  if (fDebug) {
    fOutputTSV.open(Form("%s_%d.txt",fOutputPrefix.c_str(),fRunId)); 
    AddHeaderToTSV();
  }
  //! root output
  fOutputRoot = new TFile(Form("%s_%d.root",fOutputPrefix.c_str(),fRunId),"recreate");
  
  //! seed "patient 0"
  fPersons.front()->Expose(0,fDisease);  
  fPersons.front()->SetInfectedBy(-3);
  inRotation.push_back(fPersons.front()); 
  timeOrderedListOfInfectedIDs.push_back(0);
  int totalReported = 0; 
  int totalSick = 0;
  //! Loop over days in the outbreak
  for(int iday = 0; iday < fNDays-1; iday++) {  
    //! Infrastructure for .dot output.
    if (fPeopleInDotFile < fMaxPeopleInDotFile) fLastDayForDot = iday;  
    if (fPeopleInDotFile < 2.*fMaxPeopleInDotFile) fLastDayForDotSimple = iday;
    if (iday%50 == 0) cout << "Passed day " << iday << endl;
    if(fDebug) cout << endl << "**** Working on day " << iday << endl;
    fNExposed = 0; fNInfectious = 0; fNTraced = 0; fNReported = 0; fNQuarantineToday=0;  
      
    //! remove people who have recovered
    for(auto kv : toErase) inRotation.remove(kv);
    toErase.clear();
    if(inRotation.size() < 1) { 
      if (fDebug) cout << "No sick people left. Exiting day loop early." << endl; 
      //! if no one is infected we don't have to keep going
      break; 
    }   
      
    int ninfector = 0; //! < Needed for R_e; how many people were sick today?
    int ninfectee = 0; //! < Needed for R_e; how many people were newly exposed today
    if (fDebug) cout << "Total reported " << totalReported << endl;
    //! Link the start of interventions to the number of reported cases
    if (float(totalSick)/float(fNPersons) > fStartTracingTestingInfectedFraction && fStartTracingOnDay > 9998) { 
      fStartTracingOnDay = iday; fStartTestingOnDay = iday; fSocialDistancingFrom = iday;
      if (fDebug) cout << "Started tracing, testing, and social distancing on day " << iday << " with " << totalReported << "/" << fNPersons  << "=" << float(totalReported)/float(fNPersons)<< "  > " << fStartTracingTestingInfectedFraction << endl;
    }
    
    //! iterate on all people who are exposed and not yet recovered
    for(auto kv : inRotation) {
      if(fDebug) cout << endl << "=> person " << kv->GetId() << " exposed on day " << kv->GetExposedOn() << endl;         

      //! get current status
      CVInfectionStatus currentInfectionStatus  = kv->GetInfectionStatus(iday);
      CVTracingStatus   currentTracingStatus    = kv->GetTracingStatus(iday);
      bool              currentQuarantineStatus = kv->GetQuarantineStatus(iday);
      if(fDebug) cout << " Day " << iday <<  " status:  infection " << currentInfectionStatus << ", tracing " << currentTracingStatus << ", quarantine " << currentQuarantineStatus << endl;
      IncrementPopulationStatistics(currentInfectionStatus, currentTracingStatus);
      //! Do some accounting if person recovered today
      if(currentInfectionStatus == S_Recovered) { 
        if(fDebug) {  cout << " 1 person removed from rotation " << endl;  fPeopleInDotFile++;  }
        //! use them to calcualte R_e (we can't do it sooner as only here they have infected everyone they could have)
        ninfector++;
        ninfectee = ninfectee + kv->GetNExposed();
        FillDiagnostics(kv, iday >= fStartTestingOnDay);
        //! add to remove pile
        toErase.push_back(kv); 
        continue;  
      }
      
      //! find new victims if infectious
      if((currentInfectionStatus == S_ExposedInfectious || currentInfectionStatus == S_Infectious) 
        //! .. but not quarantined
        && currentQuarantineStatus == false) {
        //! get todays infectiousness
        float infProb = kv->GetDisease()->GetInfectiousness(iday-kv->GetSymptomOnset(),kv->GetHasSymptoms());
        //! .. and the number of constacts
        int meettoday = GetPeopleMetToday(iday);
        if(fDebug) cout << " will meet " << meettoday << " and infect with probability " << infProb <<  endl;
        //! draw number of victims
        int ninfected = fRandom.Binomial(meettoday,infProb); 
        if(fDebug) cout << " will infect " << ninfected << " others." <<  endl;
        //! .. and remember how many people we met and did not infect (so we know how many to trace)
        kv->AddUninfectedContacts(iday, meettoday-ninfected);
        //! Now infect the victims ...
        for(int iinfect = ninfected; iinfect--;) {
          //! randomly pick the victim's ID
          int pp = fRandom.Integer(fNPersons);
          //! .. but not ourselves
          if (pp == kv->GetId()) pp = (pp+1)%fNPersons;
          //! only infect the person if not infected before and not in quarantine
          if(fPersons.at(pp)->GetExposedOn()<0 && !fPersons.at(pp)->GetQuarantineStatus(iday)) {
            if (fDebug) { cout << " adding exposed person " << pp << endl; AddConnectionToDot(kv->GetId(), pp, 0, iday); }
            timeOrderedListOfInfectedIDs.push_back(pp); 
            fPersons.at(pp)->Expose(iday,kv);
            totalSick++;
            inRotation.push_back(fPersons.at(pp));
      } } }
      
      //! see if we can report a person .. only if not already reported
      if(kv->GetReportedOn()<0) {
        bool doesReport = false;
        int reportday = iday;
        //! .. if showing symptoms and person willing to see doctor/has access to test
        if(currentInfectionStatus == S_Infectious && kv->GetHasSymptoms() && kv->GetDoesReport()) { 
          doesReport = true; 
        }
        //! .. or if traced and then tested positive
        if(currentTracingStatus == S_Traced  && iday > fStartTestingOnDay) { 
          //! get days since the most recent time we got traced
          int daysSinceTraced = iday - kv->GetTracedOn(kv->GetNTracedOn()-1);
          //! .. and tested
          int daysSinceTested = kv->GetDayLastTestedOn() > -1 ? iday - kv->GetDayLastTestedOn() : 9999;          
          if (fDebug) cout << "Days since traced: " << daysSinceTraced << " and days since tested: " << daysSinceTested << endl;
          //! if tracing is recent (today or yesterday) and no test was performed ...
          if(daysSinceTraced <=1 && daysSinceTested > 2) {
            //! .. test them
            if(kv->GetTestsPositive(iday) ) { 
              doesReport = true; 
              //! .. but wait for the test result
              reportday=reportday+fDaysToTestResult;
              if (fDebug) cout << "Person " << kv->GetId() << " will be reported from testing positive on day " << reportday << endl;
          } }
          //! .. or re-test after some time in quarantine, this should cover the latent period
          else if (daysSinceTested == fdTTest) { 
            if (kv->GetTestsPositive(iday) ) { 
              doesReport = true; 
              reportday=reportday+fDaysToTestResult;
              if (fDebug) cout << "Person " << kv->GetId() << " will be reported from testing positive on re-test on day " << reportday << endl;
            }
            //! .. if this turns out negative, release the person
            else {
              kv->ReleaseFromQuarantine(iday+fDaysToTestResult,fDaysInQuarantine);
              if (fDebug) cout << "Person " << kv->GetId() << " released from quarantine due to negative second test." << endl;
        } } }
        //! .. or if selected for a random test
        else if(fRandomTesting && iday > fStartTestingOnDay && fRandom.Uniform() < fRandomTestingRate) {
          //! .. that is performed once
          if(kv->GetTestsPositive(iday)) {
            doesReport = true; 
            reportday=reportday+fDaysToTestResult;  
            if (fDebug) cout << "Person " << kv->GetId() << " will be reported from testing positive by random test on day " << reportday << endl;          
        } }        
        if(doesReport) {
          kv->Report(reportday,fDaysInQuarantine); 
          totalReported++;
      } } 

      //! reporting triggers tracing of
      if (kv->GetReportedOn()+fTracingDelay == iday && iday > fStartTracingOnDay && kv->GetHasApp()) {
        if(fDebug) cout << " tracing " << kv->GetNExposed()  << " people " << endl;
        //! .. contacts you infected
        for(auto kvv : kv->GetExposed()) {
            Trace(fPersons.at(kvv),iday,kv->GetId(),1);
        } 
        //! .. the contact you got infected from (unless you are patient 0)
        if(fBackwardTracing && kv->GetInfectedBy() >= 0) { 
          Trace(fPersons.at(kv->GetInfectedBy()),iday,kv->GetId(),2);  
        }
        //! .. and a bunch of people you may have met but not infected
        TraceUninfected(kv,iday);        
      } 
               
    } //! end of sick people on this day
    // check who is in quarantine
    for (auto qkv : fPersons) { if (qkv->GetQuarantineStatus(iday)) fNQuarantineToday++; }
    //! figure out R_t. 
    fNumberRecoveredByDay.push_back(ninfector);
    fNumberInfectiousByDay.push_back(fNInfectious);
    numberNewlyInfectedByDay.push_back(ninfectee);
    AddDayToTSVAndROOT(iday);
    fLastDayWithPatients = iday;
  } //! end of days
   fNExposedTotal = 0;
  for(auto kp : fPersons) {
    if (kp->GetExposedOn() > 0) { fNExposedTotal ++; }    
  }
  fPeopleInDotFile = 0;
  if (fDebug) {
    for (auto pp :  timeOrderedListOfInfectedIDs) {
      if (fPersons.at(pp)->GetExposedOn() > -1) AddPersonToDot(fPersons.at(pp));
    }
  }
  stopwatch.Print();
  if (fDebug) {
    fOutputDot.open(Form("%s_%d.dot",fOutputPrefix.c_str(),fRunId));
    fOutputDot <<"digraph infection { \n \
                        newrank=true; \n \
                        compound=true; \n \
                        node [fontname=\"Helvetica\" fontsize=18]; \n" << endl;
    fOutputDotSimple.open(Form("%s_Simple_%d.dot",fOutputPrefix.c_str(),fRunId));
    fOutputDotSimple <<"digraph infectionsimple {  \n \
                        newrank=true; \n \
                        node [fontname=\"Helvetica\" fontsize=18]; \n" << endl;
  
    AddDaysToDot();
    //AddLegendToDot();
    fOutputDot << fDotStringPeople << endl;
  
    fOutputDot <<"}" << endl;  
    fOutputDot.close();

    fOutputDotSimple << fDotStringPeopleSimple << endl;  
    fOutputDotSimple <<"}" << endl;  
    fOutputDotSimple.close();
    
    fOutputTSV.close();
    WriteGnuplotScript();    
  }

  WriteROOTFile();
  cout << "Exposed total: " << fNExposedTotal << "( "<< float(fNExposedTotal)/float(fNPersons) * 100. <<"% of population)" << endl;
}

void CVMC::TraceUninfected(CVPerson *aperson, int day) {
  if (!fTraceUninfected) { return; }
  if(aperson->GetHasApp() == 0) { return; }
  if(fDebug) cout << " person " << aperson->GetId() << " met " << pplMet << " others in " << fDaysBackwardTrace << " days of tracing lookback time " << endl;
  for (int iday = max(0, day - fDaysBackwardTrace); iday <= day; iday++) { // loop over days in the backward trace
    int pplMet = aperson->GetUninfectedContactsDay(iday);
    if (pplMet < 0) { //! We have not yet drawn how many people we meet today, so do it now
       pplMet = GetPeopleMetToday(iday);
    }
    for (int ippl = 0; ippl < pplMet; ippl++) {
      if (fRandom.Uniform() > fTracingEfficiency) continue; //! tracing was not successful with this person
      int pp = fRandom.Integer(fNPersons); //! randomly pick the ID of the contact       
      //! .. but not ourselves
      if (pp == aperson->GetId()) pp = (pp+1)%fNPersons;
      //! we did not really meet this person because they were in Q on that day, or the trace is unsuccessful because the contact doesn't have the app
      if (fPersons.at(pp)->GetQuarantineStatus(iday) || fPersons.at(pp)->GetHasApp() == 0) continue;
      //! WARNING: instead of the usuall quarantine time uninfected contacts enter 
      //! for the time until the second test, this emulates 100% test efficiency and no false positives
      //! and disregards the possibility that the person was infected by someone else in the meantime.
      fPersons.at(pp)->Quarantine(day,fdTTest); // Note: People sent to quarantine will ignore the order if they have previously been reported (ie they know they've been infected and recovered already)
    }
  }
}
            
//! Helper function for Trace. Return value is whether or not we did the trace
bool CVMC::DoTrace(CVPerson *aperson, int day, int fromPersonID, int tlevel, int direction) {
  if (aperson->GetHasApp() == 0 || fPersons.at(fromPersonID)->GetHasApp() == 0) {
    AddConnectionToDot(fromPersonID, aperson->GetId(),98,day);
    return false;
  }
  //! figure out if person has been reported before today
  bool isPersonAlreadyKnown = (aperson->GetReportedOn() != -1)   && (aperson->GetReportedOn() <= day);
                            
  //! Person has already been reported; do not trace again. 
  //! This misses the edge case where the trace came due to a second order trace, 
  //! and we now trace with a first order trace (that would trigger a second order)
  if (isPersonAlreadyKnown) return false; 

  //! figure out if we're in time range for tracing
  int contactDifference;
  //! fromPersonID infected aperson; forward trace
  if (direction == 1) contactDifference = day - aperson->GetExposedOn();
  else contactDifference = day - fPersons.at(fromPersonID)->GetExposedOn();    

  //! contact too far in the past to be traced
  if(contactDifference > fDaysBackwardTrace) { 
    AddConnectionToDot(fromPersonID, aperson->GetId(),99,day); 
    return false; 
  }
  //! figure out if we miss the contact by chance
  if(fRandom.Uniform() > fTracingEfficiency) {
    AddConnectionToDot(fromPersonID, aperson->GetId(),98,day);
    return false;
  }
  //! Set the tracing status
  bool traced = aperson->Trace(day, day+fDaysInQuarantine); // they keep 'traced' status as long as one might be in quarantine
  if (traced) { 
    AddConnectionToDot(fromPersonID, aperson->GetId(),tlevel,day);
    if (fDebug) cout << "Dotrace: traced person " << aperson->GetId() << "from person " << fromPersonID << endl;
  }
  return traced;
}

//! 'aperson' is someone who was exposed by a person whose contacts we start tracing on day 'day'
//! we trace only if 'aperson' was exposed at most setting.daysBackwardTrace days before 'day'
//! Direction=1 : forward trace; direction=2: backward trace
void CVMC::Trace(CVPerson *aperson, int day, int fromPersonID, int direction) {
  if (direction == 2 && !fBackwardTracing) return;
  if(DoTrace(aperson, day, fromPersonID, 1, direction)) { 
   aperson->Quarantine(day,fDaysInQuarantine); //! if we traced, we also quarantine
  }
  else return; //! we did not trace this person, so we don't trace 2nd order from them
  
  if (fTracingOrder < 2) return; //! we also don't trace 2nd order if that is turned off in settings
  else if (fTracingOrder > 2 || fTracingOrder < 1) { 
    cout << "Trace:Parameter error: tracingOrder must be 1 or 2; is currently " << fTracingOrder << endl; 
    return; 
  }
  //! We confirmed that we trace 2nd order contacts now
  //! If we're on a backward trace, continue to back-trace infector
  //! If we are on a forward trace, this is not necessary as the infector is the original index case
  if (direction==2 && fBackwardTracing && aperson->GetInfectedBy() >=0)  {
    if (DoTrace( fPersons.at(aperson->GetInfectedBy()) , day, aperson->GetId(), 2, 2) ) {
      fPersons.at(aperson->GetInfectedBy())->Quarantine(day,fDaysInQuarantine); 
    }
  }
  //! For the second order trace, aperson is now the index case. Therefore, trace their uninfected contacts ...
  TraceUninfected(aperson,day);
  //! And their remaining contacts
  for (int iexp = 0; iexp < aperson->GetNExposed(); iexp++) {
    if (aperson->GetExposedPersonID(iexp) == fromPersonID) { continue; } //! if this is coming from a backward trace, we don't want to go back -> forth -> back
    CVPerson *pExp2ndOrder = fPersons.at(aperson->GetExposedPersonID(iexp));
    if (DoTrace(pExp2ndOrder, day, aperson->GetId(), 2, 1)) {
      pExp2ndOrder->Quarantine(day,fDaysInQuarantine);
    }
  }  
}//trace end


//!Diagnostic histograms
void CVMC::FillDiagnostics(CVPerson *kv, bool posttracing) {
  fhIncubationPeriod->Fill(kv->GetSymptomOnset() - kv->GetExposedOn());        
  for(int iexp = 0; iexp < kv->GetNExposedDays(); iexp++) {
    if (kv->GetHasSymptoms()) { 
      fhInfectiousnessAgeSymptom->Fill(kv->GetExposedDays(iexp) - kv->GetExposedOn());            
    }
    else {
      fhInfectiousnessAgeNoSymptom->Fill(kv->GetExposedDays(iexp) - kv->GetExposedOn());
      fhInfectiousness->Fill(kv->GetExposedDays(iexp) - kv->GetSymptomOnset());
    }
  }
  if (kv->GetNExposed() > 0) {
    fhLatentPeriod->Fill(kv->GetExposedDays(0) - kv->GetExposedOn()); 
  }
  if (kv->GetNExposed() == 1) {
    if (kv->GetHasSymptoms()) fhInfectiousnessDurationSymptom->Fill(1);
    else fhInfectiousnessDurationNoSymptom->Fill(1);      
  }          
  else if (kv->GetNExposed() > 1) {
    int fillwith  = kv->GetExposedDays(kv->GetNExposedDays()-1) - kv->GetExposedDays(0) + 1;
    if (fillwith == 0) { cout << "fill 0 " << kv->GetExposedDays(kv->GetNExposedDays()-1) << ", " << kv->GetExposedDays(0) << endl;}
    if (kv->GetHasSymptoms()) fhInfectiousnessDurationSymptom->Fill(fillwith); // need to add one because we want the inclusive number of days of infectivity ...
    else fhInfectiousnessDurationNoSymptom->Fill(fillwith); // ... if the first person is infected on day n and the second on day n+1, the index case is infectious for 2 days  (n+1) - n + 1 = 2
   }
 
   int daytoquarantine = 30;
   if (kv->GetNQuarantinedOn() > 0) { daytoquarantine = kv->GetQuarantinedOn(0) - TMath::Max(kv->GetSymptomOnset()-2, kv->GetExposedOn());  }
   int daytoreported = -2;
   if (kv->GetReportedOn() > -1) daytoreported = kv->GetReportedOn() - kv->GetExposedOn();
   int daytotested = -2;
   if (kv->GetDayLastTestedOn() > -1) daytotested =kv->GetDayLastTestedOn() - kv->GetExposedOn();
   if (posttracing) {
     fhDaysToQuarantinePostIntervention->Fill(daytoquarantine);
     fhDaysToReportedPostIntervention->Fill(daytoreported);
     fhDaysToTestedPostIntervention->Fill(daytotested);
   }
   else {
     fhNumberInfectedPreIntervention->Fill(kv->GetNExposed());
     fhDaysToQuarantinePreIntervention->Fill(daytoquarantine);
     fhDaysToReportedPreIntervention->Fill(daytoreported);
     fhDaysToTestedPreIntervention->Fill(daytotested);
   }

}
//******* DOT output *******/
void CVMC::AddPersonToDot(CVPerson *kv) {
  if (fPeopleInDotFile > 2.*fMaxPeopleInDotFile) { return; }
  string colorOfAsymptomatic = "fillcolor=\"#ECF7F9\""; // light blue
  string colorOfSymptomatic = "fillcolor=\"#D0C4EE\""; // mid purple
  int penwidth = 6;
  string colorOfTraced = "color=\"#61304e\" style=\"dashed, filled\" "; // dark purple
  string colorOfReported = "color=\"#61304e\" style=\"filled\" "; // dark purple
  string colorOfNoStatus = "color=\"#C8DCDE\" style=\"filled\" "; // dark purple
  string DotOptionsPerson = "";
  if (kv->GetNTracedOn() > 0) DotOptionsPerson = colorOfTraced;
  if (kv->GetReportedOn() > 0) DotOptionsPerson = colorOfReported; // reported overwrites traced
  if ( !(kv->GetNTracedOn() > 0) && !(kv->GetReportedOn() > 0) ) DotOptionsPerson = colorOfNoStatus;
  if (kv->GetHasSymptoms() > 0 && kv->GetExposedOn() > 0) DotOptionsPerson = DotOptionsPerson + colorOfSymptomatic;
  else { DotOptionsPerson = DotOptionsPerson + colorOfAsymptomatic; }


  if (fPeopleInDotFile < fMaxPeopleInDotFile) {
    string tracedstring = "T ";
    string reportedstring = "R ";
    string qdstring = "Q ";
    string appstring = "";    
    if (kv->GetNTracedOn() == 0) tracedstring = "T--";
    if (kv->GetReportedOn() == -1) reportedstring = "R--";
    else reportedstring = Form("R %d", kv->GetReportedOn());
    if (kv->GetNQuarantinedOn() == 0) qdstring = "Q--";
    for (int itrace = 0; itrace < kv->GetNTracedOn(); itrace++) { 
      if (itrace == 0) {tracedstring = Form("%s %d",tracedstring.c_str(), kv->GetTracedOn(itrace)); }
      else tracedstring = Form("%s, %d",tracedstring.c_str(), kv->GetTracedOn(itrace)); 
    }
    for (int iq = 0; iq < kv->GetNQuarantinedOn(); iq++) { 
      if (iq == 0) {qdstring = Form("%s %d",qdstring.c_str(), kv->GetQuarantinedOn(iq)); }
      else qdstring = Form("%s, %d",qdstring.c_str(), kv->GetQuarantinedOn(iq)); 
    }    
    if (kv->GetHasApp() == 1) appstring = "\\nApp";
    fDotStringPeople = fDotStringPeople+Form("\t subgraph cluster_%d { \n \
                                              \t %s  penwidth=%d    \n \
                                              \t %d [shape=record label=\"P%d%s|{E %d|I %d|R %d}|{%s|%s|%s}\" penwidth=0] \n } \n",
                                                            kv->GetId(),
                                                            DotOptionsPerson.c_str(),
                                                            penwidth,
                                                                        kv->GetId(),                                                                        
                                                                        kv->GetId(),
                                                                        appstring.c_str(),
                                                                        kv->GetExposedOn(),
                                                                        TMath::Max(kv->GetSymptomOnset()-2,kv->GetExposedOn()), // infectious starts 2 days before symptoms, unless incubation period is less than 2 days
                                                                        kv->GetRecoveredOn(),
                                                                        tracedstring.c_str(), 
                                                                        reportedstring.c_str(), 
                                                                        qdstring.c_str()            ) ;  
    if ( kv->GetInfectedBy() > -1) { 
      if (kv->GetExposedOn() == fPersons.at(kv->GetInfectedBy())->GetExposedOn() ) 
      fDotStringPeople = fDotStringPeople+Form("\t { rank=same %d Dh%d } \n",kv->GetId(), kv->GetExposedOn()); 
      else fDotStringPeople = fDotStringPeople+Form("\t { rank=same %d D%d } \n",kv->GetId(), kv->GetExposedOn());
    }
    else fDotStringPeople = fDotStringPeople+Form("\t { rank=same %d D%d } \n",kv->GetId(), kv->GetExposedOn());
  }
  else if (fPeopleInDotFile < 2.*fMaxPeopleInDotFile) {
    fDotStringPeople = fDotStringPeople+Form("\t subgraph cluster_%d { \n \
                                              \t %s  penwidth=%d    \n \
                                              \t %d [shape=circle label=\"\"] \n \t } \n", 
                                              kv->GetId(), 
                                              DotOptionsPerson.c_str(),
                                              penwidth,
                                              kv->GetId());
     fDotStringPeople = fDotStringPeople+Form("\t { rank=same %d D%d } \n", kv->GetId(), kv->GetExposedOn()); 
  }
  
  if (fPeopleInDotFile < 2.*fMaxPeopleInDotFile) {
    penwidth = 8;
   fDotStringPeopleSimple = fDotStringPeopleSimple+Form("\t %d [shape=circle %s penwidth=%d label=\"\"] \n", kv->GetId(), DotOptionsPerson.c_str(),penwidth) ;  
   fDotStringPeopleSimple = fDotStringPeopleSimple+Form("\t { rank=same %d D%d } \n", kv->GetId(), kv->GetExposedOn());
  }
  
  fPeopleInDotFile++;    
}

void CVMC::AddLegendToDot() {
 string colorOfInfect = "color=\"#616161\" penwidth = 2"; //
  string colorOfTrace1 =  "color=\"#B83C45\" style = dashed penwidth = 4"; //
  string colorOfTrace2 =  "color=\"#E38D2A\" style = dashed penwidth = 4" ; // 
  string colorOfMissed =  "color=\"#E8E06E\" style = dotted penwidth = 2" ; //  
  string colorOfTooFarBack =  "color=\"#FCF485\" style = dotted penwidth = 2" ; //  
  
  string colorOfAsymptomatic = "fillcolor=\"#E2F5FA\""; // light blue
  string colorOfSymptomatic = "fillcolor=\"#BCB2E2\""; // dark blue
  string colorOfTraced = "penwidth = 5 color=\"#61304e\" "; // dark purple
  string colorOfReported = "penwidth = 5 color=\"#61304e\" "; // dark purple
  string colorOfNoStatus = "penwidth = 5 color=\"#C2DCDE\""; // dark purple
 
    
  fOutputDot << "\t subgraph cluster_legend {" << endl;
  fOutputDot << "\t subgraph cluster_clean \n \t { penwidth=1 \n pclean [shape=record label=\"Person ID|{Exposed|Infectious|Recovered}|{Traced|Reported|Quarantined}\" penwidth=0] }" << endl;
  fOutputDot << Form("\t psymptom [shape=box label=\"Person symptomatic\" style=filled %s]",colorOfSymptomatic.c_str()) << endl;  
  fOutputDot << Form("\t pnosymptom [shape=box label=\"Person asymptomatic\" style=filled %s]",colorOfAsymptomatic.c_str()) << endl;  
  fOutputDot << Form("\t prep [shape=box label=\"Person reported\" %s ]", colorOfReported.c_str()) << endl;
  fOutputDot << Form("\t ptrac [shape=box label=\"Person traced\" %s]",colorOfTraced.c_str()) << endl;  
  fOutputDot << Form("\t pnostatus [shape=box label=\"Person unreported and untraced\" %s]",colorOfNoStatus.c_str()) << endl;  
  fOutputDot << "\t pclean -> pnostatus [style=\"invis\"]" << endl;  
  fOutputDot << "\t pnostatus -> prep [style=\"invis\"]" << endl;  
  fOutputDot << "\t prep -> psymptom [style=\"invis\"]" << endl;  
  fOutputDot << "\t prep -> ptrac [style=\"invis\"]" << endl;  
  fOutputDot << "\t psymptom -> dummy1 [style=\"invis\"]" << endl;    
    fOutputDot << "\t {rank=same prep ptrac }" << endl;  
    fOutputDot << "\t {rank=same psymptom pnosymptom }" << endl;  
  for (int idummy = 1; idummy < 11; idummy++) {
    fOutputDot << Form("\t dummy%d [shape=box label=\"\" style=\"invis\"] ",idummy) << endl;
  }
  for (int idummy = 1; idummy < 10; idummy = idummy+2) { 
   fOutputDot << Form("\t { rank=same dummy%d dummy%d }",idummy, idummy+1) << endl;
  }
  for (int idummy = 1; idummy < 9; idummy = idummy+2) { 
   fOutputDot << Form("\t dummy%d -> dummy%d [style=\"invis\"]",idummy, idummy+2) << endl;
  }  
  fOutputDot << "\t pnostatus -> dummy1 [style=\"invis\"]" << endl;  

  fOutputDot << Form("\t dummy1 -> dummy2 [label=\"infection\" %s]",colorOfInfect.c_str()) << endl;  
  fOutputDot << Form("\t dummy3 -> dummy4 [label=\"day [first order trace]\" %s]",colorOfTrace1.c_str()) << endl;  
  fOutputDot << Form("\t dummy5 -> dummy6 [label=\"day [second order trace]\" %s]",colorOfTrace2.c_str()) << endl;  
  fOutputDot << Form("\t dummy7 -> dummy8 [label=\"day [missed trace]\" %s]",colorOfMissed.c_str()) << endl;  
  fOutputDot << Form("\t dummy9 -> dummy10 [label=\"day [trace too far back]\" %s]",colorOfTooFarBack.c_str()) << endl;    

  fOutputDot << "\t { rank=same pclean D5}" << endl; // put the legend somewhere near the top
  fOutputDot << "\t }" << endl;

}
/** IDfrom and IDto are the ids of the people we trace from and to
 *  tlevel 0: we infect
 * tevel 1:  level 1 trace
 * tlevel 2: level 2 trace
**/
void CVMC::AddConnectionToDot(int IDfrom, int IDto, int tlevel, int tday) {
  if (!fDebug) return;
  if (fPeopleInDotFile > 2.*fMaxPeopleInDotFile) return;
  string colorOfInfect = "color=\"#616161\" "; //
  int penwidth = 6;
  string colorOfTrace1 =  "color=\"#CE2C1F\" style = dashed fontcolor=\"#A13647\" "; // dark red
  string colorOfTrace2 =  "color=\"#E38D2A\" style = dashed fontcolor=\"#D9802F\" "; // orange
  string colorOfMissed =  "color=\"#CBC56F\" style = dotted fontcolor=\"#CBC56F\" " ; //  dirty yellow
  string colorOfTooFarBack =  "color=\"#BBD575\" style = dotted fontcolor=\"#BBD575\" " ; //  
  string DotOptions = "";
  if (tlevel == 0) { DotOptions = colorOfInfect; penwidth = 1; }
  else if (tlevel == 1) { DotOptions = colorOfTrace1; penwidth = 5; }
  else if (tlevel == 2) {DotOptions = colorOfTrace2;  penwidth = 5; }
  else if (tlevel == 98) {DotOptions = colorOfMissed;penwidth = 4; }
  else if (tlevel == 99) {DotOptions = colorOfTooFarBack;penwidth = 4; }
  else cout << "CVMC::AddConnectionToDot: Warning, invalid trace level." << endl;
  
  if (fPeopleInDotFile < fMaxPeopleInDotFile)  {
    if (tlevel == 0) { fDotStringPeople = Form("%s \t %d -> %d [%s ltail=cluster_%d arrowsize=2  lhead=cluster_%d penwidth=%d] \n",fDotStringPeople.c_str(),IDfrom, IDto,DotOptions.c_str(),IDfrom, IDto,penwidth); } // do not add day if we just infect someone
    else { fDotStringPeople = Form("%s \t %d -> %d [%s label=\"D%d\" fontsize=18 fontname =\"Helvetica\" ltail=cluster_%d  lhead=cluster_%d arrowsize=2 penwidth=%d] \n",fDotStringPeople.c_str(),IDfrom, IDto,DotOptions.c_str(),tday,IDfrom, IDto, penwidth); }// only add tracing day
  }
  
  // for the simple plot
  colorOfTrace1 =  "color=\"#FF180F\" style = dashed fontcolor=\"#A13647\" "; // dark red
  colorOfTrace2 =  "color=\"#E38D2A\" style = dashed fontcolor=\"#D9802F\" "; // orange  
  if (tlevel == 0) { DotOptions = colorOfInfect; penwidth = 1; }
  else if (tlevel == 1) {DotOptions = colorOfTrace1; penwidth = 8; }
  else if (tlevel == 2) { DotOptions = colorOfTrace2; penwidth = 8; }
  if (tlevel < 10) {
    fDotStringPeopleSimple = Form("%s \t %d -> %d [%s penwidth=%d] \n",fDotStringPeopleSimple.c_str(),IDfrom, IDto,DotOptions.c_str(),penwidth); // don't add missed connections to simple output
    if (fPeopleInDotFile >= fMaxPeopleInDotFile) fDotStringPeople = Form("%s \t %d -> %d [%s penwidth=%d] \n",fDotStringPeople.c_str(),IDfrom, IDto,DotOptions.c_str(),penwidth); // don't add missed connections to reduced output
    }
  if (fDebug) { cout<<"\t " << IDfrom << " -> " << IDto << Form(" [%s]", DotOptions.c_str()) << endl;  }
}

void CVMC::AddDaysToDot() {
  fOutputDot <<"\t subgraph days {" << endl;
  fOutputDotSimple <<"\t subgraph dayssimple { \n \
                      \t   node [fontname=\"Helvetica\" fontsize=30 penwidth=0] \n" << endl;
  string option2 = "";
  for (int iday = 0; iday < fLastDayForDotSimple; iday++) {
      option2 = "style=\"invis\"";
      if (iday == fStartTracingOnDay) option2 = "style=filled, fillcolor=\"red\"";
      if (iday%5 == 0 && iday != fStartTracingOnDay) { option2 = "style=filled, fillcolor=\"gray\"";}
      if (iday == 0) option2 = option2 + "label=\"Days\" ";
        fOutputDot << "\t  D" << iday << " -> Dh" << iday << " [constraint=true style=\"invis\"] " << endl;
        fOutputDot << "\t  Dh" << iday << " -> D" << iday+1 << " [constraint=true style=\"invis\"] " << endl;        
        fOutputDot << "\t  D" << iday << Form(" [ label = \"Day\\n%d\" shape=square %s]", iday, option2.c_str()) << endl;
        fOutputDot << "\t  Dh" << iday << " [style=\"invis\"]"<< endl;                
      fOutputDotSimple << "\t  D" << iday << " -> D" << iday+1 << " [constraint=true] " << endl;
      fOutputDotSimple << "\t  D" << iday << Form(" [shape=square %s]", option2.c_str()) << endl;
  }
  fOutputDot <<"\t }" << endl; 
  fOutputDotSimple <<"\t }" << endl; 
}


//******* text output *******/
void CVMC::AddHeaderToTSV() {
  // Variables this class knows
  fOutputTSV<< "# fNDays = " << fNDays<< endl;
  fOutputTSV<< "# fNPersons = " << fNPersons << endl;
  fOutputTSV<< "# fStartTracingOnDay = " << fStartTracingOnDay << endl;  
  fOutputTSV<< "# fTracingOrder = " << fTracingOrder << endl;  
  fOutputTSV<< "# fStartTestingOnDay = " << fStartTestingOnDay << endl;    
  fOutputTSV<< "# fDaysToTestResult = " << fDaysToTestResult << endl;  
  fOutputTSV<< "# fDaysBackwardTrace = " << fDaysBackwardTrace << endl;  
  fOutputTSV<< "# fDaysInQuarantine = " << fDaysInQuarantine << endl;  
  fOutputTSV<< "# fAppProbability = " << fAppProbability << endl;  
  fOutputTSV<< "# fPeopleMetPerDay = " << fPeopleMetPerDay << endl;

/*        
  fOutputTSV<< "# SymptomProbability = " << setting.SymptomProbability << endl;
  fOutputTSV<< "# TestPositiveProbability = " << setting.TestPositiveProbability << endl;  
  fOutputTSV<< "# ExposureToInfectiousMeanDays = " << setting.ExposureToInfectiousMeanDays << endl;
  fOutputTSV<< "# ExposureToInfectiousMinDay = " << setting.ExposureToInfectiousMinDay << endl;
  fOutputTSV<< "# InfectiousToSymptomsMeanDays = " << setting.InfectiousToSymptomsMeanDays << endl;
  fOutputTSV<< "# InfectiousToSymptomsMinDay = " << setting.InfectiousToSymptomsMinDay << endl;
  fOutputTSV<< "# SymptomsToRecoveryMeanDays = " << setting.SymptomsToRecoveryMeanDays << endl;
  fOutputTSV<< "# SymptomsToRecoveryMinDay = " << setting.SymptomsToRecoveryMinDay << endl;
  fOutputTSV<< "# InfectionProbability = " << setting.InfectionProbability << endl;
  */
  fOutputTSV<< "# (1) day (2) nSusceptible  (3) nExposed  (4) nInfectious  (5) nRecovered  \
               (6) nTraced+nReported  (7) nQuarantine  (8) R_eff (5 day average)  (9)  \
                R_eff uncertainty  (10) doubling time  (11) doubling time uncertainty" << endl;
}

void CVMC::IncrementPopulationStatistics(CVInfectionStatus is, CVTracingStatus ts) {
  switch (is){
    case S_Exposed:           fNExposed++;    break;
    case S_ExposedInfectious:  fNInfectious++; break;
    case S_Infectious:        fNInfectious++; break;      
    case S_Recovered:         fNRecovered++;  break;  
    default: /*do nothing*/ break;
  }
  switch (ts){
    case S_Reported: fNReported++; break;
    case S_Traced:   fNTraced++;   break;
    default: /*do nothing*/ break;
  }
}


void CVMC::AddDayToTSVAndROOT(int day) {
  fNSusceptible = fNPersons - fNExposed - fNRecovered - fNInfectious;
  int ninfector = fNumberRecoveredByDay.at(day); int ninfectee = numberNewlyInfectedByDay.at(day);
   feffectiveR = -99.0; feffectiveRUncertainty = 0.;
  
  // calculate the doubling time
  float doublingTime = -99., doublingTimeUncertainty = -99.; //! Days needed to double number of new infections  
  float doublingDay = -99.;
  float doublingDayUncertainty = -99.0;
  // take the average of last 2 days to minimize fluctuations
  if (day > 4) {
    float lin_a, lin_m; // y = a + m*x   
    float dlin_a, dlin_m; // uncertainties. Use sqrt(y) for y as uncertainty
    int infectedTodayYesterday = float(ninfectee + fNumberRecoveredByDay.at(day-1))/2.0;
    int infected2to4daysAgo = float(fNumberRecoveredByDay.at(day-2) + fNumberRecoveredByDay.at(day-3))/2.0;
    lin_m = infectedTodayYesterday - infected2to4daysAgo;
    dlin_m = sqrt(infectedTodayYesterday+infected2to4daysAgo); // (x2-x1) = 1, hence m = (y2-y1)/(x2-x1) = (y2-y1)
    lin_a = infectedTodayYesterday - lin_m  * (day-0.5);
    dlin_a = sqrt(infectedTodayYesterday + TMath::Power(float(day)*dlin_m,2)    );
    if (abs(lin_m) > 0.005) {
    doublingDay = (2.*infectedTodayYesterday - lin_a )/lin_m;
    doublingDayUncertainty = sqrt(  TMath::Power(2./lin_m,2)*infectedTodayYesterday + TMath::Power(dlin_a/lin_m, 2) + TMath::Power( dlin_m*(2.*infectedTodayYesterday - lin_a)/(lin_m*lin_m) , 2)  );
    }
  }
  if ( doublingDay>0) doublingTime = doublingDay - day;
  if ( doublingDay>0) doublingTimeUncertainty = abs( doublingTime - (doublingDay+doublingDayUncertainty -day) );

  if (ninfector > 0) {
    feffectiveR = float(ninfectee)/float(ninfector);
    feffectiveRUncertainty = sqrt( 
                            float(ninfectee)/(float(ninfector)*float(ninfector)) 
                            + TMath::Power(float(ninfectee),2)/TMath::Power(float(ninfector),3) 
                            );
  }

  if (fDebug) fOutputTSV << day << "\t" << fNSusceptible << "\t" << fNExposed << "\t" << fNInfectious 
                    << "\t" <<  fNRecovered << "\t" << fNTraced+fNReported << "\t" 
                    << fNQuarantine.at(day) << "\t" << feffectiveR << "\t" << feffectiveRUncertainty 
                    << "\t" << doublingTime << "\t" << doublingTimeUncertainty << endl;
  if (fDebug) cout << day << "\t" << fNSusceptible << "\t" << fNExposed << "\t" << fNInfectious 
                    << "\t" <<  fNRecovered << "\t" << fNTraced+fNReported << "\t" 
                    << fNQuarantine.at(day) << "\t" << feffectiveR << "\t" << feffectiveRUncertainty << "\t" << doublingTime << endl;
  fDayForROOTTree = day;
  fPopulationLevelInformation->Fill();
                      
}
 
void CVMC::WriteGnuplotScript() {

  fOutputGnuplot.open(Form("plot%s_%d.gnu",fOutputPrefix.c_str(),fRunId));
  fOutputGnuplot << "#gnuplot" << endl;
  fOutputGnuplot << "set term png size 900,900" << endl;
  fOutputGnuplot << Form("set output \"CovidMCPopulation_%d.png\"",fRunId) << endl;

  fOutputGnuplot << "set yrange [0:1]" << endl;
  fOutputGnuplot << Form("set xrange [0:%d]",fLastDayWithPatients+50) << endl;

  fOutputGnuplot << "# (1) day , (2) nSusceptible , (3) nExposed , (4) nInfectious ,  " << endl;
  fOutputGnuplot << "# (5) nRecovered \t , nTraced+nReported , (7) nQuarantine , " << endl;
  fOutputGnuplot << "# (8) R_eff (10 day average) , (9) R_eff uncertainty" << endl;
  fOutputGnuplot << "# (10) doubling time , (11) doubling time uncertainty" << endl;

  fOutputGnuplot << "set key bottom right" << endl;
  fOutputGnuplot << "set grid" << endl;
  fOutputGnuplot << "set xlabel \"\"" << endl;
  fOutputGnuplot << "set ylabel \"Fraction of total population\"" << endl;
  fOutputGnuplot << "set multiplot layout 4,1" << endl;
  fOutputGnuplot << "set bmargin 0" << endl;
  fOutputGnuplot << "set lmargin 10" << endl;
  fOutputGnuplot << "set format x \"\"" << endl;


  fOutputGnuplot << Form("plot \"CovidMCResult_%d.txt\"  using 1:($2/($2+$3+$4+$5)) w l title \"Susceptible\", \\",fRunId) << endl;
  fOutputGnuplot << Form("     \"CovidMCResult_%d.txt\" using 1:($3/($2+$3+$4+$5))  w l title \"Exposed\", \\",fRunId) << endl;
  fOutputGnuplot << Form("     \"CovidMCResult_%d.txt\"  using 1:(($4-$6)/($2+$3+$4+$5)) w l lw 2 title \"Inf(free)\",\\",fRunId) << endl;
  fOutputGnuplot << Form("     \"CovidMCResult_%d.txt\"  using 1:($6/($2+$3+$4+$5)) w l lw 2 title \"Inf(qarant.)\",\\",fRunId) << endl;
  fOutputGnuplot << Form("     \"CovidMCResult_%d.txt\" using 1:($5/($2+$3+$4+$5)) w l title \"Recovered\" ,\\",fRunId) << endl;
  fOutputGnuplot << Form("     \"CovidMCResult_%d.txt\" using 1:(($7-$4)/($2+$3+$4+$5)) w l title \"Healthy(qarant.)\" ",fRunId) << endl;

  fOutputGnuplot << "set ylabel \"R_{eff}\"" << endl;
  fOutputGnuplot << "set xlabel \"\"" << endl;
  fOutputGnuplot << "unset logscale y" << endl;
  fOutputGnuplot << "set tmargin 1" << endl;
  fOutputGnuplot << "set yrange [0:6]" << endl;
  fOutputGnuplot << "set bmargin 0" << endl;
  fOutputGnuplot << Form("plot \"CovidMCResult_%d.txt\" using 1:8:9 w yerrorbars notitle",fRunId) << endl;


  fOutputGnuplot << "set autoscale y" << endl;
  int fitfrom = 0., fitto = 0.;
  float smoothingPeriod = 9;
  float runninaverageA = 0, runninaverageB = 0;
  for (int iday = 0; iday < fLastDayWithPatients-smoothingPeriod; iday++) {
    if (fitfrom < 1. && float(fNumberInfectiousByDay.at(iday))/float(fNPersons) > 0.0005) fitfrom = iday;
    if (iday <= smoothingPeriod) { 
      runninaverageA = runninaverageA +fNumberInfectiousByDay.at(iday); 
      runninaverageB = runninaverageB + fNumberInfectiousByDay.at(iday+smoothingPeriod); 
      }
    if (iday > smoothingPeriod) { 
      runninaverageA = runninaverageA - fNumberInfectiousByDay.at(iday-smoothingPeriod) + fNumberInfectiousByDay.at(iday);
      runninaverageB = runninaverageB - fNumberInfectiousByDay.at(iday) + fNumberInfectiousByDay.at(iday+smoothingPeriod);
      if (runninaverageB*0.6 < runninaverageA ) {fitto = iday; break;}
    }
    fitto = iday;
  }

  fOutputGnuplot << Form("set xrange [1:%d]",fitto+50) << endl;
  try { fOutputGnuplot << Form("set yrange [0:%.3f]",fNumberInfectiousByDay.at(fitto+10)*1.5) << endl; }
  catch (...) { fOutputGnuplot << Form("set yrange [0:%.3f]",fNumberInfectiousByDay.at(fitto)*1.5) << endl; }
  fOutputGnuplot << "f(x, l, a) = a*exp(x*log(2)/l)" << endl;
  fOutputGnuplot << "l = 5" << endl;
  fOutputGnuplot << "a = 17.5" << endl;
  fOutputGnuplot << "set ylabel \"People\"" << endl;
  fOutputGnuplot << Form("fit [%d:%d] f(x,l,a) \"CovidMCResult_%d.txt\" using 1:4 via l,a", fitfrom, fitto, fRunId) << endl;
  fOutputGnuplot << "ti = sprintf(\"T2 = %.2f\", l)" << endl;
  fOutputGnuplot << Form("plot \"CovidMCResult_%d.txt\"  using 1:4 w l lc 4 title \"Infected\", \\", fRunId) << endl;
  fOutputGnuplot << Form("     (x > %d && x < %d) ? f(x,l,a) : 1/0 lc 3 lw 2 title ti",fitfrom, fitto)<< endl;
  fOutputGnuplot.close();
}


void CVMC::WriteROOTFile() {
  fOutputRoot->cd();  
  fhIncubationPeriod->SetLineColor(1);
  fhInfectiousness->SetLineColor(kBlue);
  fhLatentPeriod->SetLineColor(kRed);
  fhIncubationPeriod->Scale(1.0/fhIncubationPeriod->Integral());
  fhIncubationPeriod->Write();
  fhInfectiousness->Scale(fDisease->GetInfectiousness()->Integral(-2,20)/fhInfectiousness->Integral()); // turn into probability
  fhInfectiousness->Write();
  fhLatentPeriod->Scale(1.0/fhLatentPeriod->Integral());
  fhLatentPeriod->Write();
  fhInfectiousnessAgeSymptom->Scale(1.0/fhInfectiousnessAgeSymptom->Integral());
  fhInfectiousnessAgeSymptom->Write();
  fhInfectiousnessAgeNoSymptom->Scale(1.0/fhInfectiousnessAgeNoSymptom->Integral());
  fhInfectiousnessAgeNoSymptom->Write();
  
  fhInfectiousnessDurationNoSymptom->Scale(1.0/fhInfectiousnessDurationNoSymptom->Integral());
  fhInfectiousnessDurationSymptom->Scale(1.0/fhInfectiousnessDurationSymptom->Integral());
    
  fhInfectiousnessDurationSymptom->Write();
  fhInfectiousnessDurationNoSymptom->Write();
  
  fhDaysToQuarantinePostIntervention->Scale(1.0/float(fNExposedTotal));
  fhDaysToQuarantinePreIntervention->Scale(1.0/float(fNExposedTotal));
  fhDaysToReportedPostIntervention->Scale(1.0/float(fNExposedTotal));
  fhDaysToReportedPreIntervention->Scale(1.0/float(fNExposedTotal));
  fhDaysToTestedPostIntervention->Scale(1.0/float(fNExposedTotal));  
  fhDaysToTestedPreIntervention->Scale(1.0/float(fNExposedTotal));
          
  fhDaysToQuarantinePostIntervention->Write();
  fhDaysToQuarantinePreIntervention->Write();  
  
  fhDaysToReportedPostIntervention->Write();  
  fhDaysToReportedPreIntervention->Write();
  
  fhDaysToTestedPostIntervention->Write();      
  fhDaysToTestedPreIntervention->Write();   
  fhNumberInfectedPreIntervention->Write();     
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);  
  TCanvas *tc =new TCanvas("tc","tc");
  TH1F *hdummy = new TH1F("hdummy",";day;probability",33,-3,30);
  hdummy->Draw();
  hdummy->GetYaxis()->SetRangeUser(0.0, 1.0);
    // Draw the test positive probability
  float xvals[40];
  float yvals[40];
  int xycount = 0;
  for (int iday = -3; iday < 30; iday++) {
    xvals[xycount] = iday;
    yvals[xycount] = fDisease->GetTestPositiveProb(iday);
      xycount++;
  }
  TGraph *tcc  = new TGraph(xycount, &xvals[0], &yvals[0]);
  tcc->SetMarkerStyle(2);
  tcc->SetLineColor(1);
  tcc->SetName("TestPositiveProb");
  tcc->Write();
  tcc->Draw("ALPsame");
  
  fDisease->GetIncubationPeriod()->SetLineColor(kGray);
  fhIncubationPeriod->Draw("same");
  fDisease->GetIncubationPeriod()->Draw("same");
  fDisease->GetInfectiousness()->SetLineColor(kBlue-2);

  fhInfectiousness->SetLineWidth(2);
  


  fDisease->GetInfectiousness()->Draw("same");
  fhInfectiousness->Draw("same");
  fhLatentPeriod->Draw("same");
  TLegend *fl = new TLegend();
  fl->AddEntry(fhIncubationPeriod, "Incubation Period" , "l");
  fl->AddEntry(fhInfectiousness, "Infectiousness relative to symptom start" , "l");
  fl->AddEntry(fhLatentPeriod, "Latent Period" , "l");
  fl->AddEntry(tcc, "Test positive probability" , "lp");      
  makeLegend(fl, tc, 1);
  fl->Draw();
  tc->Modified();
  tc->Update();  
  if (fDebug) tc->SaveAs(Form("diagnostic%s_%d.pdf",fOutputPrefix.c_str(),fRunId));
  tc->Write();
  fPopulationLevelInformation->Write();
  
  //
  TTree* settings = new TTree("settings","");
  settings->Branch("index",&fRunId);
  settings->Branch("randomSeed",&fRandomSeed);
  settings->Branch("nPersons",&fNPersons);
  settings->Branch("nDays",   &fNDays);
  settings->Branch("LastDayWithPatients", &fLastDayWithPatients);
  settings->Branch("peopleMetPerDay", &fPeopleMetPerDay);
  settings->Branch("daysInQuarantine",&fDaysInQuarantine);
  settings->Branch("tracingOrder",    &fTracingOrder);
  settings->Branch("startTracingTestingInfectedFraction", &fStartTracingTestingInfectedFraction);
  settings->Branch("startTracingOnDay", &fStartTestingOnDay);
  settings->Branch("daysBackwardTrace", &fDaysBackwardTrace);
  settings->Branch("backwardTracing", &fBackwardTracing);
  settings->Branch("traceUninfected", &fTraceUninfected);
  settings->Branch("SocialDistancingMaxPeople", &fSocialDistancingMaxPeople);
  settings->Branch("SocialDistancingFrom", &fSocialDistancingFrom);
  settings->Branch("SocialDistancingTo ", &fSocialDistancingTo);
  settings->Branch("SocialDistancingFactor", &fSocialDistancingFactor);
  settings->Branch("dTTest", &fdTTest);
  float TestThreshold = fDisease->GetTestThreshold();
  settings->Branch("testThreshold", &TestThreshold);
  settings->Branch("AppProbability", &fAppProbability);
  settings->Branch("ReportingProbability", &fReportingProbability);  
  settings->Branch("tracingEfficiency", &fTracingEfficiency);
  settings->Branch("TracingDelay", &fTracingDelay);
  settings->Branch("startTestingOnDay", &fStartTestingOnDay);
  settings->Branch("daysToTestResult",  &fDaysToTestResult);
  settings->Branch("randomTesting", &fRandomTesting);
  settings->Branch("randomTestingRate", &fRandomTestingRate);
  float symptomProbability      = fDisease->GetSymptomProbability();
  float testPositiveProbability = fDisease->GetTestPositiveProbability();
  float falsePositiveRate       = fDisease->GetFalsePositiveRate();
  float transmissionProbability = fDisease->GetTransmissionProbability();
  settings->Branch("symptomProbability",     &symptomProbability);      
  settings->Branch("testPositiveProbability",&testPositiveProbability);
  float AsymptomaticTransmissionScaling = fDisease->GetAsymptomaticTransmissionScaling();
  settings->Branch("AsymptomaticTransmissionScaling",&AsymptomaticTransmissionScaling);
  settings->Branch("falsePositiveRate",      &falsePositiveRate);
  settings->Branch("transmissionProbability",&transmissionProbability);  
  float incubationGamma = fDisease->GetIncubationPeriod()->GetParameter(0);
  float incubationMu    = fDisease->GetIncubationPeriod()->GetParameter(1);
  float incubationBeta  = fDisease->GetIncubationPeriod()->GetParameter(2);
  float infectionGamma  = fDisease->GetInfectiousness()->GetParameter(0);
  float infectionMu     = fDisease->GetInfectiousness()->GetParameter(1);
  float infectionBeta   = fDisease->GetInfectiousness()->GetParameter(2); 
  settings->Branch("incubationGamma",&incubationGamma); 
  settings->Branch("incubationMu",   &incubationMu); 
  settings->Branch("incubationBeta", &incubationBeta); 
  settings->Branch("infectionGamma", &infectionGamma); 
  settings->Branch("infectionMu",    &infectionMu); 
  settings->Branch("infectionBeta",  &infectionBeta); 
  
  //! Now we want to figure out the R0 with the settings we put  in
  float R00=0., R0e=0.; float R0IfInHospital=0.; float R0NoSymptoms=0.; // figure out R0
  for (int iday = -3; iday < 30; iday++) {
  //! Note: We first infect others, then check if we can report a person. Thus, people stop 
  //! infecting the day after they are reported. Hence it is "iday <= 0" not "iday < 0".
  //! presymptomatics should scale the same as asymptomatics
    if (iday <= 0) R0IfInHospital = R0IfInHospital + AsymptomaticTransmissionScaling * fPeopleMetPerDay*fDisease->GetInfectiousness()->Eval(iday);
    if (iday > 0)  R0NoSymptoms = R0NoSymptoms + fPeopleMetPerDay*fDisease->GetInfectiousness()->Eval(iday);  
  }
  //! R without removing reported symptomatics (also account for possiblity of asymptomatics infecting less
  //! The symptomatics have the same R as the asymptomatics if we don't report them.
  R00 = (R0IfInHospital+R0NoSymptoms) * symptomProbability + (R0IfInHospital+R0NoSymptoms*AsymptomaticTransmissionScaling) * (1. - symptomProbability);
 
  //! R after removing reported symptomatics  
  R0e = symptomProbability *   //! For those that develop symptoms:
        ( R0IfInHospital * fReportingProbability + //! the fraction of symptomatics that becomes reported is in the hospital or in self-quarantine, either way they stop infecting when symptoms start ... 
          (R0IfInHospital+R0NoSymptoms)*(1.-fReportingProbability)  //! but the ones that don't become reported also don't self-quarantine and keep infecting
        )   
        + (1. - symptomProbability) * (R0IfInHospital+R0NoSymptoms*AsymptomaticTransmissionScaling); // people without symptoms do not self quarantine, but they may be less infectious

  settings->Branch("R00",  &R00);
  settings->Branch("R0e",  &R0e);  

  settings->Fill();
  settings->Write();  
    
  fOutputRoot->Close();  
}

// because ROOT can't make a decent legend automatically
void CVMC::makeLegend(TLegend *ll, TPad *pad, int location)
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
  TLegendEntry* entry;  // a pointer to each entry
  int longestLabel = 0;  // count how long the labels are
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
    case 0:   // top right corner aligned with histogram top right corner
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
