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
 * Infrastructure for running the outbreak simulation.
 * 
 * @author Tina Pollmann, Christoph Wiesinger
 * @date 2020
 */
 

//c++
#include "getopt.h"

//external (included in project)
#include "json.hpp"

//root
#include "TSystem.h"

#include "CVMC.h"

using namespace std;

//! accessible via optarg
bool   gDebugMode      = false;
int    gMaxPeopleInDot =   400;
int    gNSimulations   =     1; //! times
int    gIndex          =     0; //! start index for output
int    gRandomSeed     =     0; //! random seed (0=random)

//! accessible via json 
string gOutputPrefix     =  "CovidMCResult";
int    gNPersons         =   400; //! persons
int    gNDays            =   500; //! days
float  gPeopleMetPerDay  =     5; //! persons/day
int    gDaysInQuarantine =    14; //! days
float  gStartTracingTestingInfectedFraction = 0.02; //!fraction of total population
//! .. social distancing
float  gSocialDistancingMaxPeople =    10; //! persons/day
float  gSocialDistancingFrom      = 99999; //! day
float  gSocialDistancingTo        = 99999; //! day
float  gSocialDistancingFactor    =    1.; //! 
//! .. tracing
int    gTracingOrder      =     1; //!      
int    gDaysBackwardTrace =    14; //! days
bool   gBackwardTracing   = true; //!
bool   gTraceUninfected   = true; //!
float  gAppProbability    = 0.447; //!
float  gReportingProbability    = 0.8; //!
float  gTracingEfficiency =    1.; //!
int    gTracingDelay      =     0;
//! .. testing
int    gDaysToTestResult  = 0; //! days
int    gdTTest            = 5; //! days
float  gTestThreshold     = 0.01; //! Infectivity per day at which fast test loses sensitivity to viral RNA
int    gRandomTesting     = false;
float  gRandomTestingRate = 0.01; //! 1/day
//! .. disease
float  gSymptomProbability      =  0.50;  //!
float  gTestPositiveProbability =  0.95;  //!
float  gFalsePositiveRate       =  0.00;  //!
float  gTransmissionProbability =  0.064; //!
float  gAsymptomaticTransmissionScaling = 0.1; //!
float  gIncubationGamma         =  3.06;  //!
float  gIncubationMu            =  0.00;  //!
float  gIncubationBeta          =  2.44;  //!
float  gInfectionGamma          =  2.08;  //!
float  gInfectionMu             = -2.42;  //!
float  gInfectionBeta           =  1.56;  //!

void Usage(const char* const exe);
int  GetOptions(int argc, char** argv);
void ParseJSON(string);

int main(int argc, char** argv)
{
  //! get option input
  const int nOptions = GetOptions(argc, argv);
  if(nOptions < 0 || argc - nOptions < 1) {
    Usage(argv[0]);
    return 1;
  };
  //! get json input
  string input_filename(argv[nOptions]);
  ParseJSON(input_filename.c_str());

  //! initialize and set
  CVMC* sim = new CVMC(gNPersons,gNDays,gAppProbability,gReportingProbability,gOutputPrefix);
  if(gDebugMode) sim->SetDebug();
  //! ... general settings
  sim->SetPeopleMetPerDay(gPeopleMetPerDay);
  sim->SetDaysInQuarantine(gDaysInQuarantine);
  //! ... social distancing settings
  sim->SetSocialDistancingMaxPeople(gSocialDistancingMaxPeople);
  sim->SetSocialDistancingFrom(gSocialDistancingFrom); 
  sim->SetSocialDistancingTo(gSocialDistancingTo);      
  sim->SetSocialDistancingFactor(gSocialDistancingFactor);
  //! ... tracing settings
  sim->SetTracingOrder(gTracingOrder);
  sim->SetStartTracingTestingInfectedFraction(gStartTracingTestingInfectedFraction);
  sim->SetDaysBackwardTrace(gDaysBackwardTrace);
  sim->SetBackwardTracing(gBackwardTracing);
  sim->SetTraceUninfected(gTraceUninfected);
  sim->SetTracingEfficiency(gTracingEfficiency);
  sim->SetTracingDelay(gTracingDelay);
  //! ... testing settings
  sim->SetDaysToTestResult(gDaysToTestResult);
  sim->SetDTTest(gdTTest);
  sim->SetRandomTesting(gRandomTesting);
  sim->SetRandomTestingRate(gRandomTestingRate);
  //! ... disease settings
  sim->GetDisease()->SetSymptomProbability(gSymptomProbability);          
  sim->GetDisease()->SetTestPositiveProbability(gTestPositiveProbability);   
  sim->GetDisease()->SetFalsePositiveRate(gFalsePositiveRate);
  sim->GetDisease()->SetIncubationParameters(gIncubationGamma, gIncubationMu, gIncubationBeta);    
  sim->GetDisease()->SetTransmissionProbability(gTransmissionProbability);
  sim->GetDisease()->SetAsymptomaticTransmissionScaling(gAsymptomaticTransmissionScaling); 
  sim->GetDisease()->SetInfectiousnessParameters(gInfectionGamma,gInfectionMu,gInfectionBeta);
  sim->GetDisease()->SetTestThreshold(gTestThreshold);
  //! ... and run it repeatedly
  for(int irun=gIndex;irun<gIndex+gNSimulations;irun++)
    sim->Run(irun,gRandomSeed);

  delete sim;
  //! Add all the output files together if more than one simulation was run, and then delete the individual files.
  if (gNSimulations > 1) {
    string filenames = Form("%s_*.root",gOutputPrefix.c_str());
    string filenamesum = Form("Sum_%s.root",gOutputPrefix.c_str());
    gSystem->Exec(Form("hadd -f %s %s",filenamesum.c_str(),filenames.c_str()));
    gSystem->Exec(Form("rm %s",filenames.c_str()));  
  }
  return 0;
}

void Usage(const char* const exe)
{
   cout << "usage: " << exe << " [<options>] json \n"
    "\n"
    " options:  "/*-o (or --output):   output prefix                 (default: " << gOutputPrefix   << ")\n" 
    "        */"-n (or --nsim):     simulation number             (default: " << gNSimulations   << ")\n"     
    "           -i (or --index):    start index                   (default: " << gIndex          << ")\n"    
    "           -s (or --seed):     random number seed            (default: " << gRandomSeed     << ")\n"  
    "           -d (or --debug):    run with increased verbostiy \n" 
    "           -m (or --maxdots):  maximum people in dotfile     (default: " << gMaxPeopleInDot << ")\n"  
   << endl;
}

int GetOptions(int argc, char** argv)
{
   static struct option long_options[] = {
     {"nsim",    required_argument, 0,'n'}, 
     {"index",   required_argument, 0,'i'}, 
     {"seed",    required_argument, 0,'s'}, 
     {"debug",   no_argument,       0,'d'},  
     {"maxdots", required_argument, 0,'m'},         
     {"help",    no_argument,       0,'h'},
     {0, 0, 0, 0}
   };

   int option_index = 0;
   int c;
   while ((c = getopt_long (argc, argv,
    ":n:i:s:dm:h",
    long_options, &option_index)) != -1) {
     switch (c) {
       case 'n': gNSimulations   = stoi(optarg); break;       
       case 'i': gIndex          = stoi(optarg); break; 
       case 's': gRandomSeed     = stoi(optarg); break;      
       case 'd': gDebugMode      = true;         break;
       case 'm': gMaxPeopleInDot = stoi(optarg); break;     
       case 'h': return -2;
       default:  return -2;
     }
   }
   return optind;
}

void ParseJSON(std::string inputfile)
{
  using json = nlohmann::json;
  std::ifstream i(inputfile.c_str());
  json j;
  i >> j;

  //! general
  gOutputPrefix     = j["OutputPrefix"];
  gNPersons         = j["NPersons"];
  gNDays            = j["NDays"];
  gPeopleMetPerDay  = j["PeopleMetPerDay"];
  gDaysInQuarantine = j["DaysInQuarantine"];
  gStartTracingTestingInfectedFraction = j["StartTracingTestingInfectedFraction"];
  //! social distancing
  gSocialDistancingMaxPeople = j["SocialDistancingMaxPeople"];
  gSocialDistancingFrom      = j["SocialDistancingFrom"];
  gSocialDistancingTo        = j["SocialDistancingTo"];
  gSocialDistancingFactor    = j["SocialDistancingFactor"];
  //! tracing
  gTracingOrder      = j["TracingOrder"];
  gDaysBackwardTrace = j["DaysBackwardTrace"];
  gBackwardTracing   = j["BackwardTracing"];
  gTraceUninfected   = j["TraceUninfected"];
  gAppProbability    = j["AppProbability"];
  gReportingProbability    = j["ReportingProbability"]; 
  gTracingEfficiency = j["tracingEfficiency"];
  gTracingDelay      = j["TracingDelay"];
  //! testing
  gDaysToTestResult  = j["DaysToTestResult"];
  gdTTest            = j["dTTest"];
  gTestThreshold     = j["TestThreshold"];
  gRandomTesting     = j["RandomTesting"];
  gRandomTestingRate = j["RandomTestingRate"];
  //! disease
  gSymptomProbability      = j["SymptomProbability"];
  gTestPositiveProbability = j["TestPositiveProbability"];
  gFalsePositiveRate       = j["FalsePositiveRate"];
  gTransmissionProbability = j["TransmissionProbability"];
  gAsymptomaticTransmissionScaling = j["AsymptomaticTransmissionScaling"];
  gIncubationGamma         = j["IncubationPeriod"]["gamma"];
  gIncubationMu            = j["IncubationPeriod"]["mu"];
  gIncubationBeta          = j["IncubationPeriod"]["beta"];
  gInfectionGamma          = j["Infectivity"]["gamma"];
  gInfectionMu             = j["Infectivity"]["mu"];
  gInfectionBeta           = j["Infectivity"]["beta"];
}
