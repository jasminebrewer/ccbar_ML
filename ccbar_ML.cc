#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "functions_ccbar.h"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;
using namespace std;



int main(int argc, char* argv[])
{

  string header = string("jets");
  bool D0decay=false; // don't allow the D0 to decay
  int doD0decay = atoi(argv[2]);
  string decaystring = string("nodecay");
  if (doD0decay==1) {
    D0decay=true;
    decaystring = string("decay");
  }
  
  string paramfile_name = string("params.dat");
  Pythia pythia;
  struct trackCuts track_cuts = {0};
  
  int maxnevents = initialize_pythia(paramfile_name, pythia, track_cuts, D0decay);

  cout << "Using R = " << track_cuts.jetR << std::endl;
  ofstream hadronfile(header + "_hadron_" + decaystring + "_" + argv[1] + ".dat");
  ofstream partonfile(header + "_parton_" + decaystring + "_" + argv[1] + ".dat");

  if (!hadronfile.is_open() || !partonfile.is_open())
    cout << "Unable to open file";

  // define cuts for tracks, jets, and heavy flavor particles
  track_cuts.JetEtaMin = -track_cuts.trackEtaCut + track_cuts.jetR; // jet eta range
  track_cuts.JetEtaMax = -track_cuts.JetEtaMin;
  track_cuts.HFLowPtCut = 5; // only tag heavy-flavor particles with pT > 5 GeV

  vector<fastjet::PseudoJet> final_hadrons, final_partons; // only final state particles (at hadron and parton-level, respectively)
  vector<fastjet::PseudoJet> tagged_hadrons, tagged_partons;
  vector<int> hadron_ids{constants::D0, constants::D0BAR};
  vector<int> parton_ids{constants::CHARM, constants::ANTICHARM};

  // Begin event loop. Generate event. Skip if error. List first one. 
  for (int iEvent = 0; iEvent < maxnevents; ++iEvent) {
    	  
    // ignore events where pythia aborts
    if (!pythia.next())
      continue;
    Event& event = pythia.event;

    // initialize per event
    final_hadrons.clear();
    final_partons.clear();	 
    tagged_hadrons.clear();
    tagged_partons.clear();

    // read in particles from the event. Save final-state hadrons, final-state partons, and "tagged_partons" and "tagged hadrons" corresponding to particles with the particle ids contained in hadron_ids and parton_ids
    read_event( event, final_hadrons, final_partons, tagged_hadrons, tagged_partons, hadron_ids, parton_ids, track_cuts );

    // identify g->cc splitting
    int flagged=0;
    if ( tagged_partons.size()==2 ) {
      
      splitting tagged_split = find_common_splitting(event, tagged_partons[0].user_info<MyUserInfo>().global_index(), tagged_partons[1].user_info<MyUserInfo>().global_index());
      
      if ( valid_splitting(tagged_split) ) flagged=1;
    }
    
    
    // cluster jets regardless of whether they have a ccbar pair. Later, searcj
    cluster_jets( event, final_partons, tagged_partons, parton_ids, partonfile, track_cuts, iEvent, flagged);
    cluster_jets( event, final_hadrons, tagged_hadrons, hadron_ids, hadronfile, track_cuts, iEvent, flagged);
    
  } // End of event loop.

  return 0;
}
