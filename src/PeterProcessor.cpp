#include "PeterProcessor.hpp"
#include <limits>
#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"

#include "marlin/Global.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/VerbosityLevels.h"
#include "marlinutil/GeometryUtil.h"
#include "marlinutil/DDMarlinCED.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/Randomize.h"
#include "HelixClass.h"
#include "marlinutil/CalorimeterHitType.h"
#include <unordered_map>

using namespace EVENT;
using dd4hep::rec::Vector3D;
using std::cout, std::endl, std::vector, std::string;
using UTIL::LCRelationNavigator;
PeterProcessor aPeterProcessor ;

PeterProcessor::PeterProcessor() : marlin::Processor("PeterProcessor"), EventDisplayer(this) {}


void PeterProcessor::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);                                                
    DDMarlinCED::init(this);                                                                                                                                                                                                                                                                                                                                                               
    _AngSpecfile.open ("/nfs/dust/ilc/user/mckeownp/Angular_spectra/ILD_ECAL_Angular_Spectra_data/Tau_pairs_2f_leptonic_eL_pR_photons_angular_spectra_MC_only_CORRECT.csv");
    _AngSpecfile << "Barrel Flag, Endcap Flag, Sim Flag, Overlay Flag, Phi angle to Norm, Theta angle to Norm, Energy (GeV), Parent PDG\n";
    _invMassfile.open ("/nfs/dust/ilc/user/mckeownp/Angular_spectra/ILD_ECAL_Angular_Spectra_data/Tau_pairs_2f_leptonic_eL_pR_photons_angular_spectra_MC_only_INV_MASS.csv");
    _invMassfile << "tau1_PDG, tau1_E, tau1_px, tau1_py, tau1_pz, tau2_PDG, tau2_E, tau2_px, tau2_py, tau2_pz\n";

}
    

void draw(EVENT::MCParticle* mc, std::vector<EVENT::Cluster*> clusters){
    int type = 0;
    int size = 5;
    unsigned long color = 0xd91c2f;

    Vector3D mcPos (mc->getVertex());
    Vector3D mcMom (mc->getMomentum());
    mcMom = 3000.*(mcMom.unit());
    ced_line( mcPos.x(), mcPos.y(), mcPos.z(), mcMom.x(), mcMom.y(), mcMom.z() , type , size, color);

    color = 0x09ed1c;
    for(auto* cluster : clusters){
        for( auto* hit : cluster->getCalorimeterHits() ){
            Vector3D hitPos (hit->getPosition());
            ced_hit_ID(hitPos.x(), hitPos.y(), hitPos.z(), type, 0, size, color, 0);        
        }
    }
}

void PeterProcessor::processEvent(EVENT::LCEvent * event){
    CLHEP::RandGauss::setTheSeed( marlin::Global::EVENTSEEDER->getSeed(this) ); 
    //++_nEvent;
    cout <<endl;
    cout<<"***********************************************"<<endl;                                       
    cout<<"****************Event "<<(++_nEvent)<<"**************************"<<endl;                     
    cout<<"***********************************************"<<endl;                                       
                                                                           

    //LCCollection* pfos = event->getCollection("PandoraPFOs");
    LCCollection* mc_col = event->getCollection("MCParticle");
    LCRelationNavigator nav( event->getCollection("RecoMCTruthLink") );
    // Get barrel and endcap hits, ignore Si/Sc technology
    LCCollection* ECAL_Barrel_col = event->getCollection("EcalBarrelCollection");
    LCCollection* ECAL_Endcap_col = event->getCollection("EcalEndcapsCollection");

    // If the collection type is wrong you end up with a null pointer here.    
    // Always check it !
    
    if(nullptr == mc_col) {
        streamlog_out(ERROR) << "Wrong object type in collection '" << "MCParticle" << "'" << std::endl;
      }
    

    if(nullptr == ECAL_Barrel_col) {
        streamlog_out(ERROR) << "Wrong object type in collection '" << "EcalBarrelCollection" << "'" << std::endl;
      }

    if(nullptr == ECAL_Endcap_col) {
        streamlog_out(ERROR) << "Wrong object type in collection '" << "ECAL_Endcap_Hits" << "'" << std::endl;
      }  


    // Check if this event contains taus- perhaps not the most efficient, but it does the job :-)
    bool tau_evt = false;

    // Create vector for taus
    std::vector<MCParticle*> taus_vec; 

    for (int i=0; i<mc_col->getNumberOfElements(); ++i){

        // based on assumption that the first two taus are from initial hard scattering
        MCParticle* mc_part = static_cast <MCParticle*> ( mc_col->getElementAt(i) );
        if( abs(mc_part->getPDG()) == 15 ) {
	
       	   std::cout<<"Looping Tau check "<<std::endl;   
	
            int vec_size = static_cast<int>(taus_vec.size());

            if ( taus_vec.empty() == true ){
		std::cout<<"Looping Tau check: about to add to vector "<<std::endl;
                taus_vec.push_back( mc_part );
		std::cout<<"Looping Tau check: about to added to vector "<<std::endl;

            }

            else if ( vec_size == 1 ){
                taus_vec.push_back( mc_part );
            }

            else if ( vec_size == 2 ) break;

            else continue;
	    std::cout<<"Done Tau check "<<std::endl;

        }
    }


    // have to check if length of tau vector is equal to 2 before this check. Otherwise, not tau event

    int vec_size = static_cast<int>(taus_vec.size());

    bool matching_tau_parents = false;

    if ( vec_size == 2 ){
	std::cout<<"Taus found, starting cross-checks "<<std::endl;


        // get parents of taus, and check if they are the same

        std::vector<MCParticle*> tau_1_parents = taus_vec[0] -> getParents(); //auto
        std::vector<MCParticle*> tau_2_parents = taus_vec[1] -> getParents(); //auto

        /*
        MCParticle* tau_1_parent = static_cast <MCParticle*> ( taus_vec[0] );
        MCParticle* tau_2_parent = static_cast <MCParticle*> ( taus_vec[1] );
        */

        // check if parents of MC are the same- -> then it is a tau event.

        int parent_size_tau_1 = static_cast<int>(tau_1_parents.size());
        int parent_size_tau_2 = static_cast<int>(tau_2_parents.size());
        if(parent_size_tau_1 > 1) {std::cout<<"Tau1 with more than one parent!"<<std::endl;}
        if(parent_size_tau_2 > 1) {std::cout<<"Tau2 with more than one parent!"<<std::endl;}

        if (parent_size_tau_1 == parent_size_tau_2){
            // loop over parents, and check if they are the same
            for (int j=0; j<parent_size_tau_1; ++j){
                int Parent_type_tau_1 = tau_1_parents[j]->getPDG();
                int Parent_type_tau_2 = tau_2_parents[j]->getPDG();

                if ( Parent_type_tau_1 ==  Parent_type_tau_2 ){
                        matching_tau_parents = true ;
                    }
                else { matching_tau_parents = false ;}
            }
        }
        else {std::cout<<"WARNING: Taus with different numbers of parents!"<<std::endl;}
	std::cout<<"Tau crosscheck done "<<std::endl;
    }

    // if the previous checks have passed, this is an event with tau pairs- write to file
    if ( matching_tau_parents == true ){
        std::cout<<"Write taus to file "<<std::endl;
	tau_evt = true;
        float tau1_E = taus_vec[0]->getEnergy();
        float tau2_E = taus_vec[1]->getEnergy();
        Vector3D mom_tau_1( taus_vec[0]->getMomentum() );
        Vector3D mom_tau_2( taus_vec[1]->getMomentum() );
        float tau1_px = mom_tau_1.x();
        float tau1_py = mom_tau_1.y();
        float tau1_pz = mom_tau_1.z();
        float tau2_px = mom_tau_2.x();
        float tau2_py = mom_tau_2.y();
        float tau2_pz = mom_tau_2.z();
        _invMassfile << std::scientific << taus_vec[0]->getPDG() << ",";
        _invMassfile << std::scientific << tau1_E << ",";
        _invMassfile << std::scientific << tau1_px << ",";
        _invMassfile << std::scientific << tau1_py << ",";
        _invMassfile << std::scientific << tau1_pz << ",";
        _invMassfile << std::scientific << taus_vec[1]->getPDG() << ",";
        _invMassfile << std::scientific << tau2_E << ",";
        _invMassfile << std::scientific << tau2_px << ",";
        _invMassfile << std::scientific << tau2_py << ",";
        _invMassfile << std::scientific << tau2_pz << "\n";
    	std::cout<<"Done write taus to file "<<std::endl;
    }




    /*
    for (int i=0; i<mc_col->getNumberOfElements(); ++i){

        MCParticle* mc_part = static_cast <MCParticle*> ( mc_col->getElementAt(i) );
        if( abs(mc_part->getPDG()) == 15 ) {
                tau_evt = true;
                std::cout<<endl<<"****************EVENT: "<<_nEvent<<" Contains Taus ****************"<<std::endl;
                continue;
            }
    }
    */

    if (tau_evt == true){

        // Create (unique) map for mc particles that hit the calorimeter
        struct CaloInfo{
            int n_barrel {0};
            int n_forwardEndcap {0};
            int n_backwardEndcap {0};
            float E_hit_barrel {0};
            float E_hit_forwardEndcap {0};
            float E_hit_backwardEndcap {0};
            int nhits {0};
        };

        std::unordered_map<MCParticle*, CaloInfo> mcCalo;
        
        // Loop over SimCaloHits
        // First the Barrel Collection
        for (int i=0; i<ECAL_Barrel_col->getNumberOfElements(); ++i){
            
            SimCalorimeterHit* sim_barrel_hit = static_cast <SimCalorimeterHit*> ( ECAL_Barrel_col->getElementAt(i) );
            
            // getParticleCont has to take an integer i which is ith contribution to hit- loop over contributions
            // and get MCParticle responsible for triggering shower
            for (int j=0; j<sim_barrel_hit->getNMCContributions(); ++j){
            
                MCParticle* barrel_mc = static_cast <MCParticle*> ( sim_barrel_hit-> getParticleCont(j) );

                if ( barrel_mc -> isBackscatter() == true) continue;

                // add energy of hit to energy deposited in barrel

                mcCalo[barrel_mc].E_hit_barrel = mcCalo[barrel_mc].E_hit_barrel + sim_barrel_hit -> getEnergy();

                mcCalo[barrel_mc].n_barrel ++;
                mcCalo[barrel_mc].nhits ++;
            }
        }

        // Second the Endcap Collection
        for (int i=0; i<ECAL_Endcap_col->getNumberOfElements(); ++i){
            
            SimCalorimeterHit* sim_endcap_hit = static_cast <SimCalorimeterHit*> ( ECAL_Endcap_col->getElementAt(i) );
            
            // getParticleCont has to take an integer i which is ith contribution to hit- loop over contributions
            // and get MCParticle responsible for triggering shower
            for (int j=0; j<sim_endcap_hit->getNMCContributions(); ++j){
            
                MCParticle* endcap_mc = static_cast <MCParticle*> ( sim_endcap_hit-> getParticleCont(j) );

                if ( endcap_mc -> isBackscatter() ==true ) continue;

                // check if forward or backward endcap
                // Use getMomentumAtEndpoint: getMomentum only gives momentum at production vertex
                //Vector3D endcap_mom( endcap_mc->getMomentumAtEndpoint() );
                Vector3D endcap_mom( endcap_mc->getMomentum() );

                double theta_endcap = endcap_mom.theta();
                        
                if ( theta_endcap < M_PI/2 )
                    { mcCalo[endcap_mc].n_forwardEndcap ++; 
                      mcCalo[endcap_mc].E_hit_forwardEndcap = mcCalo[endcap_mc].E_hit_forwardEndcap + sim_endcap_hit -> getEnergy();
                    }
                else
                    { mcCalo[endcap_mc].n_backwardEndcap ++; 
                      mcCalo[endcap_mc].E_hit_backwardEndcap =  mcCalo[endcap_mc].E_hit_backwardEndcap + sim_endcap_hit -> getEnergy();
                    }


                mcCalo[endcap_mc].nhits ++;
            }

        }

        // loop over all mc particles in the map
        for( const auto& [mc, caloinfo]: mcCalo ){
            
            // Filter on your desired PDGs here
            // In this case, just Photons
            if (mc->getPDG() != 22) continue;

            // Use getMomentumAtEndpoint: getMomentum only gives momentum at production vertex
            //Vector3D mom( mc->getMomentumAtEndpoint() );
            Vector3D mom( mc->getMomentum() );

            Vector3D mom_pt( mom.x(), mom.y(), 0.0 );

            bool Sim_flag = mc -> isCreatedInSimulation();
            bool Overlay_flag = mc -> isOverlay();

            std::vector<MCParticle*> Parents = mc -> getParents();  

            std::cout<<"***PHOTON MOMENTUM***"<<std::endl;
            std::cout<<mom<<std::endl;
            /*
            // cut on number of hits in the ecal to avoid backscatter
            if ( caloinfo.nhits <=50. ){
                continue;
            }
            */

            // If you don't want to look at events in the corners between endcap and barrel
            // you can exclude here
            // Here showers in multiple parts of the detector are include- beware spooky things could be happening
            if ( caloinfo.n_barrel != 0 && (caloinfo.n_forwardEndcap != 0 ||  caloinfo.n_backwardEndcap != 0) ){

                // calcualte fraction of shower energy deposited in each of barrel, forward and backward Endcap
                float tot_E_shower = caloinfo.E_hit_barrel + caloinfo.E_hit_backwardEndcap + caloinfo.E_hit_forwardEndcap;

                float E_frac_barrel = ((caloinfo.E_hit_barrel)/tot_E_shower)*100.;

                float E_frac_forward = ((caloinfo.E_hit_forwardEndcap)/tot_E_shower)*100.;

                float E_frac_backward = ((caloinfo.E_hit_backwardEndcap)/tot_E_shower)*100.;

                std::cout<<"E_frac_barrel "<<E_frac_barrel<<std::endl;
                std::cout<<"E_frac_forward "<<E_frac_forward<<std::endl;
                std::cout<<"E_frac_backward "<<E_frac_backward<<std::endl;

                if ( E_frac_barrel > 80. ){
                    // consider barrel
                    BarrelFunc( mom, mom_pt, Sim_flag, Overlay_flag, Parents,mc, caloinfo );
                }

                else if ( E_frac_forward > 80. || E_frac_backward > 80.){
                    // consider endcap
                    EndcapFunc( mom, Sim_flag, Overlay_flag, Parents, mc, caloinfo );
                }

                else continue; // most likely in a corner, which we don't want in any case
                
                /*
                int Barrel_flag = 1; // in barrel ...
                int Endcap_flag = 1; // and in endcap

               // in this case, don't record anything for angles, but do record energy
                double energy_barrel = mc->getEnergy();

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                // set Phi and Theta to zero for particle showering in both barrel and endcap
                // safe because this would be directly down the beam axis, so could never hit
                // the calorimeter

                double phi_angle = 0.;
                double theta_angle = 0.;

                // Write to file
                writeOutput(Barrel_flag, Endcap_flag, Sim_flag, Overlay_flag, phi_angle, theta_angle, energy_barrel, Parent_type, _AngSpecfile);
           
                */
            }

            // Barrel only
            if ( caloinfo.n_barrel != 0 && (caloinfo.n_forwardEndcap == 0 && caloinfo.n_backwardEndcap == 0) ){
                BarrelFunc( mom, mom_pt, Sim_flag, Overlay_flag, Parents, mc, caloinfo );
                

                /*
                Vector3D ecalNorm = getBarrelNorm( mom.phi() );

                //double phi_angle_barrel = std::acos( mom.unit()*ecalNorm );
        
                // Phi angle is defined in the transverse plane!
                double phi_angle_barrel = std::acos( mom_pt.unit()*ecalNorm ); 

                //std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double theta_angle_barrel = mom.theta();

                if(theta_angle_barrel > M_PI/2){ theta_angle_barrel = M_PI - theta_angle_barrel; }
                else{
                    theta_angle_barrel = theta_angle_barrel;
                }

                //std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double energy_barrel = mc->getEnergy();
                //std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 1; // In barrel
                int Endcap_flag = 0; // not in endcap

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                int nhits_barrel = caloinfo.nhits;

                // Write to file
                writeOutput(Barrel_flag, Endcap_flag, Sim_flag, Overlay_flag, phi_angle_barrel, theta_angle_barrel, energy_barrel, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if( energy_barrel > 100. ){
                    std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                    std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                }


                std::vector<MCParticle*> Parents_parents = Parents[0] -> getParents();


                int parent_parent_size = static_cast<int>(Parents_parents.size());
                if(parent_parent_size > 1) {std::cout<<"More than one parent of parent!"<<std::endl;}

                if (Parents_parents.empty() == true) continue;  // skip the case of empty parent list

                int Parent_parent_type = Parents_parents[0]->getPDG();

                // Print out if parent was a pi0, which itself came from a tau...
                if( abs(Parent_type) == 111 && Overlay_flag == 0 && Sim_flag == 0){
                    std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Parent's parent PDG: "<<Parent_parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                    std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                }
                */

            }

            // Endcap
            else if ( caloinfo.n_barrel == 0 && (caloinfo.n_forwardEndcap != 0 ||  caloinfo.n_backwardEndcap != 0) ){
                EndcapFunc( mom, Sim_flag, Overlay_flag, Parents, mc, caloinfo );
                /*
                //same for endcap 
                double phi_angle_endcap = mom.phi();
                //std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;            
                double theta_angle_endcap = mom.theta();
                
                if(theta_angle_endcap > M_PI/2){ theta_angle_endcap = M_PI/2. -( M_PI - theta_angle_endcap); }
                else{
                    theta_angle_endcap = M_PI/2. - theta_angle_endcap;
                }
                //std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                double energy_endcap = mc->getEnergy();
                //std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 0; // not in barrel
                int Endcap_flag = 1; // In endcap

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                int nhits_endcap = caloinfo.nhits;

                // Write to file
                writeOutput(Barrel_flag, Endcap_flag, Sim_flag, Overlay_flag, phi_angle_endcap, theta_angle_endcap, energy_endcap, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if(energy_endcap > 100.){
                    std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_endcap<<std::endl;
                    std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                }


                std::vector<MCParticle*> Parents_parents = Parents[0] -> getParents();


                int parent_parent_size = static_cast<int>(Parents_parents.size());
                if(parent_parent_size > 1) {std::cout<<"More than one parent of parent!"<<std::endl;}


                if (Parents_parents.empty() == true) continue;  // skip the case of empty parent list

                int Parent_parent_type = Parents_parents[0]->getPDG();


                // Print out if parent was a pi0, which itself came from a tau...
                if( abs(Parent_type) == 111 &&  abs(Parent_parent_type) == 15 ){
                    std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Parent's parent PDG: "<<Parent_parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_endcap<<std::endl;
                    std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                }
                */


            }

            else continue;
        }

    }

}

    /*

        for (int i=0; i<mc_col->getNumberOfElements(); ++i){ //i<pfos->getNumberOfElements(); ++i){
            std::cout<<endl<<"****************EVENT: "<<_nEvent<<"    PFO: "<<i+1<<"****************"<<std::endl;
            //ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
            //MCParticle* mc = getLinkedMCParticle(pfo, nav);
            
            MCParticle* mc = static_cast <MCParticle*> ( mc_col->getElementAt(i) );

            //auto clusters = pfo->getClusters();
            if (mc->getPDG() != 22) continue;
            //if (clusters.size() != 1) std::cout<<" BEWAAARE, SOMETHING IS SPPOOKY HERE. EVENT: "<<_nEvent<<"    PFO: "<<i+1<<std::endl;
            //if( hasEndcapHits(clusters) ) continue;

            Vector3D mom( mc->getMomentum() );

            Vector3D mom_pt( mom.x(), mom.y(), 0.0 );

            bool Sim_flag = mc -> isCreatedInSimulation();
            bool Overlay_flag = mc -> isOverlay();

            std::vector<MCParticle*> Parents = mc -> getParents();  

            std::cout<<"***PHOTON MOMENTUM***"<<std::endl;
            std::cout<<mom<<std::endl;

            //Barrel
            if( hasEndcapHits(clusters) == false &&  hasBarrelHits(clusters) == true ){
                //assuming photon doesn't change flight direction (it really shouldn't...)
                Vector3D ecalNorm = getBarrelNorm( mom.phi() );
                //std::cout<<"***ECAL NORM***"<<std::endl;
                //std::cout<<ecalNorm<<std::endl;


                //double phi_angle_barrel = std::acos( mom.unit()*ecalNorm );
        
                // Phi angle is defined in the transverse plane!
                double phi_angle_barrel = std::acos( mom_pt.unit()*ecalNorm ); 

                //std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double theta_angle_barrel = mom.theta();

                if(theta_angle_barrel > M_PI/2){ theta_angle_barrel = M_PI - theta_angle_barrel; }
                else{
                    theta_angle_barrel = theta_angle_barrel;
                }

                //std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double energy_barrel = mc->getEnergy();
                //std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 1; // In barrel

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                // Write to file
                writeOutput(Barrel_flag, Sim_flag, Overlay_flag, phi_angle_barrel, theta_angle_barrel, energy_barrel, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if(energy_barrel > 100.){
                    std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                }

                
                //AngSpecfile << Barrel_flag << ",";
                //AngSpecfile << Sim_flag << ",";
                //AngSpecfile << Overlay_flag << ",";
                //AngSpecfile << phi_angle_barrel << ",";
                //AngSpecfile << energy_barrel << ",";
                //AngSpecfile << Parent_type << ",";
                

            }
            // Endcap
            else if( hasEndcapHits(clusters) == true &&  hasBarrelHits(clusters) == false  ){
                //same for endcap 
                double phi_angle_endcap = mom.phi();
                //std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;            
                double theta_angle_endcap = mom.theta();
                
                if(theta_angle_endcap > M_PI/2){ theta_angle_endcap = M_PI/2. -( M_PI - theta_angle_endcap); }
                else{
                    theta_angle_endcap = M_PI/2. - theta_angle_endcap;
                }
                //std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                double energy_endcap = mc->getEnergy();
                //std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 0; // In endcap

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;


                // Write to file
                writeOutput(Barrel_flag, Sim_flag, Overlay_flag, phi_angle_endcap, theta_angle_endcap, energy_endcap, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if(energy_barrel > 100.){
                    std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                }

                
                //AngSpecfile << Barrel_flag << ",";
                //AngSpecfile << Sim_flag << ",";
                //AngSpecfile << Overlay_flag << ",";
                //AngSpecfile << theta_angle_endcap << ",";
                //AngSpecfile << energy_endcap << ",";
                //AngSpecfile << Parent_type << "\n";
               

            }
            else continue;
            //drawDisplay(this, event, draw, mc, clusters);
        }

   }
   else continue;
}

*/

void PeterProcessor::writeOutput(int Barrel_flag, int Endcap_flag, bool Sim_flag, bool Overlay_flag, double Phi_Angle_Norm, double Theta_Angle_Norm, double Energy, double Parent_type, std::ofstream& AngSpecfile){
    AngSpecfile << std::scientific << Barrel_flag << ",";
    AngSpecfile << std::scientific <<  Endcap_flag << ",";
    AngSpecfile << std::scientific << Sim_flag << ",";
    AngSpecfile << std::scientific << Overlay_flag << ",";
    AngSpecfile << std::scientific << Phi_Angle_Norm << ",";
    AngSpecfile << std::scientific << Theta_Angle_Norm << ",";
    AngSpecfile << std::scientific << Energy << ",";
    AngSpecfile << std::scientific << Parent_type << "\n";
}

/*
EVENT::MCParticle* PeterProcessor::getLinkedMCParticle(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav){
    auto objects = nav.getRelatedToObjects(pfo);
    auto weights = nav.getRelatedToWeights(pfo);
    if ( objects.empty() ) return nullptr;
    // get highest track weight, if it is 0, get highest cluster weight
    int max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)%10000)/1000. < (int(rhs)%10000)/1000.;}) - weights.begin();
    if ( ( int(weights[max_i])%10000 )/1000. == 0 ){
        max_i = std::max_element(weights.begin(), weights.end(), [](float lhs, float rhs){return (int(lhs)/10000)/1000. < (int(rhs)/10000)/1000.;}) - weights.begin();
    }
    return static_cast<MCParticle*> (objects[max_i]);
}
*/

void PeterProcessor::BarrelFunc( Vector3D mom, Vector3D mom_pt, bool Sim_flag, bool Overlay_flag, std::vector<MCParticle*> Parents, const auto& mc, const auto& caloinfo ){

                Vector3D ecalNorm = getBarrelNorm( mom.phi() );

                //double phi_angle_barrel = std::acos( mom.unit()*ecalNorm );
        
                // Phi angle is defined in the transverse plane!
                double phi_angle_barrel = std::acos( mom_pt.unit()*ecalNorm ); 

                //std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double theta_angle_barrel = mom.theta();

                if(theta_angle_barrel > M_PI/2){ theta_angle_barrel = M_PI - theta_angle_barrel; }
                else{
                    theta_angle_barrel = theta_angle_barrel;
                }

                int nhits_barrel = caloinfo.nhits;

                // from geometry, exclude very low angle photons- these must come from backscatter
                if ( theta_angle_barrel*180./M_PI <= 15.0 && Overlay_flag == 0 && Sim_flag == 0){
                    std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                    return;
                }
                else if ( theta_angle_barrel*180./M_PI <= 10.0 && (Overlay_flag == 1 || Sim_flag == 1) ){
                    std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                    return;
                }

                //std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                double energy_barrel = mc->getEnergy();
                //std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 1; // In barrel
                int Endcap_flag = 0; // not in endcap

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                // Write to file
                writeOutput(Barrel_flag, Endcap_flag, Sim_flag, Overlay_flag, phi_angle_barrel, theta_angle_barrel, energy_barrel, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if( energy_barrel > 100. ){
                    std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                    std::cout<<"caloinfo.barrel: "<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap: "<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap: "<<caloinfo.n_backwardEndcap<<std::endl;
                }


                std::vector<MCParticle*> Parents_parents = Parents[0] -> getParents();


                int parent_parent_size = static_cast<int>(Parents_parents.size());
                if(parent_parent_size > 1) {std::cout<<"More than one parent of parent!"<<std::endl;}

                // skip the case of empty parent list
                if (Parents_parents.empty() == false){  

                    int Parent_parent_type = Parents_parents[0]->getPDG();

                    // Print out if parent was a pi0, which itself came from a tau...
                    if( abs(Parent_type) == 111 && Overlay_flag == 0 && Sim_flag == 0){
                        std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                        std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                        std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                        std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                        std::cout<<"Parent's parent PDG: "<<Parent_parent_type<<std::endl;
                        std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                        std::cout<<"caloinfo.barrel: "<<caloinfo.n_barrel<<std::endl;
                        std::cout<<"caloinfo.forwardEndcap: "<<caloinfo.n_forwardEndcap<<std::endl;
                        std::cout<<"caloinfo.backwardEndcap: "<<caloinfo.n_backwardEndcap<<std::endl;
                    }

                    // Print out if theta angle is less than 20 degrees
                    if( theta_angle_barrel*180./M_PI < 20. ) {
                        std::cout<<"Barrel incident photon energy: "<<energy_barrel<<" (GeV)"<<std::endl;
                        std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                        std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                        std::cout<<"Barrel Theta angle from flight direction: "<<theta_angle_barrel<<" rad ("<<theta_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Barrel Phi angle between flight direction and ECAL normal: "<<phi_angle_barrel<<" rad ("<<phi_angle_barrel*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                        std::cout<<"Nhits in ECAL: "<<nhits_barrel<<std::endl;
                        std::cout<<"caloinfo.barrel: "<<caloinfo.n_barrel<<std::endl;
                        std::cout<<"caloinfo.forwardEndcap: "<<caloinfo.n_forwardEndcap<<std::endl;
                        std::cout<<"caloinfo.backwardEndcap: "<<caloinfo.n_backwardEndcap<<std::endl;
                    }
                }
}


void PeterProcessor::EndcapFunc( Vector3D mom, bool Sim_flag, bool Overlay_flag, std::vector<MCParticle*> Parents, const auto& mc, const auto& caloinfo ){
                
                double phi_angle_endcap = mom.phi();
                //std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;            
                double theta_angle_endcap = mom.theta();
                
                if(theta_angle_endcap > M_PI/2){ theta_angle_endcap = M_PI/2. -( M_PI - theta_angle_endcap); }
                else{
                    theta_angle_endcap = M_PI/2. - theta_angle_endcap;
                }

                int nhits_endcap = caloinfo.nhits;

                if ( theta_angle_endcap*180./M_PI <= 10.0 && (Overlay_flag == 1 || Sim_flag == 1) ){
                    std::cout<<"Nhits in ECAL: "<<nhits_endcap<<std::endl;
                    return;
                }

                //std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                double energy_endcap = mc->getEnergy();
                //std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                //std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                //std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;

                int Barrel_flag = 0; // not in barrel
                int Endcap_flag = 1; // In endcap

                int parent_size = static_cast<int>(Parents.size());
                if(parent_size > 1) {std::cout<<"More than one parent!"<<std::endl;}

                int Parent_type = Parents[0]->getPDG();
                //std::cout<<"Parent PDG: "<<Parent_type<<std::endl;

                // Write to file
                writeOutput(Barrel_flag, Endcap_flag, Sim_flag, Overlay_flag, phi_angle_endcap, theta_angle_endcap, energy_endcap, Parent_type, _AngSpecfile);

                // Print out if high energy Photon....
                if(energy_endcap > 100.){
                    std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                    std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                    std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                    std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                    std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                    std::cout<<"Nhits in ECAL: "<<nhits_endcap<<std::endl;
                    std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                    std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                    std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                }


                std::vector<MCParticle*> Parents_parents = Parents[0] -> getParents();


                int parent_parent_size = static_cast<int>(Parents_parents.size());
                if(parent_parent_size > 1) {std::cout<<"More than one parent of parent!"<<std::endl;}

                // skip the case of empty parent list
                if (Parents_parents.empty() == false){   

                    int Parent_parent_type = Parents_parents[0]->getPDG();


                    // Print out if parent was a pi0, which itself came from a tau...
                    if( abs(Parent_type) == 111 &&  abs(Parent_parent_type) == 15 ){
                        std::cout<<"Endcap incident photon energy: "<<energy_endcap<<" (GeV)"<<std::endl;
                        std::cout<<"Created in Simulation?: " << Sim_flag << std::endl;
                        std::cout<<"Is Overlay?: " << Overlay_flag << std::endl;
                        std::cout<<"Endcap Theta angle: "<<theta_angle_endcap<<" rad ("<<theta_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Endcap Phi angle: "<<phi_angle_endcap<<" rad ("<<phi_angle_endcap*180./M_PI<<" deg)"<<std::endl;
                        std::cout<<"Parent PDG: "<<Parent_type<<std::endl;
                        std::cout<<"Parent's parent PDG: "<<Parent_parent_type<<std::endl;
                        std::cout<<"Nhits in ECAL: "<<nhits_endcap<<std::endl;
                        std::cout<<"caloinfo.barrel"<<caloinfo.n_barrel<<std::endl;
                        std::cout<<"caloinfo.forwardEndcap"<<caloinfo.n_forwardEndcap<<std::endl;
                        std::cout<<"caloinfo.backwardEndcap"<<caloinfo.n_backwardEndcap<<std::endl;
                    }
                }
}

dd4hep::rec::Vector3D PeterProcessor::getBarrelNorm(double phi){
    /*Return the normal Vector3D of the octant ECAL plane selected by the phi angle.*/
    double rEcal = 1804.8; //Radius to the ECAL surface in mm
    int nSides = 8;
    double step = M_PI/nSides;

    // phi is in range[-pi, pi]. Check the side with singularity point.
    if( phi < (- M_PI + step) || phi > (M_PI - step) )
        return Vector3D(rEcal, M_PI, M_PI/2., Vector3D::spherical).unit();

    double ecal_phi = -M_PI + 2*step;
    for (int i=0; i < nSides-1; ++i){
        if (ecal_phi-step <= phi && phi < ecal_phi+step)
            return Vector3D(rEcal, ecal_phi, M_PI/2., Vector3D::spherical).unit();
        else ecal_phi += 2*step;
    }
    // this will never happen
    return Vector3D();
}


/*

bool PeterProcessor::hasEndcapHits(std::vector<EVENT::Cluster*> clusters){
    for(auto* cluster : clusters){
        for(auto* hit : cluster->getCalorimeterHits() ){
            CHT cht( hit->getType() );
            if ( cht.layout() == CHT::endcap ){
                std::cout<<"Hit "<<hit<<" is in the endcap here:"<<std::endl;
                std::cout<<cht<<std::endl;
                return true;
            }
        }
    }
    return false;
}

bool PeterProcessor::hasBarrelHits(std::vector<EVENT::Cluster*> clusters){
    for(auto* cluster : clusters){
        for(auto*hit : cluster->getCalorimeterHits() ){
            CHT cht( hit->getType() );
            if (cht.layout() == CHT::barrel){
                std::cout<<"Hit "<<hit<<" is in the barrel here:"<<std::endl;
                std::cout<<cht<<std::endl;
                return true;
            }            
        }
    }
    return false;
}

*/


void PeterProcessor::end(){
    _AngSpecfile.close();
    _invMassfile.close();
}
