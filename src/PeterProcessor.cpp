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

using namespace EVENT;
using dd4hep::rec::Vector3D;
using std::cout, std::endl, std::vector, std::string;
using UTIL::LCRelationNavigator;
PeterProcessor aPeterProcessor ;

PeterProcessor::PeterProcessor() : marlin::Processor("PeterProcessor"), EventDisplayer(this) {}


void PeterProcessor::init(){
    marlin::Global::EVENTSEEDER->registerProcessor(this);                                                
    DDMarlinCED::init(this);                                                                                                                                                                                                                                                                                                                                                               
    _AngSpecfile.open ("/nfs/dust/ilc/user/mckeownp/Angular_spectra/ILD_ECAL_Angular_Spectra_data/4jet_photons_angular_spectra.csv");
    _AngSpecfile << "Barrel (1)/ Endcap (0) Flag, Sim Flag, Overlay Flag, Phi angle to Norm, Theta angle to Norm, Energy (GeV), Parent PDG\n";

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
                                                                           

    LCCollection* pfos = event->getCollection("PandoraPFOs");
    LCRelationNavigator nav( event->getCollection("RecoMCTruthLink") );

    for (int i=0; i<pfos->getNumberOfElements(); ++i){
        std::cout<<endl<<"****************EVENT: "<<_nEvent<<"    PFO: "<<i+1<<"****************"<<std::endl;
        ReconstructedParticle* pfo = static_cast <ReconstructedParticle*> ( pfos->getElementAt(i) );
        MCParticle* mc = getLinkedMCParticle(pfo, nav);

        auto clusters = pfo->getClusters();
        if (mc->getPDG() != 22) continue;
        if (clusters.size() != 1) std::cout<<" BEWAAARE, SOMETHING IS SPPOOKY HERE. EVENT: "<<_nEvent<<"    PFO: "<<i+1<<std::endl;
        //if( hasEndcapHits(clusters) ) continue;

        Vector3D mom( mc->getMomentum() );

        bool Sim_flag = mc -> isCreatedInSimulation();
        bool Overlay_flag = mc -> isOverlay();

        std::vector<MCParticle*> Parents = mc -> getParents();  

        std::cout<<"***PHOTON MOMENTUM***"<<std::endl;
        std::cout<<mom<<std::endl;
        
        if( hasEndcapHits(clusters) == false &&  hasBarrelHits(clusters) == true ){
            //assuming photon doesn't change flight direction (it really shouldn't...)
            Vector3D ecalNorm = getBarrelNorm( mom.phi() );
            //std::cout<<"***ECAL NORM***"<<std::endl;
            //std::cout<<ecalNorm<<std::endl;
            double phi_angle_barrel = std::acos( mom.unit()*ecalNorm );
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

            /*
            AngSpecfile << Barrel_flag << ",";
            AngSpecfile << Sim_flag << ",";
            AngSpecfile << Overlay_flag << ",";
            AngSpecfile << phi_angle_barrel << ",";
            AngSpecfile << energy_barrel << ",";
            AngSpecfile << Parent_type << ",";
            */

        }
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

            /*
            AngSpecfile << Barrel_flag << ",";
            AngSpecfile << Sim_flag << ",";
            AngSpecfile << Overlay_flag << ",";
            AngSpecfile << theta_angle_endcap << ",";
            AngSpecfile << energy_endcap << ",";
            AngSpecfile << Parent_type << "\n";
            */

        }
        else continue;

        //drawDisplay(this, event, draw, mc, clusters);
    }
}

void PeterProcessor::writeOutput(int Barrel_flag, bool Sim_flag, bool Overlay_flag, double Phi_Angle_Norm, double Theta_Angle_Norm, double Energy, double Parent_type, std::ofstream& AngSpecfile){
    AngSpecfile << Barrel_flag << ",";
    AngSpecfile << Sim_flag << ",";
    AngSpecfile << Overlay_flag << ",";
    AngSpecfile << Phi_Angle_Norm << ",";
    AngSpecfile << Theta_Angle_Norm << ",";
    AngSpecfile << Energy << ",";
    AngSpecfile << Parent_type << "\n";
}

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


void PeterProcessor::end(){
    _AngSpecfile.close();
}