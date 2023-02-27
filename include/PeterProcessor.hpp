#ifndef PeterProcessor_h
#define PeterProcessor_h 1

#include "EventDisplayer.h"
#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "UTIL/LCRelationNavigator.h"
#include "DDRec/Vector3D.h"
#include "DD4hep/Detector.h"
#include <vector>

class PeterProcessor : public marlin::Processor, EventDisplayer {
    friend class EventDisplayer;
    public:
        marlin::Processor* newProcessor() { return new PeterProcessor; }
        PeterProcessor();
        void init();
        void processEvent(EVENT::LCEvent* event);
        void end();
        void writeOutput(int Barrel_flag, bool Sim_flag, bool Overlay_flag, double Phi_Angle_Norm, double Theta_Angle_Norm, double Energy, double Parent_type, std::ofstream& AngSpecfile);

        EVENT::MCParticle* getLinkedMCParticle(EVENT::ReconstructedParticle* pfo, UTIL::LCRelationNavigator nav);
        dd4hep::rec::Vector3D getBarrelNorm(double phi);
        bool hasEndcapHits(std::vector<EVENT::Cluster*> clusters);
        bool hasBarrelHits(std::vector<EVENT::Cluster*> clusters); 
    private:
        std::ofstream _AngSpecfile;
        dd4hep::Detector& _detector = dd4hep::Detector::getInstance();
        int _nEvent{};
};

#endif
