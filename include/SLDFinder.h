#ifndef SLDFinder_h
#define SLDFinder_h 1
#include <iostream>
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"
#include "EVENT/LCStrVec.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCParameters.h>

#include "TFile.h"
#include "TH1F.h"
#include "TH1F.h"
#include "TTree.h"

using namespace lcio ;
using namespace marlin ;

class SLDFinder : public Processor
{
	public:

		virtual Processor*  newProcessor()
		{
			return new SLDFinder;
		}
		SLDFinder();
		virtual ~SLDFinder() = default;
		SLDFinder(const SLDFinder&) = delete;
		SLDFinder& operator=(const SLDFinder&) = delete;
		virtual void init();
		virtual void processRunHeader();
		virtual void processEvent( EVENT::LCEvent *pLCEvent );
		virtual void findSLDecay( EVENT::MCParticle *testMCP , bool &isSLD , int &parentFlavour , int &leptonFlavour );
		virtual void check( EVENT::LCEvent *pLCEvent );
		virtual void end();
		void Clear();

	private:

		typedef std::vector<int>		IntVector;
		typedef std::vector<double>		DoubleVector;
		typedef std::vector<float>		FloatVector;

		std::string				m_MCParticleCollection{};
		std::string				m_MCTruthRecoLinkCollection{};
		std::string				m_RecoMCTruthLinkCollection{};
		std::string				m_SLDLeptonCollection{};
		std::string				m_rootFile{};
		bool					m_fillRootTree = true;

		int					n_CSLD;
		int					n_BSLD;

		int					m_nRun;
		int					m_nEvt;
		int					m_nRunSum;
		int					m_nEvtSum;
		int					m_nSLDecayOfBHadron;
		int					m_nSLDecayOfCHadron;
		int					m_nSLDecayTotal;
		int					m_nSLDecayToElectron;
		int					m_nSLDecayToMuon;
		int					m_nSLDecayToTau;
		IntVector				m_SLDType{};
		IntVector				m_SLDMode{};
		TFile					*m_pTFile{};
		TTree					*m_pTTree{};
		TH1F					*h_CSLD{};
		TH1F					*h_BSLD{};

};

#endif
