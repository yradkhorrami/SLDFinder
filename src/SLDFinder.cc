#include "SLDFinder.h"
using namespace lcio ;
using namespace marlin ;
SLDFinder aSLDFinder;

SLDFinder::SLDFinder() :

	Processor("SLDFinder"),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_nSLDecayOfBHadron(0),
	m_nSLDecayOfCHadron(0),
	m_nSLDecayOfTauLepton(0),
	m_nSLDecayTotal(0),
	m_nSLDecayToElectron(0),
	m_nSLDecayToMuon(0),
	m_nSLDecayToTau(0),
	m_pTFile(NULL),
	m_pTTree(NULL)
{
	registerInputCollection(	LCIO::MCPARTICLE,
					"MCParticleCollection" ,
					"Name of the MCParticle collection"  ,
					m_MCParticleCollection,
					std::string("MCParticle")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"MCTruthRecoLink" ,
					"Name of the MCParticle-ReconstructedParticle Relations collection"  ,
					m_MCTruthRecoLinkCollection ,
					std::string("MCTruthRecoLink")
				);

	registerInputCollection(	LCIO::LCRELATION,
					"RecoMCTruthLink" ,
					"Name of the ReconstructedParticle-MCParticle Relations collection"  ,
					m_RecoMCTruthLinkCollection ,
					std::string("RecoMCTruthLink")
				);

	registerOutputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"CorrectedPfoCollection",
					"Name of output pfo collection",
					m_SLDLeptonCollection,
					std::string("CorrectedPfoCollection")
				);

	registerProcessorParameter(	"fillRootTree",
					"Fill root tree to check processor performance",
					m_fillRootTree,
					bool(true)
				);

	registerProcessorParameter(	"RootFile",
	                                "Name of the output root file",
					m_rootFile,
					std::string("SLDFinder.root")
				);

}

void SLDFinder::init()
{
	streamlog_out(DEBUG) << "   init called  " << std::endl;
	printParameters();
	m_nRun = 0 ;
	m_nEvt = 0 ;
	m_nRunSum = 0;
	m_nEvtSum = 0;
	this->Clear();
	m_pTFile = new TFile(m_rootFile.c_str(), "recreate");

	m_pTTree = new TTree("SLDecays", "SLDecays");
	m_pTTree->SetDirectory(m_pTFile);
	m_pTTree->Branch("run", &m_nRun, "run/I");
	m_pTTree->Branch("event", &m_nEvt, "event/I");
	m_pTTree->Branch("nSLDecayOfBHadron", &m_nSLDecayOfBHadron, "nSLDecayOfBHadron/I");
	m_pTTree->Branch("nSLDecayOfCHadron", &m_nSLDecayOfCHadron, "nSLDecayOfCHadron/I");
	m_pTTree->Branch("nSLDecayOfTauLepton", &m_nSLDecayOfTauLepton, "nSLDecayOfTauLepton/I");
	m_pTTree->Branch("nSLDecayTotal", &m_nSLDecayTotal, "nSLDecayTotal/I");
	m_pTTree->Branch("nSLDecayToElectron", &m_nSLDecayToElectron, "nSLDecayToElectron/I");
	m_pTTree->Branch("nSLDecayToMuon", &m_nSLDecayToMuon, "nSLDecayToMuon/I");
	m_pTTree->Branch("nSLDecayToTau", &m_nSLDecayToTau, "nSLDecayToTau/I");

}

void SLDFinder::processRunHeader()
{
	m_nRun = 0;
	m_nEvt = 0;
	++m_nRunSum;
}

void SLDFinder::Clear()
{
	m_nRun = 0;
	m_nEvt = 0;
	m_nSLDecayOfBHadron = 0;
	m_nSLDecayOfCHadron = 0;
	m_nSLDecayOfTauLepton = 0;
	m_nSLDecayTotal = 0;
	m_nSLDecayToElectron = 0;
	m_nSLDecayToMuon = 0;
	m_nSLDecayToTau = 0;
}

void SLDFinder::processEvent( EVENT::LCEvent *pLCEvent )
{
	m_nRun = pLCEvent->getRunNumber();
	m_nEvt = pLCEvent->getEventNumber();
	++m_nEvtSum;
	bool isSLD = false;
	int parentFlavour = 0;
	int leptonFlavour = 0;

	streamlog_out(MESSAGE) << "" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////	Processing event 	" << m_nEvt << "	////////////////////" << std::endl;
	streamlog_out(MESSAGE) << "	////////////////////////////////////////////////////////////////////////////" << std::endl;

	LCCollection *MCParticleCollection{};
	try
	{
		MCParticleCollection = pLCEvent->getCollection( m_MCParticleCollection );
		int nMCP = MCParticleCollection->getNumberOfElements();
		for ( int i_mcp = 0 ; i_mcp < nMCP ; ++i_mcp )
		{
			MCParticle *testMCP = dynamic_cast<EVENT::MCParticle*>( MCParticleCollection->getElementAt( i_mcp ) );
			if ( testMCP->getParents().size() == 0 ) continue;
			isSLD = false;
			parentFlavour = 0;
			leptonFlavour = 0;
			findSLDecay( testMCP , isSLD , parentFlavour , leptonFlavour );
			if ( isSLD )
			{
				if ( abs( parentFlavour ) == 4 ) ++m_nSLDecayOfCHadron;
				if ( abs( parentFlavour ) == 5 ) ++m_nSLDecayOfBHadron;
				if ( abs( parentFlavour ) == 15 ) ++m_nSLDecayOfTauLepton;
				if ( abs( leptonFlavour ) == 11 ) ++m_nSLDecayToElectron;
				if ( abs( leptonFlavour ) == 13 ) ++m_nSLDecayToMuon;
				if ( abs( leptonFlavour ) == 15 ) ++m_nSLDecayToTau;
				streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Found one semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG3) << "		Flavour of Parent Particle:	" << parentFlavour << std::endl;
				streamlog_out(DEBUG3) << "		Flavour of Daughter Particle:	" << leptonFlavour << std::endl;
				++m_nSLDecayTotal;
			}
		}
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
	}
	m_pTTree->Fill();
}

void SLDFinder::findSLDecay( EVENT::MCParticle *testMCP , bool &isSLD , int &parentFlavour , int &leptonFlavour )
{
	if ( ( abs( testMCP->getPDG() ) == 11 || abs( testMCP->getPDG() ) == 13 || abs( testMCP->getPDG() ) == 15 ) && testMCP->getGeneratorStatus() == 1 && !( testMCP->isOverlay() ) )
	{
		int leptonCharge = ( int ) testMCP->getCharge();
		for ( unsigned int i_parent = 0 ; i_parent < testMCP->getParents().size() ; ++i_parent )
		{
			MCParticle *parentMCP = testMCP->getParents()[ i_parent ];
			for ( unsigned int i_daughter = 0 ; i_daughter < parentMCP->getDaughters().size() ; ++i_daughter )
			{
				int expectedNeutrinoPDG = -1 * ( testMCP->getPDG() - leptonCharge );
				MCParticle *daughterMCP = parentMCP->getDaughters()[ i_daughter ];
				if ( daughterMCP->getPDG() == expectedNeutrinoPDG )
				{
					isSLD = true;
					if ( parentMCP->getGeneratorStatus() == 2 && ( floor( abs( parentMCP->getPDG() ) / 100 ) == 5 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 5 ) ) ) parentFlavour = 5;
					if ( parentMCP->getGeneratorStatus() == 2 && ( floor( abs( parentMCP->getPDG() ) / 100 ) == 4 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 4 ) ) ) parentFlavour = 4;
					if ( parentMCP->getGeneratorStatus() == 2 && abs( parentMCP->getPDG() ) == 15 ) parentFlavour = 15;
					if ( abs( testMCP->getPDG() ) == 11 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 11;
					if ( abs( testMCP->getPDG() ) == 13 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 13;
					if ( abs( testMCP->getPDG() ) == 15 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 15;
				}
			}
		}
	}
}

void SLDFinder::check( EVENT::LCEvent *pLCEvent )
{
	LCCollection *SLDLeptonCollection{};
	try
	{
		SLDLeptonCollection = pLCEvent->getCollection( m_SLDLeptonCollection );
		int nSLDs = SLDLeptonCollection->getNumberOfElements();
		streamlog_out(DEBUG) << " CHECK : processed events: " << m_nEvtSum << " (Number of found semi-leptonic decays: " << m_nSLDecayTotal << " , Number of extracted recoLeptons in output collection: " << nSLDs <<" )" << std::endl;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Out collection not found in event " << m_nEvt << std::endl;
	}
}

void SLDFinder::end()
{

	m_pTFile->cd();
	m_pTTree->Write();
	m_pTFile->Close();
	delete m_pTFile;
	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
