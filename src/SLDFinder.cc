#include "SLDFinder.h"
using namespace lcio ;
using namespace marlin ;
SLDFinder aSLDFinder;

SLDFinder::SLDFinder() :

	Processor("SLDFinder"),
	n_CSLD(0),
	n_BSLD(0),
	m_nRun(0),
	m_nEvt(0),
	m_nRunSum(0),
	m_nEvtSum(0),
	m_nSLDecayOfBHadron(0),
	m_nSLDecayOfCHadron(0),
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

	registerInputCollection(	LCIO::RECONSTRUCTEDPARTICLE,
					"JetCollection",
					"Name of input jet collection",
					m_inputJetCollection,
					std::string("Durham_nJets")
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
					"SLDLeptons",
					"Name of leptons from semi-leptonic decays collection",
					m_SLDLeptonCollection,
					std::string("SLDLeptons")
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
	m_pTTree->Branch("nSLDecayTotal", &m_nSLDecayTotal, "nSLDecayTotal/I");
	m_pTTree->Branch("nSLDecayPerJet", &m_nSLDecayPerJet);
	m_pTTree->Branch("nSLDecayToElectron", &m_nSLDecayToElectron, "nSLDecayToElectron/I");
	m_pTTree->Branch("nSLDecayToMuon", &m_nSLDecayToMuon, "nSLDecayToMuon/I");
	m_pTTree->Branch("nSLDecayToTau", &m_nSLDecayToTau, "nSLDecayToTau/I");
	m_pTTree->Branch("SLDType", &m_SLDType);
	m_pTTree->Branch("SLDMode", &m_SLDMode);
	h_CSLD = new TH1F( "SLDOfCHadron" , ";" , 3 , 0.0 , 3.0 ); n_CSLD = 0;
	h_CSLD->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_CSLD->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_CSLD->GetXaxis()->SetBinLabel(3,"#tau^{#pm}");
	h_BSLD = new TH1F( "SLDOfBHadron" , ";" , 3 , 0.0 , 3.0 ); n_BSLD = 0;
	h_BSLD->GetXaxis()->SetBinLabel(1,"e^{#pm}");
	h_BSLD->GetXaxis()->SetBinLabel(2,"#mu^{#pm}");
	h_BSLD->GetXaxis()->SetBinLabel(3,"#tau^{#pm}");
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
	m_nSLDecayTotal = 0;
	m_nSLDecayPerJet.clear();
	m_nSLDecayToElectron = 0;
	m_nSLDecayToMuon = 0;
	m_nSLDecayToTau = 0;
	m_SLDType.clear();
	m_SLDMode.clear();
}

void SLDFinder::processEvent( EVENT::LCEvent *pLCEvent )
{
	this->Clear();
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
	LCCollection *jetCollection{};
	try
	{
		MCParticleCollection = pLCEvent->getCollection( m_MCParticleCollection );
		jetCollection = pLCEvent->getCollection( m_inputJetCollection );
		LCRelationNavigator RecoMCParticleNav( pLCEvent->getCollection( m_RecoMCTruthLinkCollection ) );
		LCRelationNavigator MCParticleRecoNav( pLCEvent->getCollection( m_MCTruthRecoLinkCollection ) );
		int nMCP = MCParticleCollection->getNumberOfElements();
		int nJets = jetCollection->getNumberOfElements();
		float weightPFOtoMCP = 0.0;
		float weightMCPtoPFO = 0.0;
		pfoVector jetVector{};
		bool recoLeptonIsInJet = false;
		int assignedJetIndex = -1;
		std::vector<int> nSLDinJet{};
		for ( int i_jet = 0 ; i_jet < nJets ; ++i_jet )
		{
			ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt( i_jet ) );
			jetVector.push_back( jet );
			nSLDinJet.push_back( 0 );
		}
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
				if ( abs( parentFlavour ) == 4 )
				{
					++m_nSLDecayOfCHadron;
					m_SLDType.push_back( 4 );
					++n_CSLD;
					if ( abs( leptonFlavour ) == 11 ) h_CSLD->Fill( 0.5 );
					if ( abs( leptonFlavour ) == 13 ) h_CSLD->Fill( 1.5 );
					if ( abs( leptonFlavour ) == 15 ) h_CSLD->Fill( 2.5 );
				}
				else if ( abs( parentFlavour ) == 5 )
				{
					++m_nSLDecayOfBHadron;
					m_SLDType.push_back( 5 );
					++n_BSLD;
					if ( abs( leptonFlavour ) == 11 ) h_BSLD->Fill( 0.5 );
					if ( abs( leptonFlavour ) == 13 ) h_BSLD->Fill( 1.5 );
					if ( abs( leptonFlavour ) == 15 ) h_BSLD->Fill( 2.5 );
				}
				if ( abs( leptonFlavour ) == 11 )
				{
					++m_nSLDecayToElectron;
					m_SLDMode.push_back( 11 );
				}
				else if ( abs( leptonFlavour ) == 13 )
				{
					++m_nSLDecayToMuon;
					m_SLDMode.push_back( 13 );
				}
				else if ( abs( leptonFlavour ) == 15 )
				{
					++m_nSLDecayToTau;
					m_SLDMode.push_back( 13 );
				}
				streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG3) << "	<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Found one semi-leptonic decay >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG3) << "		Flavour of Parent Particle:	" << parentFlavour << std::endl;
				streamlog_out(DEBUG3) << "		Flavour of Daughter Particle:	" << leptonFlavour << std::endl;
				if ( abs( parentFlavour ) == 4 || abs( parentFlavour ) == 5 || abs( parentFlavour ) == 15 ) ++m_nSLDecayTotal;
				weightPFOtoMCP = 0.0;
				weightMCPtoPFO = 0.0;
				ReconstructedParticle* linkedRecoLepton = getLinkedPFO( testMCP , RecoMCParticleNav , MCParticleRecoNav , true , false , weightPFOtoMCP , weightMCPtoPFO );
				assignedJetIndex = -1;
				recoLeptonIsInJet = false;
				getJetAssignedToParticle( linkedRecoLepton , jetVector , recoLeptonIsInJet , assignedJetIndex );
				if ( recoLeptonIsInJet ) ++nSLDinJet[ assignedJetIndex ];
				for ( int i_jet = 0 ; i_jet < nJets ; ++i_jet ) m_nSLDecayPerJet.push_back( nSLDinJet[ i_jet ] );
			}
		}
		streamlog_out(DEBUG3) << "		Number of semi-leptonic decays of B-hadron:	" << m_nSLDecayOfBHadron << std::endl;
		streamlog_out(DEBUG3) << "		Number of semi-leptonic decays of C-hadron:	" << m_nSLDecayOfCHadron << std::endl;
		streamlog_out(DEBUG3) << "		Number of semi-leptonic decays to electron:	" << m_nSLDecayToElectron << std::endl;
		streamlog_out(DEBUG3) << "		Number of semi-leptonic decays to muon:	" << m_nSLDecayToMuon << std::endl;
		streamlog_out(DEBUG3) << "		Number of semi-leptonic decays to tau-lepton:	" << m_nSLDecayToTau << std::endl;
	}
	catch(DataNotAvailableException &e)
	{
		streamlog_out(MESSAGE) << "	Input collection not found in event " << m_nEvt << std::endl;
	}
	m_pTTree->Fill();
}

void SLDFinder::findSLDecay( EVENT::MCParticle *testMCP , bool &isSLD , int &parentFlavour , int &leptonFlavour )
{
	if ( ( abs( testMCP->getPDG() ) == 11 || abs( testMCP->getPDG() ) == 13 || abs( testMCP->getPDG() ) == 15 ) && !( testMCP->isOverlay() ) )// && testMCP->getGeneratorStatus() == 1 )
	{
		int leptonCharge = ( int ) testMCP->getCharge();
		for ( unsigned int i_parent = 0 ; i_parent < testMCP->getParents().size() ; ++i_parent )
		{
			MCParticle *parentMCP = testMCP->getParents()[ i_parent ];
			if ( parentMCP->getGeneratorStatus() == 2 && ( floor( abs( parentMCP->getPDG() ) / 100 ) == 4 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 4 ) || floor( abs( parentMCP->getPDG() ) / 100 ) == 5 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 5 ) ) )
			{
				for ( unsigned int i_daughter = 0 ; i_daughter < parentMCP->getDaughters().size() ; ++i_daughter )
				{
					int expectedNeutrinoPDG = -1 * ( testMCP->getPDG() - leptonCharge );
					MCParticle *daughterMCP = parentMCP->getDaughters()[ i_daughter ];
					if ( daughterMCP->getPDG() == expectedNeutrinoPDG )
					{
						isSLD = true;
						if ( floor( abs( parentMCP->getPDG() ) / 100 ) == 4 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 4 ) ) parentFlavour = 4;
						if ( floor( abs( parentMCP->getPDG() ) / 100 ) == 5 || ( floor( abs( parentMCP->getPDG() ) / 1000 ) == 5 ) ) parentFlavour = 5;
						if ( abs( testMCP->getPDG() ) == 11 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 11;
						if ( abs( testMCP->getPDG() ) == 13 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 13;
//						if ( abs( testMCP->getPDG() ) == 15 && testMCP->getGeneratorStatus() == 1 ) leptonFlavour = 15;
						if ( abs( testMCP->getPDG() ) == 15 ) leptonFlavour = 15;
					}
				}
			}
		}
	}
}

EVENT::ReconstructedParticle* SLDFinder::getLinkedPFO( EVENT::MCParticle *mcParticle , LCRelationNavigator RecoMCParticleNav , LCRelationNavigator MCParticleRecoNav , bool getChargedPFO , bool getNeutralPFO , float &weightPFOtoMCP , float &weightMCPtoPFO )
{
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "	Look for PFO linked to visible MCParticle:" << std::endl;
	streamlog_out(DEBUG1) << *mcParticle << std::endl;
	ReconstructedParticle* linkedPFO{};
	bool foundlinkedPFO = false;
	const EVENT::LCObjectVec& PFOvec = MCParticleRecoNav.getRelatedToObjects( mcParticle );
	const EVENT::FloatVec&  PFOweightvec = MCParticleRecoNav.getRelatedToWeights( mcParticle );
	streamlog_out(DEBUG0) << "	Visible MCParticle is linked to " << PFOvec.size() << " PFO(s)" << std::endl;
	weightPFOtoMCP = 0.0;
	weightMCPtoPFO = 0.0;
	double maxweightPFOtoMCP = 0.;
	double maxweightMCPtoPFO = 0.;
	int iPFOtoMCPmax = -1;
	int iMCPtoPFOmax = -1;
	for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
	{
		double pfo_weight = 0.0;
		double trackWeight = ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0;
		double clusterWeight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
		if ( getChargedPFO && !getNeutralPFO )
		{
			pfo_weight = trackWeight;
		}
		else if ( getNeutralPFO && !getChargedPFO )
		{
			pfo_weight = clusterWeight;
		}
		else
		{
			pfo_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
		}
		streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << trackWeight << " , Cluster: " << clusterWeight << ")" << std::endl;
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
		if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
		{
			maxweightMCPtoPFO = pfo_weight;
			iMCPtoPFOmax = i_pfo;
			streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
		}
	}
	if ( getChargedPFO && maxweightMCPtoPFO < 0.8 )
	{
		streamlog_out(DEBUG1) << "	MCParticle has link weight lower than 0.8 ( " << maxweightMCPtoPFO << " ), looking for linked PFO in clusters" << std::endl;
		for ( unsigned int i_pfo = 0; i_pfo < PFOvec.size(); i_pfo++ )
		{
			double pfo_weight = ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0;
			streamlog_out(DEBUG0) << "	Visible MCParticle linkWeight to PFO: " << PFOweightvec.at( i_pfo ) << " (Track: " << ( int( PFOweightvec.at( i_pfo ) ) % 10000 ) / 1000.0 << " , Cluster: " << ( int( PFOweightvec.at( i_pfo ) ) / 10000 ) / 1000.0 << ")" << std::endl;
			ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( i_pfo );
			if ( pfo_weight > maxweightMCPtoPFO )//&& track_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightMCPtoPFO = pfo_weight;
				iMCPtoPFOmax = i_pfo;
				streamlog_out(DEBUG0) << "	PFO at index: " << testPFO->id() << " has TYPE: " << testPFO->getType() << " and MCParticle to PFO link weight is " << pfo_weight << std::endl;
			}
		}
	}
	if ( iMCPtoPFOmax != -1 )
	{
		ReconstructedParticle *testPFO = (ReconstructedParticle *) PFOvec.at( iMCPtoPFOmax );
		const EVENT::LCObjectVec& MCPvec = RecoMCParticleNav.getRelatedToObjects( testPFO );
		const EVENT::FloatVec&  MCPweightvec = RecoMCParticleNav.getRelatedToWeights( testPFO );
		for ( unsigned int i_mcp = 0; i_mcp < MCPvec.size(); i_mcp++ )
		{
			double mcp_weight = 0.0;
			double trackWeight = ( int( MCPweightvec.at( i_mcp ) ) % 10000 ) / 1000.0;
			double clusterWeight = ( int( MCPweightvec.at( i_mcp ) ) / 10000 ) / 1000.0;
			if ( getChargedPFO && !getNeutralPFO )
			{
				mcp_weight = trackWeight;
			}
			else if ( getNeutralPFO && !getChargedPFO )
			{
				mcp_weight = clusterWeight;
			}
			else
			{
				mcp_weight = ( trackWeight > clusterWeight ? trackWeight : clusterWeight );
			}
			MCParticle *testMCP = (MCParticle *) MCPvec.at( i_mcp );
			if ( mcp_weight > maxweightPFOtoMCP )//&& mcp_weight >= m_MinWeightTrackMCTruthLink )
			{
				maxweightPFOtoMCP = mcp_weight;
				iPFOtoMCPmax = i_mcp;
				streamlog_out(DEBUG0) << "	MCParticle at index: " << testMCP->id() << " has PDG: " << testMCP->getPDG() << " and PFO to MCParticle link weight is " << mcp_weight << std::endl;
			}
		}
		if ( iPFOtoMCPmax != -1 )
		{
			if ( MCPvec.at( iPFOtoMCPmax ) == mcParticle )
			{
				linkedPFO = testPFO;
				foundlinkedPFO = true;
			}
		}
	}

	if( foundlinkedPFO )
	{
		streamlog_out(DEBUG1) << "	Linked PFO to MCParticle found successfully " << std::endl;
		streamlog_out(DEBUG1) << *linkedPFO << std::endl;
		weightPFOtoMCP = maxweightPFOtoMCP;
		weightMCPtoPFO = maxweightMCPtoPFO;
		return linkedPFO;
	}
	else
	{
		streamlog_out(DEBUG1) << "	Couldn't Find a PFO linked to MCParticle" << std::endl;
		return NULL;
	}
}

void SLDFinder::getJetAssignedToParticle( EVENT::ReconstructedParticle *particle , std::vector<EVENT::ReconstructedParticle*> jetVector , bool &foundParticleInJet , int &assignedJetIndex )
{
	foundParticleInJet = false;
//	EVENT::ReconstructedParticle *assignedJet = NULL;
	streamlog_out(DEBUG1) << "" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << "------------- Looking for particle (" << particle->id() << ") in jets --------------" << std::endl;
	streamlog_out(DEBUG1) << "----------------------------------------------------------------------" << std::endl;
	streamlog_out(DEBUG1) << *particle << std::endl;
	for ( unsigned int i_jet = 0 ; i_jet < jetVector.size() ; ++i_jet )
	{
		EVENT::ReconstructedParticle* jet = jetVector[ i_jet ];
		streamlog_out(DEBUG0) << "	Looking for particle in jet " << i_jet << std::endl;
		streamlog_out(DEBUG0) << *jet << std::endl;
		int nParticles = ( jet->getParticles() ).size();
		for ( int i_particle = 0 ; i_particle < nParticles ; ++i_particle )
		{
			ReconstructedParticle* testParticle = jet->getParticles()[ i_particle ];
			streamlog_out(DEBUG0) << "	Checking particle " << i_particle << " (" << testParticle->id() << ") in jet" << std::endl;
			if ( testParticle == particle )
			{
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<  Found particle in jet  >>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
				streamlog_out(DEBUG2) << *testParticle << std::endl;
				streamlog_out(DEBUG2) << "" << std::endl;
//				assignedJet = jet;
				foundParticleInJet = true;
				assignedJetIndex = i_jet;
			}
		}
	}
//	return assignedJet;
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
	h_CSLD->Scale( 100.0 / n_CSLD ); h_CSLD->GetYaxis()->SetRangeUser(0.0, 100.0); h_CSLD->Write();
	h_BSLD->Scale( 100.0 / n_BSLD ); h_BSLD->GetYaxis()->SetRangeUser(0.0, 100.0); h_BSLD->Write();
	m_pTFile->Close();
	delete m_pTFile;
	std::cout << " END : processed events: " << m_nEvtSum << std::endl;

}
