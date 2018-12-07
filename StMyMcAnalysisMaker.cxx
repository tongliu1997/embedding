/*************************************************
 *
 * $Id: StMyMcAnalysisMaker.cxx,v 1.40 2015/03/13 18:44:50 perev Exp $
 * $Log: StMyMcAnalysisMaker.cxx,v $
 *
 *
 * Examples that use the structures of
 * StMcEvent and StAssociationMaker
 *
 *************************************************/
#include <assert.h>
#include <Stiostream.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TFile.h"

#include "StMyMcAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

#include "StMessMgr.h"

#include "StAssociationMaker/StAssociationMaker.h"
#include "StAssociationMaker/StTrackPairInfo.hh"

#include "StThreeVectorF.hh"

#include "StEventTypes.h"

#include "StMcEventTypes.hh"

#include "StMcEvent.hh"
#include "StMcEvent/StMcEvent.hh"
#include "StMcEvent/StMcVertex.hh"
#include "StMcEvent/StMcTrack.hh"
// Define data Members for the histograms
const Int_t   StMyMcAnalysisMaker::mNumDeltaX = 50;
const Int_t   StMyMcAnalysisMaker::mNumDeltaZ = 50;
const Float_t StMyMcAnalysisMaker::mMinDeltaX = -0.52;
const Float_t StMyMcAnalysisMaker::mMaxDeltaX =  0.52;
const Float_t StMyMcAnalysisMaker::mMinDeltaZ = -0.52;
const Float_t StMyMcAnalysisMaker::mMaxDeltaZ =  0.52;
//struct TpcHitMRPair_t {
//  Float_t sector, row, isDet,
//    xM, yM, zM, pxM, pyM, pzM, dEM, dSM, nM,
//    xR, yR, zR, dER, IdM, IdR, qR, nR;
//};
//static const Char_t *vTpcHitMRPair = "sector:row:isDet:xM:yM:zM:pxM:pyM:pzM:dEM:dSM:nM:xR:yR:zR:dER:IdM:IdR:qR:nR";
//static TpcHitMRPair_t TpcHitMRPair;
//struct SvtHitMRPair_t {
//  Float_t barrel, layer, ladder, wafer,  hybrid, index,
//    xM, yM, zM, pxM, pyM, pzM, dEM, dSM, nM,
//    xR, yR, zR, dER, IdM, IdR, qR, nR;
//};
//static const Char_t *vSvtHitMRPair = "barrel:layer:ladder:wafer:hybrid:index:xM:yM:zM:pxM:pyM:pzM:dEM:dSM:nM:xR:yR:zR:dER:IdM:IdR:qR:nR";
//static SvtHitMRPair_t SvtHitMRPair;
//struct SsdHitMRPair_t {
//  Float_t ladder, wafer,
//    xM, yM, zM, pxM, pyM, pzM, dEM, dSM, nM,
//    xR, yR, zR, dER, IdM, IdR, qR, nR;
//};
//static const Char_t *vSsdHitMRPair = "ladder:wafer:xM:yM:zM:pxM:pyM:pzM:dEM:dSM:nM:xR:yR:zR:dER:IdM:IdR:qR:nR";
//static SsdHitMRPair_t SsdHitMRPair;
ClassImp(StMyMcAnalysisMaker)

//_________________________________________________
StMyMcAnalysisMaker::StMyMcAnalysisMaker(const char *OutName, const char *name, const char *title):
  StMaker(name,title)//, mNtupleFile(0)
{
    mOutName=OutName;
    //  StMyMcAnalysisMaker Constructor
    // - zero all pointers defined in the header file
/*    mAssociationCanvas = 0;
    mMomResolution  = 0;
    mHitResolution  = 0;   
    mSvtHitResolution  = 0;   
    mSsdHitResolution  = 0;   
    coordRec        = 0;  
    coordMcPartner  = 0;
    mTrackNtuple    = 0;
    mTpcHitNtuple   = 0;
    mSvtHitNtuple   = 0;
    mSsdHitNtuple   = 0;
*/

}

//_________________________________________________
StMyMcAnalysisMaker::~StMyMcAnalysisMaker()
{
/**/
}

//_____________________________________________________________________________

void StMyMcAnalysisMaker::Clear(const char*)
{
    // StMyMcAnalysisMaker - Clear,
    // Don't delete the canvas, so that it stays if needed
    
    StMaker::Clear();
}

//_________________________________________________
Int_t StMyMcAnalysisMaker::Finish()
{
/*  if (mNtupleFile) {
    mNtupleFile->Write();
    mNtupleFile->Close();
  }
*/
  if(mOutName!=NULL){
    char mOutHistName[255];
    sprintf(mOutHistName,"%s_hists.root",mOutName);
//    TString HistName=mOutHistName;
    TFile *fout = new TFile(mOutHistName,"RECREATE");
    fout->cd();
    WriteHistograms();
    fout->Close(); 
  }
    return StMaker::Finish();
}


void StMyMcAnalysisMaker::WriteHistograms(){
mMcPt->Write();
mTrkEff->Write();
mRcPt->Write();
mMcNhits->Write();
mRcNhits->Write();
mMcPt_eta->Write();
mRcPt_eta->Write();
mMceta_phi->Write();
mRceta_phi->Write();
Match_Pt_eta->Write();
mResponseMtx->Write();
mHitResolution->Write();
mSvtHitResolution->Write();
mMomResolution->Write();
mSsdHitResolution->Write();
coordRec->Write();
coordMcPartner->Write();
}

//_________________________________________________
Int_t StMyMcAnalysisMaker::Init()
{
    // StMyMcAnalysisMaker - Init
    SetZones();  // This is my method to set the zones for the canvas.

//    mNtupleFile = GetTFile();
 //   if (mNtupleFile) {mNtupleFile->cd(); mNtupleFile = 0;}
//    else {
//    char mOutNtupleName[255];
//    sprintf(mOutNtupleName,"%s_TrackMapNtuple.root",mOutName);
//    mNtupleFile = new TFile(mOutNtupleName,"RECREATE","Track Ntuple");}
    // Book Histograms Here so they can be found and deleted by Victor's chain (I hope).

    mMcPt = new TH1F("MCPt","Monte Carlo Transverse Momentum",80,0,2);
    mMcPt->SetXTitle("Pt(GeV/c)");
    mMcPt->SetYTitle("Counts");
    mTrkEff = new TH1F("TrkEff","Monte Carlo Tracking Efficiency",80,0,2);
    mTrkEff->SetXTitle("Pt(GeV/c)");
//    mTrkEff->SetYTitle("Matched tracks");
    mRcPt = new TH1F("RCPt","Reconstructed Transverse Momentum",80,0,2);
    mRcPt->SetXTitle("Pt(GeV/c)");
    mRcPt->SetYTitle("Counts");
    mResponseMtx = new TH2F("ResMtx","Response Matrix",80,0,2,80,0,2);
    mResponseMtx->SetXTitle("Monte Carlo Pt(GeV/c)");
    mResponseMtx->SetYTitle("Reconstructed Pt(GeV/c)");

    mMcNhits = new TH1F("MCNhits","NHits of MC Tracks",60,-0.5,59.5);
    mMcNhits->GetXaxis()->SetTitle("NHits"); 
    mRcNhits = new TH1F("RCNhits","NHits of Reconstructed Tracks",60,-0.5,59.5);
    mRcNhits->GetXaxis()->SetTitle("NHits"); 
    mMcPt_eta = new TH2F("MCPtEta","Pt_eta of MC tracks",80,0,20,80,-1,1);
    mMcPt_eta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mMcPt_eta->GetYaxis()->SetTitle("#eta");
    mRcPt_eta = new TH2F("RCPtEta","Pt_eta of RC tracks",80,0,20,80,-1,1);
    mRcPt_eta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    mMcPt_eta->GetYaxis()->SetTitle("#eta");
    Match_Pt_eta = new TH2F("MatchPtEta","Pt_eta of Matched tracks",80,0,20,80,-1,1);
    Match_Pt_eta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Match_Pt_eta->GetYaxis()->SetTitle("#eta");
    mMceta_phi = new TH2F("MCetaphi","Eta_Phi distribution of MC tracks",80,-1,1,80,0,6.29);
    mRceta_phi = new TH2F("RCetaphi","Eta_phi distribution of RC tracks",80,-1,1,80,0,6.29);
    Match_eta_phi = new TH2F("Matchetaphi","Eta_phi matching distribution",80,-1,1,80,0,6.29);


    mHitResolution = new TH2F("hitRes","Delta Z Vs Delta X for Hits",
			     mNumDeltaX,mMinDeltaX,mMaxDeltaX,mNumDeltaZ,mMinDeltaZ,mMaxDeltaZ);
    mHitResolution->SetXTitle("Delta X (cm)");
    mHitResolution->SetYTitle("Delta Z (cm)");
    mSvtHitResolution = new TH2F("SvtHitRes","Delta Z Vs Delta X for SvtHits",
			     mNumDeltaX,mMinDeltaX,mMaxDeltaX,mNumDeltaZ,mMinDeltaZ,mMaxDeltaZ);
    mSvtHitResolution->SetXTitle("Delta X (cm)");
    mSvtHitResolution->SetYTitle("Delta Z (cm)");

    mSsdHitResolution = new TH2F("SsdHitRes","Delta Z Vs Delta X for SsdHits",
			     mNumDeltaX,mMinDeltaX,mMaxDeltaX,mNumDeltaZ,mMinDeltaZ,mMaxDeltaZ);
    mSsdHitResolution->SetXTitle("Delta X (cm)");
    mSsdHitResolution->SetYTitle("Delta Z (cm)");

    mMomResolution = new TH1F("momRes","(|p| - |pmc|)/|pmc|",100,-1.,1.);
    mMomResolution->SetXTitle("Resolution (%)");
    cout<<"Defining histograms"<<endl;
    coordRec = new TH2F("coordRc","X vs Y pos. of Hits", 100, -200, 200, 100, -200, 200);
    coordRec->SetXTitle("X (cm)");
    coordRec->SetYTitle("Y (cm)");
    
    coordMcPartner = new TH2F("coordMc","X vs Y pos. of Hits", 100, -200, 200, 100, -200, 200);
    coordMcPartner->SetXTitle("X (cm)");
    coordMcPartner->SetYTitle("Y (cm)");

    // Define the file for the Ntuple, otherwise it won't be available later.
    // one must define the file _after_ the histograms are booked, otherwise they are
    // not owned by the maker, but are stored in the file, breaking the code in StAssociator.

//    mNtupleFile = new TFile("TrackMapNtuple.root","RECREATE","Track Ntuple");
    
//    const char* vars = "px:py:pz:p:pxrec:pyrec:pzrec:prec:commTpcHits:hitDiffX:hitDiffY:hitDiffZ:mcTrkId:mostCommIdTruth:nHitsIdTruth:nMcHits:nFitPts:nDetPts:quality";
//    mTrackNtuple = new TNtuple("TrackNtuple","Track Pair Info",vars);
//    mTrackNtuple->SetAutoSave(100000000);
//    if (! m_Mode || m_Mode & 0x1) {
//      mTpcHitNtuple = new TNtuple("TpcHitNtuple","the TPC hit pairs Info",vTpcHitMRPair);
//      mTpcHitNtuple->SetAutoSave(100000000);
//    }
//    if (! m_Mode || m_Mode & 0x2) {
//      mSvtHitNtuple = new TNtuple("SvtHitNtuple","the SVT hit pairs Info",vSvtHitMRPair);
//      mSvtHitNtuple->SetAutoSave(100000000);
//    }
//    if (! m_Mode || m_Mode & 0x4) {
//      mSsdHitNtuple = new TNtuple("SsdHitNtuple","the SSD hit pairs Info",vSsdHitMRPair);
//      mSsdHitNtuple->SetAutoSave(100000000);
//    }
    return StMaker::Init();
}
//_________________________________________________
Int_t StMyMcAnalysisMaker::Make()
{
  // Get the pointers we need, we have to use the titles we gave them in the
  // macro.  I just used the defaults.
  

  StEvent* rEvent =  (StEvent*) GetInputDS("StEvent");
//  const StMcTrack  *mTrack = 0;
  
  // StMcEvent
  StMcEvent* mEvent = (StMcEvent*) GetDataSet("StMcEvent");
  
  // StAssociationMaker
  StAssociationMaker* assoc = 0;
  assoc = (StAssociationMaker*) GetMaker("StAssociationMaker");
  
  // the Multimaps...
  rcTpcHitMapType* theHitMap   = assoc->rcTpcHitMap();
  mcTpcHitMapType* theMcHitMap = assoc->mcTpcHitMap();
  rcSvtHitMapType* svtHitMap   = assoc->rcSvtHitMap();
  mcSvtHitMapType* svtMcHitMap = assoc->mcSvtHitMap();
  rcSsdHitMapType* ssdHitMap   = assoc->rcSsdHitMap();
  mcSsdHitMapType* ssdMcHitMap = assoc->mcSsdHitMap();
  rcTrackMapType*  theTrackMap = assoc->rcTrackMap();
  mMcTrackMap = assoc->mcTrackMap();
  mcV0MapType* theMcV0Map      = assoc->mcV0Map(); 
  
  if (!theHitMap) {
    gMessMgr->Warning() << "----------WARNING----------\n"
			<< "No Hit Map found for this event!" << endm;
    return kStWarn;
  }
  // Example: look at the position of the primary vertex
  //          Map is not needed for this, but it's a good check,
  //          tracking will not be good if primary vertex was not well placed.
  
  // First check whether the Primary Vertex is there at all.
  StThreeVectorD VertexPos(0,0,0);
  if (rEvent->primaryVertex()) {
    VertexPos = rEvent->primaryVertex()->position();
    int nprimtrk=rEvent->primaryVertex()->daughters().size();
    StTrack* rTrack;
    for(int itk=0;itk<nprimtrk;itk++){
      rTrack=rEvent->primaryVertex()->daughter(itk);
      int rcnhit=rTrack->topologyMap().numberOfHits(kTpcId);
      mRcNhits->Fill(rcnhit);
//      cout<<rTrack->pidTraits()<<endl;
    }
  }
  else {
    cout << "----------- WARNING ------------" << endl;
    cout << "No primary vertex from StEvent!!" << endl;
    cout << "Assume it is at (0,0,0)         " << endl;
    
  }

  float DeltaX;
  float DeltaZ;
  if (theHitMap){// && mTpcHitNtuple) {// TPC
    StTpcHitCollection* recHits = rEvent->tpcHitCollection();
    StMcTpcHitCollection* mcHits = mEvent->tpcHitCollection();
    assert (recHits || mcHits);
//    cout << "Making Hit Resolution Histogram..." << endl;
    // Loop over Rec Hits
    
    for (unsigned int iSector=0; iSector< recHits->numberOfSectors(); iSector++) {
      for (unsigned int iPadrow=0; iPadrow<recHits->sector(iSector)->numberOfPadrows();
	   iPadrow++) {
	for (StTpcHitIterator iter = recHits->sector(iSector)->padrow(iPadrow)->hits().begin();
	     iter != recHits->sector(iSector)->padrow(iPadrow)->hits().end();
	     iter++) {
	  const StTpcHit   *rhit = dynamic_cast<const StTpcHit   *> (*iter);
	  assert(rhit);
	  if (rhit->TestBit(StMcHit::kMatched)) 
	  {
	    pair<rcTpcHitMapIter,rcTpcHitMapIter>
	      recBounds = theHitMap->equal_range(rhit);
	    
	    for (rcTpcHitMapIter it2=recBounds.first; it2!=recBounds.second; ++it2){
	      
	      const StMcTpcHit *mhit = dynamic_cast<const StMcTpcHit *> ((*it2).second);
	      assert ( mhit);
	      DeltaX = rhit->position().x() - mhit->position().x();
	      DeltaZ = rhit->position().z() - mhit->position().z();
	      
	      mHitResolution->Fill(DeltaX,DeltaZ);
	      
/*	      memset (&TpcHitMRPair, 0, sizeof(TpcHitMRPair));
	      TpcHitMRPair.sector   = rhit->sector();
	      TpcHitMRPair.row      = rhit->padrow();
	      TpcHitMRPair.xR       = rhit->position().x();
	      TpcHitMRPair.yR       = rhit->position().y();
	      TpcHitMRPair.zR       = rhit->position().z();
	      TpcHitMRPair.dER      = rhit->charge();
	      TpcHitMRPair.IdR      = rhit->idTruth();
	      TpcHitMRPair.qR       = rhit->qaTruth();
	      TpcHitMRPair.nR       = theHitMap->count(rhit);
	      TpcHitMRPair.isDet    = mhit->isDet();
	      TpcHitMRPair.xM       = mhit->position().x();
	      TpcHitMRPair.yM       = mhit->position().y();
	      TpcHitMRPair.zM       = mhit->position().z();
	      TpcHitMRPair.pxM      = mhit->localMomentum().x();
	      TpcHitMRPair.pyM      = mhit->localMomentum().y();
	      TpcHitMRPair.pzM      = mhit->localMomentum().z();
	      TpcHitMRPair.dEM      = mhit->dE();
	      TpcHitMRPair.dSM      = mhit->dS();
*/
//	      mTrack     = mhit->parentTrack();
//	      if (mTrack) TpcHitMRPair.IdM = mTrack->key();
//	      else        TpcHitMRPair.IdM = 0;
//	      TpcHitMRPair.nM     = theMcHitMap->count(mhit);
//	      mTpcHitNtuple->Fill(&TpcHitMRPair.sector);
	    }
	  }
/*	  else {
	    memset (&TpcHitMRPair, 0, sizeof(TpcHitMRPair));
	    TpcHitMRPair.sector   = rhit->sector();
	    TpcHitMRPair.row      = rhit->padrow();
	    TpcHitMRPair.xR       = rhit->position().x();
	    TpcHitMRPair.yR       = rhit->position().y();
	    TpcHitMRPair.zR       = rhit->position().z();
	    TpcHitMRPair.dER      = rhit->charge();
	    TpcHitMRPair.IdR      = rhit->idTruth();
	    TpcHitMRPair.qR       = rhit->qaTruth();
	    TpcHitMRPair.nR       = theHitMap->count(rhit);
//	    mTpcHitNtuple->Fill(&TpcHitMRPair.sector);
	  }
*/	}
      }
    }
    for (unsigned int iSector=0;
	 iSector<mcHits->numberOfSectors(); iSector++) {
      
      if (Debug()) {cout << iSector + 1 << " "; flush(cout);}
      StMcTpcSectorHitCollection* tpcSectHitColl = mcHits->sector(iSector);
      for (unsigned int iPadrow=0;
	   iPadrow<tpcSectHitColl->numberOfPadrows();
	   iPadrow++) {
	StMcTpcPadrowHitCollection* tpcPadRowHitColl = tpcSectHitColl->padrow(iPadrow);
	for (unsigned int iHit=0;
	     iHit<tpcPadRowHitColl->hits().size();
	     iHit++){
	  const StMcTpcHit *mhit = dynamic_cast<const StMcTpcHit *> (tpcPadRowHitColl->hits()[iHit]);
	  assert (mhit);
	  if (mhit->TestBit(StMcHit::kMatched)) continue;
/*	  memset (&TpcHitMRPair, 0, sizeof(TpcHitMRPair));
	  TpcHitMRPair.sector   = mhit->sector();
	  TpcHitMRPair.row      = mhit->padrow();
	  TpcHitMRPair.isDet    = mhit->isDet();
	  TpcHitMRPair.xM       = mhit->position().x();
	  TpcHitMRPair.yM       = mhit->position().y();
	  TpcHitMRPair.zM       = mhit->position().z();
	  TpcHitMRPair.pxM      = mhit->localMomentum().x();
	  TpcHitMRPair.pyM      = mhit->localMomentum().y();
	  TpcHitMRPair.pzM      = mhit->localMomentum().z();
	  TpcHitMRPair.dEM      = mhit->dE();
	  TpcHitMRPair.dSM      = mhit->dS();
*/
//	  mTrack     = mhit->parentTrack();
//	  if (mTrack) TpcHitMRPair.IdM = mTrack->key();
//	  else        TpcHitMRPair.IdM = 0;
//	  mTpcHitNtuple->Fill(&TpcHitMRPair.sector);
	}
      }
    }
  }
  if (svtHitMap && svtMcHitMap){// && mSvtHitNtuple) {  // svt hits
    cout<<"There are hists to make"<<endl;
    StSvtHitCollection* recHits = rEvent->svtHitCollection();
    StMcSvtHitCollection* mcHits = mEvent->svtHitCollection();
    if (recHits && mcHits) {
      for (unsigned int iBarrel=0; iBarrel< recHits->numberOfBarrels(); iBarrel++) {
	for (unsigned int iLadder=0; iLadder<recHits->barrel(iBarrel)->numberOfLadders(); iLadder++) {
	  for (unsigned int iWafer = 0; iWafer < recHits->barrel(iBarrel)->ladder(iLadder)->numberOfWafers(); iWafer++) {
	    for (StSvtHitIterator iter = recHits->barrel(iBarrel)->ladder(iLadder)->wafer(iWafer)->hits().begin();
	       iter != recHits->barrel(iBarrel)->ladder(iLadder)->wafer(iWafer)->hits().end();
		 iter++) {
	      const StSvtHit   *rhit = dynamic_cast<const StSvtHit   *> (*iter);
	      assert(rhit);
	      if (rhit->TestBit(StMcHit::kMatched)) 
	      {
		pair<rcSvtHitMapIter,rcSvtHitMapIter> recBounds = svtHitMap->equal_range(rhit);
		for (rcSvtHitMapIter it2=recBounds.first; it2!=recBounds.second; ++it2){
		  const StMcSvtHit *mhit = dynamic_cast<const StMcSvtHit *> ((*it2).second);
		  assert ( mhit);
		  DeltaX = rhit->position().x() - mhit->position().x();
		  DeltaZ = rhit->position().z() - mhit->position().z();
		  mSvtHitResolution->Fill(DeltaX,DeltaZ);
/*		  memset (&SvtHitMRPair, 0, sizeof(SvtHitMRPair));
		  SvtHitMRPair.barrel   = rhit->barrel();
		  SvtHitMRPair.ladder   = rhit->ladder();
		  SvtHitMRPair.layer    = rhit->layer();
		  SvtHitMRPair.wafer    = rhit->wafer();
		  SvtHitMRPair.hybrid   = rhit->hybrid();
		  SvtHitMRPair.index    = rhit->index();
		  SvtHitMRPair.xR       = rhit->position().x();
		  SvtHitMRPair.yR       = rhit->position().y();
		  SvtHitMRPair.zR       = rhit->position().z();
		  SvtHitMRPair.dER      = rhit->charge();
		  SvtHitMRPair.IdR      = rhit->idTruth();
		  SvtHitMRPair.qR       = rhit->qaTruth();
		  SvtHitMRPair.nR       = svtHitMap->count(rhit);
		  SvtHitMRPair.xM       = mhit->position().x();
		  SvtHitMRPair.yM       = mhit->position().y();
		  SvtHitMRPair.zM       = mhit->position().z();
		  SvtHitMRPair.pxM      = mhit->localMomentum().x();
		  SvtHitMRPair.pyM      = mhit->localMomentum().y();
		  SvtHitMRPair.pzM      = mhit->localMomentum().z();
		  SvtHitMRPair.dEM      = mhit->dE();
		  SvtHitMRPair.dSM      = mhit->dS();
*/
//		  mTrack     = mhit->parentTrack();
//		  if (mTrack) SvtHitMRPair.IdM = mTrack->key();
//		  else        SvtHitMRPair.IdM = 0;
//		  SvtHitMRPair.nM     = svtMcHitMap->count(mhit);
//		  mSvtHitNtuple->Fill(&SvtHitMRPair.barrel);
		}
	      }
/*	      else {
		memset (&SvtHitMRPair, 0, sizeof(SvtHitMRPair));
		SvtHitMRPair.barrel   = rhit->barrel();
		SvtHitMRPair.ladder   = rhit->ladder();
		SvtHitMRPair.layer    = rhit->layer();
		SvtHitMRPair.wafer    = rhit->wafer();
		SvtHitMRPair.hybrid   = rhit->hybrid();
		SvtHitMRPair.index    = rhit->index();
		SvtHitMRPair.xR       = rhit->position().x();
		SvtHitMRPair.yR       = rhit->position().y();
		SvtHitMRPair.zR       = rhit->position().z();
		SvtHitMRPair.dER      = rhit->charge();
		SvtHitMRPair.IdR      = rhit->idTruth();
		SvtHitMRPair.qR       = rhit->qaTruth();
		SvtHitMRPair.nR       = svtHitMap->count(rhit);
		mSvtHitNtuple->Fill(&SvtHitMRPair.barrel);

	      }
*/
	    }
	  }
	}
      }
      for (unsigned int iBarrel=0; iBarrel< mcHits->numberOfBarrels(); iBarrel++) {
	for (unsigned int iLadder=0; iLadder<mcHits->barrel(iBarrel)->numberOfLadders(); iLadder++) {
	  for (unsigned int iWafer = 0; iWafer < mcHits->barrel(iBarrel)->ladder(iLadder)->numberOfWafers(); iWafer++) {
	    for (StMcSvtHitIterator iter = mcHits->barrel(iBarrel)->ladder(iLadder)->wafer(iWafer)->hits().begin();
	       iter != mcHits->barrel(iBarrel)->ladder(iLadder)->wafer(iWafer)->hits().end();
		 iter++) {
	      const StMcSvtHit   *mhit = dynamic_cast<const StMcSvtHit   *> (*iter);
	      assert (mhit);
	      if (mhit->TestBit(StMcHit::kMatched)) continue;
/*	      memset (&SvtHitMRPair, 0, sizeof(SvtHitMRPair));
	      SvtHitMRPair.barrel   = mhit->barrel();
	      SvtHitMRPair.ladder   = mhit->ladder();
	      SvtHitMRPair.layer    = mhit->layer();
	      SvtHitMRPair.wafer    = mhit->wafer();
	      SvtHitMRPair.hybrid   = mhit->hybrid();
	      SvtHitMRPair.index    = -1;
	      SvtHitMRPair.barrel   = mhit->barrel();
	      SvtHitMRPair.ladder   = mhit->ladder();
	      SvtHitMRPair.xM       = mhit->position().x();
	      SvtHitMRPair.yM       = mhit->position().y();
	      SvtHitMRPair.zM       = mhit->position().z();
	      SvtHitMRPair.pxM      = mhit->localMomentum().x();
	      SvtHitMRPair.pyM      = mhit->localMomentum().y();
	      SvtHitMRPair.pzM      = mhit->localMomentum().z();
	      SvtHitMRPair.dEM      = mhit->dE();
	      SvtHitMRPair.dSM      = mhit->dS();
*/
//	      mTrack     = mhit->parentTrack();
//	      if (mTrack) SvtHitMRPair.IdM = mTrack->key();
//	      else        SvtHitMRPair.IdM = 0;
//	      mSvtHitNtuple->Fill(&SvtHitMRPair.barrel);
	    }
	  }
	}
      }
    }
  }
  if (ssdHitMap && ssdMcHitMap){// && mSsdHitNtuple) {  // ssd hits
    StSsdHitCollection* recHits = rEvent->ssdHitCollection();
    StMcSsdHitCollection* mcHits = mEvent->ssdHitCollection();
    if (recHits && mcHits) {
      for (unsigned int iLadder=0; iLadder<recHits->numberOfLadders(); iLadder++) {
	for (unsigned int iWafer = 0; iWafer < recHits->ladder(iLadder)->numberOfWafers(); iWafer++) {
	  for (StSsdHitIterator iter = recHits->ladder(iLadder)->wafer(iWafer)->hits().begin();
	       iter != recHits->ladder(iLadder)->wafer(iWafer)->hits().end();
	       iter++) {
	    const StSsdHit   *rhit = dynamic_cast<const StSsdHit   *> (*iter);
	    assert(rhit);
	    if (rhit->TestBit(StMcHit::kMatched)) 
	      {
		pair<rcSsdHitMapIter,rcSsdHitMapIter> recBounds = ssdHitMap->equal_range(rhit);
		for (rcSsdHitMapIter it2=recBounds.first; it2!=recBounds.second; ++it2){
		  const StMcSsdHit *mhit = dynamic_cast<const StMcSsdHit *> ((*it2).second);
		  assert ( mhit);
		  DeltaX = rhit->position().x() - mhit->position().x();
		  DeltaZ = rhit->position().z() - mhit->position().z();
		  mSsdHitResolution->Fill(DeltaX,DeltaZ);
/*		  memset (&SsdHitMRPair, 0, sizeof(SsdHitMRPair));
		  SsdHitMRPair.ladder   = rhit->ladder();
		  SsdHitMRPair.wafer    = rhit->wafer();
		  SsdHitMRPair.xR       = rhit->position().x();
		  SsdHitMRPair.yR       = rhit->position().y();
		  SsdHitMRPair.zR       = rhit->position().z();
		  SsdHitMRPair.dER      = rhit->charge();
		  SsdHitMRPair.IdR      = rhit->idTruth();
		  SsdHitMRPair.qR       = rhit->qaTruth();
		  SsdHitMRPair.nR       = ssdHitMap->count(rhit);
		  SsdHitMRPair.xM       = mhit->position().x();
		  SsdHitMRPair.yM       = mhit->position().y();
		  SsdHitMRPair.zM       = mhit->position().z();
		  SsdHitMRPair.pxM      = mhit->localMomentum().x();
		  SsdHitMRPair.pyM      = mhit->localMomentum().y();
		  SsdHitMRPair.pzM      = mhit->localMomentum().z();
		  SsdHitMRPair.dEM      = mhit->dE();
		  SsdHitMRPair.dSM      = mhit->dS();
*/
//		  mTrack     = mhit->parentTrack();
//		  if (mTrack) SsdHitMRPair.IdM = mTrack->key();
//		  else        SsdHitMRPair.IdM = 0;
//		  SsdHitMRPair.nM     = ssdMcHitMap->count(mhit);
//		  mSsdHitNtuple->Fill(&SsdHitMRPair.ladder);
		}
	      }
/*	    else {
	      memset (&SsdHitMRPair, 0, sizeof(SsdHitMRPair));
	      SsdHitMRPair.ladder   = rhit->ladder();
	      SsdHitMRPair.wafer    = rhit->wafer();
	      SsdHitMRPair.xR       = rhit->position().x();
	      SsdHitMRPair.yR       = rhit->position().y();
	      SsdHitMRPair.zR       = rhit->position().z();
	      SsdHitMRPair.dER      = rhit->charge();
	      SsdHitMRPair.IdR      = rhit->idTruth();
	      SsdHitMRPair.qR       = rhit->qaTruth();
	      SsdHitMRPair.nR       = ssdHitMap->count(rhit);
//	      mSsdHitNtuple->Fill(&SsdHitMRPair.ladder);
	    }
*/
	  }
	}
      }
      for (unsigned int iLadder=0; iLadder<mcHits->numberOfLadders(); iLadder++) {
	for (unsigned int iWafer = 0; iWafer < mcHits->ladder(iLadder)->numberOfWafers(); iWafer++) {
	  for (StMcSsdHitIterator iter = mcHits->ladder(iLadder)->wafer(iWafer)->hits().begin();
	       iter != mcHits->ladder(iLadder)->wafer(iWafer)->hits().end();
	       iter++) {
	    const StMcSsdHit   *mhit = dynamic_cast<const StMcSsdHit   *> (*iter);
	    assert (mhit);
	    if (mhit->TestBit(StMcHit::kMatched)) continue;
/*	    memset (&SsdHitMRPair, 0, sizeof(SsdHitMRPair));
	    SsdHitMRPair.ladder   = mhit->ladder();
	    SsdHitMRPair.wafer    = mhit->wafer();
	    SsdHitMRPair.ladder   = mhit->ladder();
	    SsdHitMRPair.xM       = mhit->position().x();
	    SsdHitMRPair.yM       = mhit->position().y();
	    SsdHitMRPair.zM       = mhit->position().z();
	    SsdHitMRPair.pxM      = mhit->localMomentum().x();
	    SsdHitMRPair.pyM      = mhit->localMomentum().y();
	    SsdHitMRPair.pzM      = mhit->localMomentum().z();
	    SsdHitMRPair.dEM      = mhit->dE();
	    SsdHitMRPair.dSM      = mhit->dS();
*/
//	    mTrack     = mhit->parentTrack();
//	    if (mTrack) SsdHitMRPair.IdM = mTrack->key();
//	    else        SsdHitMRPair.IdM = 0;
//	    mSsdHitNtuple->Fill(&SsdHitMRPair.ladder);
	  }
	}
      }
    }
  }
//  cout << "Finished Making Histogram." << endl;
  
  if (!theTrackMap) {
    gMessMgr->Warning() << "----------WARNING----------\n"
			<< "No Track Map found for this event!" << endm;
    return kStWarn;
  }
  // Example: look at the magnitude of the momentum of
  //          the MC track associated with first track in track Collection
  
  StSPtrVecTrackNode& rcTrackNodes = rEvent->trackNodes();
  StTrackNode*        firstTrackNode = 0;
  const StGlobalTrack*      firstTrack = 0;
  if (! rcTrackNodes.size()) {
    cout << "Track Nodes List is empty!" << endl;
    return kStWarn;
  }
  firstTrackNode = *(rcTrackNodes.begin());
  if (! firstTrackNode) {
    cout << "No tracks in Track Nodes List!" << endl;
    return kStWarn;
  }
  firstTrack = dynamic_cast<const StGlobalTrack*>(firstTrackNode->track(global));
  if (! firstTrack) {
    cout << "GlobalTrack for first Track Node has not been found!" << endl;
    return kStWarn;
  }
//  cout << "MC Tracks associated with first Track in collection: " << theTrackMap->count(firstTrack) << endl;
  if (!theTrackMap->count(firstTrack) && theTrackMap->size()) {
//      cout << "First track in track node container was not associated.  Pick first track in map." << endl;
      firstTrack = theTrackMap->begin()->first; //first entry in the map is a pair, and first entry of pair is the global track
  }
  pair<rcTrackMapIter,rcTrackMapIter> trackBounds = theTrackMap->equal_range(firstTrack);
//  cout << "Momentum of First Track and of Associated Tracks (if there are any):" << endl;
  // Get the momentum of the track and compare it to MC Track.
  // Use primary track if available
  const StPrimaryTrack* pTrk = dynamic_cast<const StPrimaryTrack*>(firstTrack->node()->track(primary));
  StThreeVectorD recMom(0,0,0);
  if (pTrk)
    recMom = pTrk->geometry()->momentum();
  else
    recMom = firstTrack->geometry()->momentum();

  // Example: Make a histogram of the momentum resolution of the event
  //          Make an Ntuple with rec & monte carlo mom, mean hit difference, and # of common hits
  const StGlobalTrack* recTrack;
  const StPrimaryTrack* primTrk;
  const StMcTrack*     mcTrack;
  StThreeVectorD p(0,0,0);
  StThreeVectorD pmc(0,0,0);
  float diff =0;
  
  float* values = new float[19];
 
  StMcVertex* mcpVtx=mEvent->primaryVertex();
  if(mcpVtx){
    int ntracks=mcpVtx->numberOfDaughters();
    StMcTrack* mctk;
//    cout<<mcpVtx->position()<<endl;
    for(int itrk=0;itrk<ntracks;itrk++){
      mctk=mcpVtx->daughter(itrk);
      assert(mctk);
      if(!mctk) continue;
      mMcNhits->Fill(mctk->tpcHits().size());


      if(mctk->key()==0 && mctk->geantId()==0) {cout<<"啥？？";continue;} // not geant tracks 
//      if(mctk->IsPrimary()==true) cout<<"Primary ";
//      cout<<"MC track pt is "<<mctk->pt()<<" Geant Id is " <<mctk->geantId()<<endl;
//      StPrimaryTrack *rctk = 0;
      mMcPt->Fill(mctk->pt());
      mMcPt_eta->Fill(mctk->pt(),mctk->momentum().pseudoRapidity());
      mMceta_phi->Fill(mctk->momentum().pseudoRapidity(),mctk->momentum().phi()); 
//    cout<<mctk->tpcHits().size()<<"\t"<<mctk->pt()<<"\t";
      if(assoc) {
        if(mMcTrackMap) {

          StTrackPairInfo* matchedTrackPairInfo = findBestMatchedGlobal(mctk);
          if(matchedTrackPairInfo) {
            StGlobalTrack *gRctrk = (StGlobalTrack*) matchedTrackPairInfo->partnerTrack();
            if(gRctrk) {StPrimaryTrack* rctk = dynamic_cast<StPrimaryTrack *>(gRctrk->node()->track(primary));
//          tpcCommonHits = matchedTrackPairInfo->commonTpcHits();
//	      if(rctk)rctk->Print("");else {cout<<"Global ";gRctrk->Print("");} 
	      if(rctk){
		float rcpt=rctk->geometry()->momentum().perp();
		float rceta=rctk->geometry()->momentum().pseudoRapidity();
//		cout<<"Primary Momentum "<<rctk->geometry()->momentum().perp()<<endl;
		mResponseMtx->Fill(mctk->pt(),rcpt);
		mRcPt->Fill(rcpt);
		mRcPt_eta->Fill(rcpt,rceta);
		Match_Pt_eta->Fill(mctk->pt(),mctk->momentum().pseudoRapidity());
		mTrkEff->Fill(mctk->pt());
      		mRceta_phi->Fill(rctk->geometry()->momentum().pseudoRapidity(),rctk->geometry()->momentum().phi()); 
      		Match_eta_phi->Fill(mctk->momentum().pseudoRapidity(),mctk->momentum().phi()); 
//		cout<<" found"<<endl;
	      }
//	    else cout<<"\n";
//	      else
//		mTrkEff->Fill(mctk->pt(),0);
//	      else
//		cout<<"Global Momentum "<<gRctrk->geometry()->momentum().perp()<<endl;
	    }
//	    else cout<<"\n\n";
	  }
//	  else cout<<"Track not matched. \n"<<endl;
        }
//       else cout<<"No map\n"<<endl;
      }
//     else cout<<"\n\n\n"; 
    }
  }
  else{
    cout<<"There is no primary vertex in the event!妈妈服务器欺负我！"<<endl;
    cout<<"If you are wondering, the Chinese means I am mad."<<endl;
  }
  for (rcTrackMapIter tIter=theTrackMap->begin();
       tIter!=theTrackMap->end(); ++tIter){
    
    recTrack = (*tIter).first;
    //yf    if ((*tIter).second->commonTpcHits()<10) continue;
    mcTrack = (*tIter).second->partnerMcTrack();
    pmc = mcTrack->momentum();
    for (int k=0; k<3; k++) values[k] = pmc[k];  	 
    values[3]=pmc.mag(); 	 
    
    primTrk = dynamic_cast<const StPrimaryTrack*>(recTrack->node()->track(primary)); 	 
    if (primTrk) 	 
	p = primTrk->geometry()->momentum(); 	 
    else 	 
	p = recTrack->geometry()->momentum(); 	 
    for (int j=0; j<3; j++) values[j+4] = p[j]; 	 
    values[7]=p.mag(); 	 
    values[8]=(*tIter).second->commonTpcHits();
    // Fill 1d Mom. resolution Histogram
    diff = (p.mag() - pmc.mag())/pmc.mag();
    mMomResolution->Fill(diff,1.);
    // Loop to get Mean hit position diff.
    StThreeVectorF rHitPos(0,0,0);
    StThreeVectorF mHitPos(0,0,0);
    StPtrVecHit recTpcHits = recTrack->detectorInfo()->hits(kTpcId);
    multimap<int,float> idTruths;
    set<int> uniqueIdTruths;
    for (StHitIterator hi=recTpcHits.begin();
	 hi!=recTpcHits.end(); hi++) {
      StHit* hit = *hi;
      StTpcHit* rHit = dynamic_cast<StTpcHit*>(hit);
      if (!rHit) { cout << "This Hit is not a TPC Hit"<< endl; continue;}
      idTruths.insert( multimap<int,float>::value_type(rHit->idTruth(),rHit->qaTruth()));
      uniqueIdTruths.insert(static_cast<int>(rHit->idTruth()));
      pair<rcTpcHitMapIter,rcTpcHitMapIter> rBounds = theHitMap->equal_range(rHit);
      for (rcTpcHitMapIter hIter=rBounds.first; hIter!=rBounds.second; hIter++) {
	const StMcTpcHit* mHit = (*hIter).second;
	if (mHit->parentTrack() != mcTrack) continue;
	rHitPos += rHit->position();
	mHitPos += mHit->position();
      }// Associated Hits Loop.
      
    } // Hits of rec. Track Loop
    rHitPos /=(float) (*tIter).second->commonTpcHits();
    mHitPos /=(float) (*tIter).second->commonTpcHits();
    for (int jj=0; jj<3; jj++) values[9+jj] = rHitPos[jj] - mHitPos[jj];
    values[12] = mcTrack->key();
    // Figure out the most common IdTruth; the dominatrix track!
    int mostCommonIdTruth = -9; 
    int cachedNHitsIdTruth = 0;
    for (set<int>::iterator si=uniqueIdTruths.begin(); si!=uniqueIdTruths.end(); ++si) {
	int currentNHitsIdTruth = idTruths.count(static_cast<int>(*si));
	if (currentNHitsIdTruth>cachedNHitsIdTruth) {
	    mostCommonIdTruth = *si; 
	    cachedNHitsIdTruth = currentNHitsIdTruth;
	}
    }
    // at this point we know the most common IdTruth,
    // now calculate the track "quality" for this track, averaging
    // the hit qualities
    float idQuality = 0;
    pair<multimap<int,float>::iterator,multimap<int,float>::iterator> mostCommRange = idTruths.equal_range(mostCommonIdTruth);
    for (multimap<int,float>::iterator mi=mostCommRange.first; mi!=mostCommRange.second; ++mi) {
	idQuality+=mi->second;
    }
    idQuality/=cachedNHitsIdTruth;
    values[13] = mostCommonIdTruth;
    values[14] = cachedNHitsIdTruth;
    values[15] = mcTrack->tpcHits().size();
    values[16] = recTrack->fitTraits().numberOfFitPoints(kTpcId);
    values[17] = recTpcHits.size();
    values[18] = idQuality;
//    mTrackNtuple->Fill(values);
  } // Tracks in Map Loop
//  cout << "Finished Track Loop, Made Ntuple" << endl;
  //delete vars;
  delete [] values;
  
  // Example: Make 2 Histograms
  // - x and y positions of the hits from the reconstructed track.
  // - x and y positions of the hits from the  Monte Carlo  track.
  
  unsigned int maxCommonTpcHits = 0;
  const StMcTrack* partner = 0;
  trackBounds = theTrackMap->equal_range(firstTrack);
  for (rcTrackMapIter rcIt = trackBounds.first;
       rcIt != trackBounds.second;
       ++rcIt) {
    if ((*rcIt).second->commonTpcHits() >  maxCommonTpcHits) {
      partner = (*rcIt).second->partnerMcTrack();
      maxCommonTpcHits = (*rcIt).second->commonTpcHits();
    }
  }
  StHitIterator rcHitIt;
  StMcTpcHitIterator mcHitIt;
  StMcSvtHitIterator mcSitIt;
  StPtrVecHit theHits = firstTrack->detectorInfo()->hits();
  for (rcHitIt  = theHits.begin();
       rcHitIt != theHits.end();
       rcHitIt++) coordRec->Fill((*rcHitIt)->position().x(),(*rcHitIt)->position().y());
  if (partner) {
    for (mcHitIt  = ((std::vector<StMcTpcHit*> *)&partner->tpcHits() )->begin();
	 mcHitIt != partner->tpcHits().end();
	 mcHitIt++) coordMcPartner->Fill((*mcHitIt)->position().x(),(*mcHitIt)->position().y());
    for (mcSitIt  = ((std::vector<StMcSvtHit*> *)&partner->svtHits())->begin();
	 mcSitIt != partner->svtHits().end();
	 mcSitIt++) coordMcPartner->Fill((*mcSitIt)->position().x(),(*mcSitIt)->position().y());
    
  }
  if (!theMcV0Map) {
    gMessMgr->Warning() << "----------WARNING----------\n"
			<< "No V0 Map found for this event!" << endm;
    return kStWarn;
  }
  //Example: Print out position of V0 vertices that have been associated.
  //         (LSB)
  StSPtrVecMcVertex& mcVertices = mEvent->vertices();
  const StV0Vertex* rcV0Partner;
  StMcVertexIterator mcVertexIt;
  
  //Loop over all MC vertices
  bool foundV0Pair = false;
  for (mcVertexIt = mcVertices.begin(); mcVertexIt != mcVertices.end();
       mcVertexIt++){
    // Get the upper and lower bounds.
    pair<mcV0MapIter,mcV0MapIter> mcV0Bounds = theMcV0Map->equal_range(*mcVertexIt);
    
    // Print out MC vertex position if there is an associated V0.
    
    if (mcV0Bounds.first != mcV0Bounds.second) {
      cout << "Printing Position of a V0 pair:\n";
      cout << "Position of MC V0 vertex: " << (*mcVertexIt)->position() << endl;
      foundV0Pair = true;
    }
    //Now loop over the bounds      
    for(mcV0MapIter mcV0MapIt = mcV0Bounds.first;
	mcV0MapIt != mcV0Bounds.second; ++mcV0MapIt){
      rcV0Partner = (*mcV0MapIt).second;
      cout << "Position of rc V0 vertex: " << rcV0Partner->position() << endl;
      
    }
    if (foundV0Pair) break; // Only print the information of 1 reconstructed vertex, to avoid a lot of output.
  }
  
  
  
  //mAssociationCanvas = new TCanvas("mAssociationCanvas", "Histograms",200,10,900,500);
  
  return kStOK;
}
//=======================================================================
StTrackPairInfo* StMyMcAnalysisMaker::findBestMatchedGlobal(StMcTrack* mcTrack)
{
  pair<mcTrackMapIter,mcTrackMapIter> mcBounds
    = mMcTrackMap->equal_range(mcTrack);
  StTrackPairInfo* candTrackPair = 0;  // used for finding the best matched track
  StGlobalTrack* candTrack = 0;
  mcTrackMapIter mcMapIter = mcBounds.first;
  for ( ; mcMapIter != mcBounds.second; ++mcMapIter){
    StTrackPairInfo* assocPair = (*mcMapIter).second;
    StGlobalTrack* globTrack =(StGlobalTrack*) assocPair->partnerTrack();
    if (Debug() > 1) {
      cout << * assocPair << endl;
      cout << "globTrack FitPoints Tpc/FtpcE/W = " << globTrack->fitTraits().numberOfFitPoints(kTpcId)
        << "/" << globTrack->fitTraits().numberOfFitPoints(kFtpcEastId)
        << "/" << globTrack->fitTraits().numberOfFitPoints(kFtpcWestId) << endl;
    }
    if (!globTrack || globTrack->flag()<=0) continue;
    if (globTrack->fitTraits().numberOfFitPoints(kTpcId)>=10 ||
        globTrack->fitTraits().numberOfFitPoints(kFtpcEastId)>=5 ||
        globTrack->fitTraits().numberOfFitPoints(kFtpcWestId)>=5) {
      if (!candTrackPair) {
        candTrackPair = assocPair;
        candTrack = globTrack;
      }
      else if (globTrack->fitTraits().numberOfFitPoints(kTpcId) > candTrack->fitTraits().numberOfFitPoints(kTpcId)) {
        candTrackPair = assocPair;
        candTrack = globTrack;
      }
      else if (globTrack->fitTraits().numberOfFitPoints(kFtpcEastId) > candTrack->fitTraits().numberOfFitPoints(kFtpcEastId)) {
        candTrackPair = assocPair;
        candTrack = globTrack;
      }
      else if (globTrack->fitTraits().numberOfFitPoints(kFtpcWestId) > candTrack->fitTraits().numberOfFitPoints(kFtpcWestId)) {
        candTrackPair = assocPair;
        candTrack = globTrack;
      }

    } // fit points requirement
  }// bounds loop
  return candTrackPair; // Note that candTrack might be zero, for example if only one track is matched and has 9 tpc fit pts.
}
