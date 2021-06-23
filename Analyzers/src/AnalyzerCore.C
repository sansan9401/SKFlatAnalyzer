#include "AnalyzerCore.h"

AnalyzerCore::AnalyzerCore(){

  outfile = NULL;
  mcCorr = new MCCorrection();
  puppiCorr = new PuppiSoftdropMassCorr();
  fakeEst = new FakeBackgroundEstimator();
  cfEst = new CFBackgroundEstimator();
  pdfReweight = new PDFReweight();
  muonGE = new GeneralizedEndpoint();
  muonGEScaleSyst = new GEScaleSyst();

}

AnalyzerCore::~AnalyzerCore(){

  //=== hist maps

  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH1D.clear();

  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH2D.clear();

  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH3D.clear();
  
  //==== output rootfile

  if(outfile){
    outfile->Close();
    delete outfile;
  }

  //==== Tools

  if(mcCorr) delete mcCorr;
  if(puppiCorr) delete puppiCorr;
  if(fakeEst) delete fakeEst;
  if(cfEst) delete cfEst;
  if(pdfReweight) delete pdfReweight;
  if(muonGE) delete muonGE;
  if(muonGEScaleSyst) delete muonGEScaleSyst;

}

//==== Attach the historams to ai different direcotry, not outfile
//==== We will write these histograms in WriteHist() to outfile
void AnalyzerCore::SwitchToTempDir(){

  gROOT->cd();
  TDirectory *tempDir = NULL;
  int counter = 0;
  while (!tempDir) {
    //==== First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "AnalyzerCore" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    //==== Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());
  }
  tempDir->cd();

}

void AnalyzerCore::SetOutfilePath(TString outname){
  outfile = new TFile(outname,"RECREATE");
};

Event AnalyzerCore::GetEvent(){

  Event ev;
  if(!IsDATA) ev.SetMCweight(gen_weight);
  ev.SetTrigger(*HLT_TriggerName);
  ev.SetMET(pfMET_Type1_pt,pfMET_Type1_phi);
  ev.SetnPV(nPV);
  ev.SetEra(GetEra());

  return ev;

}

std::vector<Muon> AnalyzerCore::GetAllMuons(){

  std::vector<Muon> out;
  if(!muon_pt) return out;
  for(unsigned int i=0; i<muon_pt->size(); i++){

    Muon mu;

    mu.SetCharge(muon_charge->at(i));
    mu.SetMiniAODPt(muon_pt->at(i));
    mu.SetMiniAODTunePPt(muon_TuneP_pt->at(i));

    double rc = muon_roch_sf->at(i);
    double rc_err = muon_roch_sf_up->at(i)-rc;
    //==== For the Rochester corection, up and down err are the same
    mu.SetMomentumScaleUpDown( muon_pt->at(i) * (rc+rc_err), muon_pt->at(i) * (rc-rc_err)  );
    mu.SetPtEtaPhiM(muon_pt->at(i)*rc, muon_eta->at(i), muon_phi->at(i), muon_mass->at(i));

    //==== TuneP
    //==== Apply scailing later with AnalyzerCore::UseTunePMuon()
    mu.SetTuneP4(muon_TuneP_pt->at(i), muon_TuneP_ptError->at(i), muon_TuneP_eta->at(i), muon_TuneP_phi->at(i), muon_TuneP_charge->at(i));

    mu.SetMVA(muon_MVA->at(i));

    mu.SetdXY(muon_dxyVTX->at(i), muon_dxyerrVTX->at(i));
    mu.SetdZ(muon_dzVTX->at(i), muon_dzerrVTX->at(i));
    mu.SetIP3D(muon_3DIPVTX->at(i), muon_3DIPerrVTX->at(i));
    mu.SetTypeBit(muon_TypeBit->at(i));
    mu.SetIDBit(muon_IDBit->at(i));
    mu.SetisPOGHighPt(muon_ishighpt->at(i));
    mu.SetChi2(muon_normchi->at(i));
    mu.SetIso(muon_PfChargedHadronIsoR04->at(i),muon_PfNeutralHadronIsoR04->at(i),muon_PfGammaIsoR04->at(i),muon_PFSumPUIsoR04->at(i),muon_trkiso->at(i));
    mu.SetTrackerLayers(muon_trackerLayers->at(i));

    //==== Should be set after Eta is set
    mu.SetMiniIso(
      muon_PfChargedHadronMiniIso->at(i), 
      muon_PfNeutralHadronMiniIso->at(i), 
      muon_PfGammaMiniIso->at(i), 
      muon_PFSumPUMiniIso->at(i),
      Rho,
      mu.EA()
    );

    out.push_back(mu);

  }
  return out;

}

std::vector<Muon> AnalyzerCore::GetMuons(TString id, double ptmin, double fetamax){

  std::vector<Muon> muons = GetAllMuons();
  std::vector<Muon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    if(!( muons.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << muons.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(muons.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(muons.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( muons.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( muons.at(i) );
  }
  return out;

}

std::vector<Electron> AnalyzerCore::GetAllElectrons(){

  std::vector<Electron> out;
  if(!electron_Energy) return out;
  for(unsigned int i=0; i<electron_Energy->size(); i++){

    Electron el;

    el.SetEnShift(  electron_Energy_Scale_Up->at(i)/electron_Energy->at(i), electron_Energy_Scale_Down->at(i)/electron_Energy->at(i) );
    el.SetResShift( electron_Energy_Smear_Up->at(i)/electron_Energy->at(i), electron_Energy_Smear_Down->at(i)/electron_Energy->at(i) );

    el.SetPtEtaPhiE(1., electron_eta->at(i), electron_phi->at(i), electron_Energy->at(i));
    double el_theta = el.Theta();
    double el_pt = electron_Energy->at(i) * TMath::Sin( el_theta );
    el.SetPtEtaPhiE( el_pt, electron_eta->at(i), electron_phi->at(i), electron_Energy->at(i));

    el.SetUncorrE(electron_EnergyUnCorr->at(i));
    el.SetSC(electron_scEta->at(i), electron_scPhi->at(i), electron_scEnergy->at(i));
    el.SetCharge(electron_charge->at(i));
    el.SetdXY(electron_dxyVTX->at(i), electron_dxyerrVTX->at(i));
    el.SetdZ(electron_dzVTX->at(i), electron_dzerrVTX->at(i));
    el.SetIP3D(electron_3DIPVTX->at(i), electron_3DIPerrVTX->at(i));
    el.SetMVA(electron_MVAIso->at(i), electron_MVANoIso->at(i));
    el.SetPassConversionVeto(electron_passConversionVeto->at(i));
    el.SetNMissingHits(electron_mHits->at(i));
    el.SetRho(Rho);
    el.SetIsGsfCtfScPixChargeConsistent(electron_isGsfCtfScPixChargeConsistent->at(i));
    el.SetR9(electron_r9->at(i));

    el.SetCutBasedIDVariables(
      electron_Full5x5_SigmaIEtaIEta->at(i),
      electron_dEtaInSeed->at(i),
      electron_dPhiIn->at(i),
      electron_HoverE->at(i),
      electron_InvEminusInvP->at(i),
      electron_e2x5OverE5x5->at(i),
      electron_e1x5OverE5x5->at(i),
      electron_trackIso->at(i),
      electron_dr03EcalRecHitSumEt->at(i),
      electron_dr03HcalDepth1TowerSumEt->at(i),
      electron_dr03HcalTowerSumEt->at(i),
      electron_dr03TkSumPt->at(i),
      electron_ecalPFClusterIso->at(i),
      electron_hcalPFClusterIso->at(i),
      electron_ecalDriven->at(i)
    );

    el.SetIDBit(electron_IDBit->at(i));
    vector<int> temp_idcutbit;
    for(unsigned int j=0; j<Electron::N_SELECTOR; j++){
      temp_idcutbit.push_back( electron_IDCutBit->at( i*Electron::N_SELECTOR + j ) );
    }
    el.SetIDCutBit(temp_idcutbit);
    el.SetRelPFIso_Rho(electron_RelPFIso_Rho->at(i));

    //==== Should be ran after SCeta is set
    el.SetMiniIso(
      electron_chMiniIso->at(i),
      electron_nhMiniIso->at(i),
      electron_phMiniIso->at(i),
      electron_puChMiniIso->at(i),
      Rho,
      el.EA()
    );

    out.push_back(el);

  }
  return out;

}

std::vector<Electron> AnalyzerCore::GetElectrons(TString id, double ptmin, double fetamax){

  std::vector<Electron> electrons = GetAllElectrons();
  std::vector<Electron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    if(!( electrons.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << electrons.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(electrons.at(i).scEta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(electrons.at(i).scEta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( electrons.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( electrons.at(i) );
  }
  return out;

}

std::vector<Lepton *> AnalyzerCore::MakeLeptonPointerVector(const std::vector<Muon>& muons, double TightIso, bool UseMini){

  std::vector<Lepton *> out;
  for(unsigned int i=0; i<muons.size(); i++){
    Lepton *l = (Lepton *)(&muons.at(i));
    if( !(l->LeptonFlavour() == Lepton::MUON) ){
      cout << "[AnalyzerCore::MakeLeptonPointerVector(std::vector<Muon>& muons)] Not muon.." << endl;
      exit(EXIT_FAILURE);
    }
    if(TightIso>0){

      double this_RelIso = l->RelIso();
      if(UseMini) this_RelIso = l->MiniRelIso();
      double ptcone = l->CalcPtCone(this_RelIso, TightIso);
      l->SetPtCone( ptcone );

    }
    out.push_back(l);
  }
  return out;

}
std::vector<Lepton *> AnalyzerCore::MakeLeptonPointerVector(const std::vector<Electron>& electrons, double TightIso, bool UseMini){

  std::vector<Lepton *> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    Lepton *l = (Lepton *)(&electrons.at(i));
    if( !(l->LeptonFlavour() == Lepton::ELECTRON) ){
      cout << "[AnalyzerCore::MakeLeptonPointerVector(std::vector<ELECTRON>& electrons)] Not electron.." << endl;
      exit(EXIT_FAILURE);
    }
    if(TightIso>0){

      double this_RelIso = l->RelIso();
      if(UseMini) this_RelIso = l->MiniRelIso();
      double ptcone = l->CalcPtCone(this_RelIso, TightIso);
      l->SetPtCone( ptcone );

    }
    out.push_back(l);
  }
  return out;

}



std::vector<Photon> AnalyzerCore::GetAllPhotons(){

  std::vector<Photon> out;
  for(unsigned int i=0; i<photon_Energy->size(); i++){
    
    Photon pho;

    pho.SetPtEtaPhiE(1., photon_eta->at(i), photon_phi->at(i), photon_Energy->at(i));
    double pho_theta = pho.Theta();
    double pho_pt = photon_Energy->at(i) * TMath::Sin( pho_theta );
    pho.SetPtEtaPhiE( pho_pt, photon_eta->at(i), photon_phi->at(i), photon_Energy->at(i));

    pho.SetSC(photon_scEta->at(i), photon_scPhi->at(i));
    pho.SetRho(Rho);

    pho.SetCutBasedIDVariables(
      photon_Full5x5_SigmaIEtaIEta->at(i),
      photon_HoverE->at(i),
      photon_ChIsoWithEA->at(i),
      photon_NhIsoWithEA->at(i),
      photon_PhIsoWithEA->at(i)
    );
    
    std::vector<bool> ids = {
      photon_passLooseID->at(i),
      photon_passMediumID->at(i),
      photon_passTightID->at(i),
      photon_passMVAID_WP80->at(i),
      photon_passMVAID_WP90->at(i)
     };
    pho.SetPOGIDs(ids);
        
    out.push_back(pho);
    
  }
  return out;
  
}

std::vector<Photon> AnalyzerCore::GetPhotons(TString id, double ptmin, double fetamax){

  std::vector<Photon> photons = GetAllPhotons();
  std::vector<Photon> out;
  for(unsigned int i=0; i<photons.size(); i++){
    if(!( photons.at(i).Pt()>ptmin )){
      continue;
    }
    if(!( fabs(photons.at(i).scEta())<fetamax )){
      continue;
    }
    if(!( photons.at(i).PassID(id) )){
      continue;
    }
    out.push_back( photons.at(i) );
  }
  return out;
}



std::vector<Jet> AnalyzerCore::GetAllJets(){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jet_pt->size(); i++){
    Jet jet;
    jet.SetPtEtaPhiM(jet_pt->at(i), jet_eta->at(i), jet_phi->at(i), jet_m->at(i));

    //==== Jet energy up and down are 1.xx or 0.99, not energy
    jet.SetEnShift( jet_shiftedEnUp->at(i), jet_shiftedEnDown->at(i) );
    if(!IsDATA){
      jet *= jet_smearedRes->at(i);
      jet.SetResShift( jet_smearedResUp->at(i)/jet_smearedRes->at(i), jet_smearedResDown->at(i)/jet_smearedRes->at(i) );
    }
    jet.SetCharge(jet_charge->at(i));

    jet.SetArea(jet_area->at(i));
    jet.SetGenFlavours(jet_partonFlavour->at(i), jet_hadronFlavour->at(i));
    std::vector<double> tvs = {
      jet_CSVv2->at(i),
      jet_DeepCSV->at(i),
      jet_DeepCvsL->at(i),
      jet_DeepCvsB->at(i),
      jet_DeepFlavour_b->at(i),
      jet_DeepFlavour_bb->at(i),
      jet_DeepFlavour_lepb->at(i),
      jet_DeepFlavour_c->at(i),
      jet_DeepFlavour_uds->at(i),
      jet_DeepFlavour_g->at(i),
      jet_CvsL->at(i),
      jet_CvsB->at(i),
    };
    jet.SetTaggerResults(tvs);
    jet.SetEnergyFractions(jet_chargedHadronEnergyFraction->at(i), jet_neutralHadronEnergyFraction->at(i), jet_neutralEmEnergyFraction->at(i), jet_chargedEmEnergyFraction->at(i), jet_muonEnergyFraction->at(i));
    jet.SetMultiplicities(jet_chargedMultiplicity->at(i), jet_neutralMultiplicity->at(i));
    jet.SetPileupJetId(jet_PileupJetId->at(i));
    jet.SetTightJetID(jet_tightJetID->at(i));
    jet.SetTightLepVetoJetID(jet_tightLepVetoJetID->at(i));

    out.push_back(jet);
  }

  return out;

}

std::vector<Jet> AnalyzerCore::GetJets(TString id, double ptmin, double fetamax){

  std::vector<Jet> jets = GetAllJets();
  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    if(!( jets.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << jets.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(jets.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(jets.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( jets.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( jets.at(i) );
  }
  return out;

}

std::vector<FatJet> AnalyzerCore::GetAllFatJets(){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<fatjet_pt->size(); i++){
    FatJet jet;
    jet.SetPtEtaPhiM(fatjet_pt->at(i), fatjet_eta->at(i), fatjet_phi->at(i), fatjet_m->at(i));
    jet.SetEnShift( fatjet_shiftedEnUp->at(i), fatjet_shiftedEnDown->at(i) );
    if(!IsDATA){
      jet *= fatjet_smearedRes->at(i);
      jet.SetResShift( fatjet_smearedResUp->at(i)/fatjet_smearedRes->at(i), fatjet_smearedResDown->at(i)/fatjet_smearedRes->at(i) );
    }
    jet.SetCharge(fatjet_charge->at(i));

    jet.SetArea(fatjet_area->at(i));
    jet.SetGenFlavours(fatjet_partonFlavour->at(i), fatjet_hadronFlavour->at(i));
    std::vector<double> tvs = {
      fatjet_CSVv2->at(i),
      fatjet_DeepCSV->at(i),
      fatjet_DeepCvsL->at(i),
      fatjet_DeepCvsB->at(i),
      fatjet_DeepFlavour_b->at(i),
      fatjet_DeepFlavour_bb->at(i),
      fatjet_DeepFlavour_lepb->at(i),
      fatjet_DeepFlavour_c->at(i),
      fatjet_DeepFlavour_uds->at(i),
      fatjet_DeepFlavour_g->at(i),
      fatjet_CvsL->at(i),
      fatjet_CvsB->at(i),
    };
    jet.SetTaggerResults(tvs);
    jet.SetEnergyFractions(fatjet_chargedHadronEnergyFraction->at(i), fatjet_neutralHadronEnergyFraction->at(i), fatjet_neutralEmEnergyFraction->at(i), fatjet_chargedEmEnergyFraction->at(i), fatjet_muonEnergyFraction->at(i));
    jet.SetMultiplicities(fatjet_chargedMultiplicity->at(i), fatjet_neutralMultiplicity->at(i));
    jet.SetLSF(fatjet_LSF->at(i), fatjet_LSFlep_PID->at(i));
    jet.SetTightJetID(fatjet_tightJetID->at(i));
    jet.SetTightLepVetoJetID(fatjet_tightLepVetoJetID->at(i));
    jet.SetPuppiTaus(fatjet_puppi_tau1->at(i), fatjet_puppi_tau2->at(i), fatjet_puppi_tau3->at(i), fatjet_puppi_tau4->at(i));
    jet.SetSDMass(fatjet_softdropmass->at(i));

    out.push_back(jet);
  }

  return out;

}

std::vector<FatJet> AnalyzerCore::GetFatJets(TString id, double ptmin, double fetamax){

  std::vector<FatJet> jets = GetAllFatJets();
  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    if(!( jets.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << jets.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(jets.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(jets.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( jets.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( jets.at(i) );
  }
  return out;

}

std::vector<Gen> AnalyzerCore::GetGens(){

  std::vector<Gen> out;
  if(IsDATA) return out;

  for(unsigned int i=0; i<gen_pt->size(); i++){

    Gen gen;

    gen.SetIsEmpty(false);
    gen.SetPtEtaPhiM(gen_pt->at(i), gen_eta->at(i), gen_phi->at(i), gen_mass->at(i));
    gen.SetCharge(gen_charge->at(i));
    gen.SetIndexPIDStatus(i, gen_PID->at(i), gen_status->at(i));
    gen.SetMother(gen_mother_index->at(i));
    gen.SetGenStatusFlag_isPrompt( gen_isPrompt->at(i) );
    gen.SetGenStatusFlag_isPromptFinalState( gen_isPromptFinalState->at(i) );
    gen.SetGenStatusFlag_isTauDecayProduct( gen_isTauDecayProduct->at(i) );
    gen.SetGenStatusFlag_isPromptTauDecayProduct( gen_isPromptTauDecayProduct->at(i) );
    gen.SetGenStatusFlag_isDirectPromptTauDecayProductFinalState( gen_isDirectPromptTauDecayProductFinalState->at(i) );
    gen.SetGenStatusFlag_isHardProcess( gen_isHardProcess->at(i) );
    gen.SetGenStatusFlag_isLastCopy( gen_isLastCopy->at(i) );
    gen.SetGenStatusFlag_isLastCopyBeforeFSR( gen_isLastCopyBeforeFSR->at(i) );
    gen.SetGenStatusFlag_isPromptDecayed( gen_isPromptDecayed->at(i) );
    gen.SetGenStatusFlag_isDecayedLeptonHadron( gen_isDecayedLeptonHadron->at(i) );
    gen.SetGenStatusFlag_fromHardProcessBeforeFSR( gen_fromHardProcessBeforeFSR->at(i) );
    gen.SetGenStatusFlag_fromHardProcessDecayed( gen_fromHardProcessDecayed->at(i) );
    gen.SetGenStatusFlag_fromHardProcessFinalState( gen_fromHardProcessFinalState->at(i) );
    gen.SetGenStatusFlag_isMostlyLikePythia6Status3( gen_isMostlyLikePythia6Status3->at(i) );

    out.push_back(gen);

  }

  return out;

}

std::vector<LHE> AnalyzerCore::GetLHEs(){

  std::vector<LHE> out;
  if(IsDATA) return out;

  for(unsigned int i=0; i<LHE_Px->size(); i++){

    LHE lhe;

    lhe.SetPxPyPzE(LHE_Px->at(i), LHE_Py->at(i), LHE_Pz->at(i), LHE_E->at(i));
    lhe.SetIndexIDStatus(i, LHE_ID->at(i), LHE_Status->at(i));

    out.push_back(lhe);

  }

  return out;

}

std::vector<Muon> AnalyzerCore::UseTunePMuon(const std::vector<Muon>& muons){

  std::vector<Muon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    //==== muons is a const vector. So in this function, we have to copy the elements like below
    Muon this_muon=muons.at(i);

    Particle this_tunep4 = this_muon.TuneP4();

    //==== Momentum scaling
    //==== 1) if tuneP Pt < 200 -> Rochester
    //==== 2) if tuneP pt >= 200 -> Generalized Endpoint

    double new_pt( this_tunep4.Pt() ), new_pt_up( this_tunep4.Pt() ), new_pt_down( this_tunep4.Pt() );
    if(this_tunep4.Pt()<200){

      //==== 19/03/24 (jskim) : For 99% of the muons, MiniAODPt and TunePPt are same
      //==== we can just use MiniAODPt * RochesterCorrection, multiplied by (TuneP Pt)/(MiniAODPt)
      double TunePOverPt = this_tunep4.Pt() / this_muon.MiniAODPt();
      new_pt      = TunePOverPt * this_muon.Pt(); // this_muon.Pt() = MiniAODPt * RochesterCorrection
      new_pt_up   = TunePOverPt * this_muon.MomentumShift(+1);
      new_pt_down = TunePOverPt * this_muon.MomentumShift(-1);

/*
      cout << "## Rochester ##" << endl;
      cout << "this_muon.MiniAODPt() = " << this_muon.MiniAODPt() << endl;
      cout << "this_muon.MiniAODTunePPt() = " << this_muon.MiniAODTunePPt() << endl;
      cout << "new_pt = " << new_pt << endl;
      cout << "new_pt_up = " << new_pt_up << endl;
      cout << "new_pt_down = " << new_pt_down << endl;
*/

    }
    else{

      //==== Unlike rochester, GE method should be only applied to MC

      if(!IsDATA){

        //==== ScaledPts defined in GeneralizedEndpointPt.h ..

        ScaledPts ptvalues;
        //==== TODO FIXME
        //==== 19/09/02 : There is no GEScaleSyst map for 2016
        if(DataYear==2016) ptvalues = muonGE->GeneralizedEndpointPt(this_tunep4.Pt(), this_tunep4.Charge(), this_tunep4.Eta(), this_tunep4.Phi()*180./M_PI, event);
        else ptvalues = muonGEScaleSyst->GEPt(DataYear, this_tunep4.Pt(), this_tunep4.Eta(), this_tunep4.Phi(), this_tunep4.Charge());

        new_pt = ptvalues.ScaledPt;
        //==== Mode == 1 : Kappa up
        //==== Mode == 2 : Kappa down
        new_pt_up = ptvalues.ScaeldPt_Up;
        new_pt_down = ptvalues.ScaeldPt_Down;

/*
        cout << "## GeneralizedEndpointPt ##" << endl;
        cout << "old_pt = " << this_tunep4.Pt() << endl;
        cout << "new_pt = " << new_pt << endl;
        cout << "new_pt_up = " << new_pt_up << endl;
        cout << "new_pt_down = " << new_pt_down << endl;
*/

      }

    }

    //==== Scale the pt
    this_muon.SetPtEtaPhiM( new_pt, this_tunep4.Eta(), this_tunep4.Phi(), this_tunep4.M() );
    this_muon.SetMomentumScaleUpDown(new_pt_up,new_pt_down);
    this_muon.SetCharge( this_tunep4.Charge() );
    this_muon.SetMiniAODPt( this_muon.MiniAODTunePPt() );

/*
    cout << "@@@@ TuneP @@@@" << endl;
    cout << "this_muon.Pt() = " << this_muon.Pt() << endl;
    cout << "this_muon.MiniAODPt() = " << this_muon.MiniAODPt() << endl;
    cout << "this_muon.MiniAODTunePPt() = " << this_muon.MiniAODTunePPt() << endl;
    cout << "this_muon.MomentumShift(0) = " << this_muon.MomentumShift(0) << endl;
    cout << "this_muon.MomentumShift(+1) = " << this_muon.MomentumShift(+1) << endl;
    cout << "this_muon.MomentumShift(-1) = " << this_muon.MomentumShift(-1) << endl;
*/

    out.push_back(this_muon);
  }

  return out;

}

std::vector<Muon> AnalyzerCore::SelectMuons(const std::vector<Muon>& muons, TString id, double ptmin, double fetamax){

  std::vector<Muon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    if(!( muons.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << muons.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(muons.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(muons.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( muons.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( muons.at(i) );
  }
  return out;

}

std::vector<Electron> AnalyzerCore::SelectElectrons(const std::vector<Electron>& electrons, TString id, double ptmin, double fetamax){

  std::vector<Electron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    if(!( electrons.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << electrons.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(electrons.at(i).scEta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(electrons.at(i).scEta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( electrons.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back(electrons.at(i));
  }
  return out;

}

std::vector<Jet> AnalyzerCore::SelectJets(const std::vector<Jet>& jets, TString id, double ptmin, double fetamax){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    if(!( jets.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << jets.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(jets.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(jets.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( jets.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( jets.at(i) );
  }
  return out;

}

std::vector<FatJet> AnalyzerCore::SelectFatJets(const std::vector<FatJet>& jets, TString id, double ptmin, double fetamax){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    if(!( jets.at(i).Pt()>ptmin )){
      //cout << "Fail Pt : pt = " << jets.at(i).Pt() << ", cut = " << ptmin << endl;
      continue;
    }
    if(!( fabs(jets.at(i).Eta())<fetamax )){
      //cout << "Fail Eta : eta = " << fabs(jets.at(i).Eta()) << ", cut = " << fetamax << endl;
      continue;
    }
    if(!( jets.at(i).PassID(id) )){
      //cout << "Fail ID" << endl;
      continue;
    }
    out.push_back( jets.at(i) );
  }
  return out;

}

std::vector<Electron> AnalyzerCore::ScaleElectrons(const std::vector<Electron>& electrons, int sys){

  std::vector<Electron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    //==== electrons is a const vector. So in this function, we have to copy the elements like below
    Electron this_electron = electrons.at(i);

    double this_sf = this_electron.EnShift(sys);
    this_electron.SetPtEtaPhiM( this_electron.Pt() * this_sf, this_electron.Eta(), this_electron.Phi(), this_electron.M() );

    out.push_back( this_electron );
  }

  return out;

}
std::vector<Electron> AnalyzerCore::SmearElectrons(const std::vector<Electron>& electrons, int sys){

  std::vector<Electron> out;
  for(unsigned int i=0; i<electrons.size(); i++){
    //==== electrons is a const vector. So in this function, we have to copy the elements like below
    Electron this_electron = electrons.at(i);

    double this_sf = this_electron.ResShift(sys);
    this_electron.SetPtEtaPhiM( this_electron.Pt() * this_sf, this_electron.Eta(), this_electron.Phi(), this_electron.M() );

    out.push_back( this_electron );
  }

  return out;

}

std::vector<Muon> AnalyzerCore::ScaleMuons(const std::vector<Muon>& muons, int sys){

  std::vector<Muon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    //==== muons is a const vector. So in this function, we have to copy the elements like below
    Muon this_muon = muons.at(i);

    //==== Even for TuneP muons, MomentumShift() are set correctly from AnalyzerCore::UseTunePMuon()
    //==== So we can just use MomentumShift()

    this_muon.SetPtEtaPhiM( this_muon.MomentumShift(sys), this_muon.Eta(), this_muon.Phi(), this_muon.M() );

    out.push_back(this_muon);

  }

  return out;

}

std::vector<Jet> AnalyzerCore::ScaleJets(const std::vector<Jet>& jets, int sys){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    Jet this_jet = jets.at(i);

    this_jet *= this_jet.EnShift(sys);

    out.push_back( this_jet );
  }

  return out;

}
std::vector<Jet> AnalyzerCore::SmearJets(const std::vector<Jet>& jets, int sys){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    Jet this_jet = jets.at(i);

    this_jet *= this_jet.ResShift(sys);

    out.push_back( this_jet );
  }

  return out;

}

std::vector<FatJet> AnalyzerCore::ScaleFatJets(const std::vector<FatJet>& jets, int sys){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    FatJet this_jet = jets.at(i);

    this_jet *= this_jet.EnShift(sys);

    out.push_back( this_jet );
  }

  return out;

}
std::vector<FatJet> AnalyzerCore::SmearFatJets(const std::vector<FatJet>& jets, int sys){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    FatJet this_jet = jets.at(i);

    this_jet *= this_jet.ResShift(sys);

    out.push_back( this_jet );
  }

  return out;

}
//Fatjet SDMass systematics (https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging#2016%20scale%20factors%20and%20correctio)
std::vector<FatJet> AnalyzerCore::ScaleSDMassFatJets(const std::vector<FatJet>& jets, int sys){
  
  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    FatJet this_jet = jets.at(i);
    double current_SDMass = this_jet.SDMass() * (1. + double(sys) * 0.0094 );
    this_jet.SetSDMass( current_SDMass );
    
    out.push_back( this_jet );
  }
  
  return out;
  
}
std::vector<FatJet> AnalyzerCore::SmearSDMassFatJets(const std::vector<FatJet>& jets, int sys){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    //==== jets is a const vector. So in this function, we have to copy the elements like below
    FatJet this_jet = jets.at(i);
    double current_SDMass = this_jet.SDMass() * (1. + double(sys) * 0.20 );
    this_jet.SetSDMass( current_SDMass );

    out.push_back( this_jet );
  }

  return out;

}

bool AnalyzerCore::PassMETFilter(){

  //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Moriond_2018

  if(!Flag_goodVertices) return false;
  if(!IsFastSim){
    if(!Flag_globalSuperTightHalo2016Filter) return false;
  }
  if(!Flag_HBHENoiseFilter) return false;
  if(!Flag_HBHENoiseIsoFilter) return false;
  if(!Flag_EcalDeadCellTriggerPrimitiveFilter) return false;
  if(!Flag_BadPFMuonFilter) return false;
  //if(!Flag_BadChargedCandidateFilter) return false; // TODO 19/05/04 twiki says this is under review, and not recommended
  if(IsDATA && !Flag_eeBadScFilter) return false;

  if(DataYear>=2017){
    //if(!Flag_ecalBadCalibReducedMINIAODFilter) return false; // TODO WIP for UL
  }

  return true;

}

void AnalyzerCore::initializeAnalyzerTools(){

  //==== MCCorrection
  mcCorr->SetMCSample(MCSample);
  mcCorr->SetEra(GetEra());
  mcCorr->SetIsDATA(IsDATA);
  mcCorr->SetEventInfo(run, lumi, event);
  mcCorr->SetIsFastSim(IsFastSim);
  if(!IsDATA){
    mcCorr->ReadHistograms();
    mcCorr->SetupJetTagging();
  }

  puppiCorr->SetEra(GetEra());
  puppiCorr->ReadHistograms();

  //==== FakeBackgroundEstimator
  fakeEst->SetEra(GetEra());
  fakeEst->ReadHistograms();

  //==== CFBackgroundEstimator
  cfEst->SetEra(GetEra());
  cfEst->ReadHistograms();

}

double AnalyzerCore::GetPrefireWeight(int sys){

  if(IsDATA) return 1.;
  else{

    if(DataYear>2017) return 1.;
    else{

      if(sys==0) return L1PrefireReweight_Central;
      else if(sys>0) return L1PrefireReweight_Up;
      else return L1PrefireReweight_Down;

      //return mcCorr->GetPrefireWeight(photons, jets, sys);

    }

  }

  cout << "[AnalyzerCore::GetPrefireWeight] wtf" << endl;
  exit(EXIT_FAILURE);
  return 1.;

}

double AnalyzerCore::GetPileUpWeight(int N_pileup, int syst){
  if(IsDATA) return 1.;
  else return mcCorr->GetPileUpWeight(N_pileup, syst);
}

double AnalyzerCore::GetPDFWeight(LHAPDF::PDF* pdf_){

  double pdf_1 = pdf_->xfxQ(genWeight_id1, genWeight_X1, genWeight_Q);
  double pdf_2 = pdf_->xfxQ(genWeight_id2, genWeight_X2, genWeight_Q);

  return pdf_1 * pdf_2;

}

double AnalyzerCore::GetPDFReweight(){

  return GetPDFWeight(pdfReweight->NewPDF)/GetPDFWeight(pdfReweight->ProdPDF);

}

double AnalyzerCore::GetPDFReweight(int member){

  return GetPDFWeight(pdfReweight->PDFErrorSet.at(member))/GetPDFWeight(pdfReweight->ProdPDF);

}

bool AnalyzerCore::IsOnZ(double m, double width){
  if( fabs(m-M_Z) < width ) return true;
  else return false;
}

double AnalyzerCore::MT(TLorentzVector a, TLorentzVector b){
  double dphi = a.DeltaPhi(b);
  return TMath::Sqrt( 2.*a.Pt()*b.Pt()*(1.- TMath::Cos(dphi) ) );
}

double AnalyzerCore::MT2(TLorentzVector a, TLorentzVector b, Particle METv, double METgap){

  TLorentzVector METa, METb;
  METa.SetPxPyPzE( 0., 0., 0., 0.);
  double MTa, MTb;
  double tempMETa =0., tempMT2 = TMath::Max(MT(a, METv), MT(b, METv));

  while(tempMETa < METv.Pt()){

    METa.SetPxPyPzE(tempMETa*TMath::Cos(METv.Phi()), tempMETa*TMath::Sin(METv.Phi()), 0., tempMETa);
    METb = METv - METa;

    MTa = MT(METa, a);
    MTb = MT(METb, b);

    tempMT2 = TMath::Min(tempMT2, TMath::Max(MTa, MTb));

    tempMETa = tempMETa + METgap;
  }  

  return tempMT2;

}

double AnalyzerCore::projectedMET(TLorentzVector a, TLorentzVector b, Particle METv){

  if( fabs(a.DeltaPhi(METv)) < fabs(b.DeltaPhi(METv)) ){
    return (METv.Pt() * TMath::Sin(fabs(a.DeltaPhi(METv))) );
  }
  else return (METv.Pt() * TMath::Sin(fabs(b.DeltaPhi(METv))) );

}

bool AnalyzerCore::HasFlag(TString flag){

  //cout << "[AnalyzerCore::HasFlag] Userflags.size() = " << Userflags.size() << endl;
  //for(unsigned int i=0; i<Userflags.size(); i++){
  //  cout << "[AnalyzerCore::HasFlag] " << Userflags.at(i) << endl;
  //}
  
  return std::find(Userflags.begin(), Userflags.end(), flag) != Userflags.end();

}

std::vector<Muon> AnalyzerCore::MuonWithoutGap(const std::vector<Muon>& muons){

  std::vector<Muon> out;
  for(unsigned int i=0; i<muons.size(); i++){
    double this_eta = fabs( muons.at(i).Eta() );
    if( 1.444 <= this_eta && this_eta < 1.566 ) continue;

    out.push_back( muons.at(i) );
  }

  return out;

}

std::vector<Muon> AnalyzerCore::MuonPromptOnly(const std::vector<Muon>& muons, const std::vector<Gen>& gens){

  if(IsDATA) return muons;

  std::vector<Muon> out;

  for(unsigned int i=0; i<muons.size(); i++){
    if(GetLeptonType(muons.at(i), gens)<=0) continue;
    out.push_back( muons.at(i) );
  }

  return out;

}

std::vector<Muon> AnalyzerCore::MuonUsePtCone(const std::vector<Muon>& muons){

  std::vector<Muon> out;

  for(unsigned int i=0; i<muons.size(); i++){
    //==== muons is a const vector. So in this function, we have to copy the elements like below
    Muon this_muon = muons.at(i);
    this_muon.SetPtEtaPhiM( muons.at(i).PtCone(), muons.at(i).Eta(), muons.at(i).Phi(), muons.at(i).M() );
    out.push_back( this_muon );
  }

  return out;

}

Muon AnalyzerCore::MuonUsePtCone(const Muon& muon){

  //==== muon is a const object. So in this function, we have to copy the object like below
  Muon this_muon = muon;
  this_muon.SetPtEtaPhiM( muon.PtCone(), muon.Eta(), muon.Phi(), muon.M() );
  return this_muon;

}

Particle AnalyzerCore::UpdateMET(const Particle& METv, const std::vector<Muon>& muons){

  float met_x = METv.Px();
  float met_y = METv.Py();

  double px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  for(unsigned int i=0; i<muons.size(); i++){

    px_orig+= muons.at(i).MiniAODPt()*TMath::Cos(muons.at(i).Phi());
    py_orig+= muons.at(i).MiniAODPt()*TMath::Sin(muons.at(i).Phi());

    px_corrected += muons.at(i).Px();
    py_corrected += muons.at(i).Py();

  }

  met_x = met_x + px_orig - px_corrected;
  met_y = met_y + py_orig - py_corrected;

  Particle METout;
  METout.SetPxPyPzE(met_x,met_y,0,sqrt(met_x*met_x+met_y*met_y));
  return METout;

}

std::vector<Muon> AnalyzerCore::MuonApplyPtCut(const std::vector<Muon>& muons, double ptcut){

  std::vector<Muon> out;

  for(unsigned int i=0; i<muons.size(); i++){
    if(!( muons.at(i).Pt() > ptcut )) continue;
    out.push_back( muons.at(i) );
  }

  return out;

}

std::vector<Electron> AnalyzerCore::ElectronPromptOnly(const std::vector<Electron>& electrons, const std::vector<Gen>& gens){

  if(IsDATA) return electrons;

  std::vector<Electron> out;

  for(unsigned int i=0; i<electrons.size(); i++){
    if(GetLeptonType(electrons.at(i), gens)<=0) continue;
    out.push_back( electrons.at(i) );
  }

  return out;

}

std::vector<Electron> AnalyzerCore::ElectronUsePtCone(const std::vector<Electron>& electrons){

  std::vector<Electron> out;

  for(unsigned int i=0; i<electrons.size(); i++){
    //==== electrons is a const vector. So in this function, we have to copy the elements like below
    Electron this_electron = electrons.at(i);
    this_electron.SetPtEtaPhiM( electrons.at(i).PtCone(), electrons.at(i).Eta(), electrons.at(i).Phi(), electrons.at(i).M() );
    out.push_back( this_electron );
  }

  return out;

}

Electron AnalyzerCore::ElectronUsePtCone(const Electron& electron){

  //==== electron is a const object. So in this function, we have to copy the object like below
  Electron this_electron = electron;
  this_electron.SetPtEtaPhiM( electron.PtCone(), electron.Eta(), electron.Phi(), electron.M() );
  return this_electron;

}

std::vector<Electron> AnalyzerCore::ElectronApplyPtCut(const std::vector<Electron>& electrons, double ptcut){

  std::vector<Electron> out;

  for(unsigned int i=0; i<electrons.size(); i++){
    if(!( electrons.at(i).Pt() > ptcut )) continue;
    out.push_back( electrons.at(i) );
  }

  return out;

}

std::vector<Jet> AnalyzerCore::JetsAwayFromFatJet(const std::vector<Jet>& jets, const std::vector<FatJet>& fatjets, double mindr){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){

    bool Overlap = false;
    for(unsigned int j=0; j<fatjets.size(); j++){
      if( ( jets.at(i) ).DeltaR( fatjets.at(j) ) < mindr ){
        Overlap = true;
        break;
      }
    }
    if(!Overlap) out.push_back( jets.at(i) );

  }

  return out;

}

std::vector<Jet> AnalyzerCore::JetsVetoLeptonInside(const std::vector<Jet>& jets, const std::vector<Electron>& els, const std::vector<Muon>& mus, double dR){

  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){

    bool HasLeptonInside = false;

    for(unsigned int j=0; j<els.size(); j++){
      if( jets.at(i).DeltaR( els.at(j) ) < dR ){
        HasLeptonInside = true;
        break;
      }
    }
    if(HasLeptonInside) continue;

    for(unsigned int j=0; j<mus.size(); j++){
      if( jets.at(i).DeltaR( mus.at(j) ) < dR ){
        HasLeptonInside = true;
        break;
      }
    }
    if(HasLeptonInside) continue;

    //==== if all fine,
    out.push_back( jets.at(i) );

  }
  return out;

}

std::vector<FatJet> AnalyzerCore::FatJetsVetoLeptonInside(const std::vector<FatJet>& jets, const std::vector<Electron>& els, const std::vector<Muon>& mus, double dR){

  std::vector<FatJet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    FatJet this_jet = jets.at(i);

    bool HasLeptonInside = false;

    for(unsigned int j=0; j<els.size(); j++){
      if( this_jet.DeltaR( els.at(j) ) < dR ){
        HasLeptonInside = true;
        break;
      }
    }
    if(HasLeptonInside) continue;

    for(unsigned int j=0; j<mus.size(); j++){
      if( this_jet.DeltaR( mus.at(j) ) < dR ){
        HasLeptonInside = true;
        break;
      }
    }
    if(HasLeptonInside) continue;

    //==== if all fine,
    out.push_back( this_jet );

  }

  return out;

}

std::vector<Jet> AnalyzerCore::JetsAwayFromPhoton(const std::vector<Jet>& jets, const std::vector<Photon>& photons, double mindr){
  
  std::vector<Jet> out;
  for(unsigned int i=0; i<jets.size(); i++){
    
    bool Overlap = false;
    for(unsigned int j=0; j<photons.size(); j++){
      if( ( jets.at(i) ).DeltaR( photons.at(j) ) < mindr ){
        Overlap = true;
        break;
      }
    }
    if(!Overlap) out.push_back( jets.at(i) );

  }

  return out;

}


Particle AnalyzerCore::AddFatJetAndLepton(const FatJet& fatjet, const Lepton& lep){

  if(fatjet.DeltaR( lep )<0.8){
    return fatjet;
  }
  else{
    return fatjet+lep;
  }

}

//=========================================================
//==== Gen Matching Tools

void AnalyzerCore::PrintGen(const std::vector<Gen>& gens){

  cout << "===========================================================" << endl;
  cout << "RunNumber:EventNumber = " << run << ":" << event << endl;
  cout << "index\tPID\tStatus\tMIdx\tMPID\tStart\tPt\tEta\tPhi\tM" << endl;
  for(unsigned int i=2; i<gens.size(); i++){
    Gen gen = gens.at(i);
    vector<int> history = TrackGenSelfHistory(gen, gens);
    cout << i << "\t" << gen.PID() << "\t" << gen.Status() << "\t" << gen.MotherIndex() << "\t" << gens.at(gen.MotherIndex()).PID()<< "\t" << history[0] << "\t";
    printf("%.2f\t%.2f\t%.2f\t%.2f\n",gen.Pt(), gen.Eta(), gen.Phi(), gen.M());
  }

}

Gen AnalyzerCore::GetGenMatchedLepton(const Lepton& lep, const std::vector<Gen>& gens){

  //==== find status 1 lepton

  int reco_PID = -999;
  if(lep.LeptonFlavour()==Lepton::ELECTRON) reco_PID = 11;
  else if(lep.LeptonFlavour()==Lepton::MUON) reco_PID = 13;
  else{
    cout << "[AnalyzerCore::GetGenMatchedLepton] input lepton flavour not set" << endl;
    exit(EXIT_FAILURE);
  }

  double min_dR = 0.1;
  Gen gen_closest;
  for(unsigned int i=2; i<gens.size(); i++){

    Gen gen = gens.at(i);

    //==== Status 1
    if( gen.Status() != 1 ) continue;
    //==== PID
    if( abs( gen.PID() ) != reco_PID ) continue;
    //==== reject ISR?
    if( gen.MotherIndex() < 0 ) continue;
    //==== dR matching
    if( gen.DeltaR( lep ) < min_dR ){
      min_dR = gen.DeltaR( lep ) ;
      gen_closest = gen;
    }

  }

  return gen_closest;

}

Gen AnalyzerCore::GetGenMatchedPhoton(const Lepton& lep, const std::vector<Gen>& gens){
//Find if there is a photon candidate for a source of external conversion (similar pt, direction of lepton)

  double min_dR = 0.2;//1)
  Gen gen_closest;
  double pt_min = 10., dPtRelmax=0.5, dRmax=0.2;//1)
  for(unsigned int i=2; i<gens.size(); i++){

    Gen gen = gens.at(i);

    if( gen.MotherIndex() < 0 ) continue;
    if( ! ( abs(gen.PID())==22 && (gen.Status()==1 || gen.Status()==23) ) ) continue;//2)
    if( gen.Pt() < pt_min ) continue;
    if( !(lep.Pt()/gen.Pt()>(1.-dPtRelmax) && lep.Pt()/gen.Pt()<(1.+dPtRelmax)) ) continue;
    if( gen.DeltaR( lep ) > dRmax ) continue;
    if( gen.Status()==23 && !IsFinalPhotonSt23_Public(gens) ) continue;//2)

    if( gen.DeltaR( lep ) < min_dR ){
      min_dR = gen.DeltaR( lep ) ;
      gen_closest = gen;
    }

  }

  return gen_closest;

//1) Validity of cuts are checked for electron with PT>10. Note that the previous dPtRel<0.2 cut in the CatNtuple analysis was optimized for PT(e)>25.
//   External conversion rate proportionate to M^{-2}, thus muon external conversion rate is negligible.
//2) As the algorithm targets stable photon near lepton, only status-1 photon should be investigated theoretically. 
//   However, in some events, last status in photon history is 23. This is believed to be due to skimming of GEN history between pythia & MiniAOD. 
//   This case is also covered here. Note that IsFinalPhotonSt23 has to analyze whole history, i.e. computationally heavier part in this function.
//   Therefore it is put at last part of criteria for efficient functioning.
//For more details about the optimization, ask J.B..
}

vector<int> AnalyzerCore::TrackGenSelfHistory(const Gen& me, const std::vector<Gen>& gens){
//returns {index of the first history of the gen, index of the last history of the gen's mother}

  int myindex = me.Index();

  if(myindex<2){
    vector<int> out = {myindex, -1};
    return out;
  }

  int mypid = gens.at(myindex).PID();

  int currentidx = myindex;
  int motherindex = me.MotherIndex();

  while(gens.at(motherindex).PID() == mypid){

    //==== Go one generation up
    currentidx = motherindex;
    motherindex = gens.at(motherindex).MotherIndex();

    if(motherindex<0) break;
  }

  vector<int> out = {currentidx, motherindex};

  return out;

}

bool AnalyzerCore::IsFromHadron(const Gen& me, const std::vector<Gen>& gens){

  bool out = false;

  int myindex = me.Index();
  if(myindex<2) return true;

  vector<int> my_history = TrackGenSelfHistory(me, gens);
  Gen          Start = gens.at( my_history[0] );
  Gen MotherOf_Start = gens.at( my_history[1] );

  //==== Status 21~29 are from hardprocess
  //==== Means it's lepton from hardprocess
  //==== e.g., leptons from Z start their lives with status 23
  if( 20 < Start.Status() && Start.Status() < 30 ) return false;

  Gen current_me = Start; // me will always be Start
  Gen current_mother = Start; // initializing
  while( current_mother.Index() >= 2 ){

    vector<int> current_history = TrackGenSelfHistory(current_me, gens);

    //==== Go one generation up
    //==== not being used after this line
    //==== not a bug!!
    current_me = gens.at(current_history[1]);

    //==== Now look at mother of previous "me"
    current_mother = gens.at(current_history[1]);

    vector<int> current_mother_history = TrackGenSelfHistory(current_mother, gens);
    Gen StartOf_current_mother = gens.at(current_mother_history[0]);
    int current_mother_PID = current_mother.PID();

    //==== From Z,W,H,t
    if( current_mother_PID==23 || current_mother_PID==24 || current_mother_PID==25 || current_mother_PID==6 ){
      out = false;
      break;
    }
    //==== From Signal
    if( IsSignalPID(current_mother_PID) ){
      out = false;
    }
    if( (current_mother_PID==11 || current_mother_PID==13 || current_mother_PID==15 || current_mother_PID==22) && (StartOf_current_mother.Status()>20 && StartOf_current_mother.Status()<30)){
      out = false;
      break;
    }
    if( current_mother_PID>50 ){
      out=true;
      break;
    }
    if( (current_mother_PID>=1 && current_mother_PID<=5) || current_mother_PID==21 ){
      out=true; 
      break;
    }

  }

  return out;

}


int AnalyzerCore::GetPrElType_InSameSCRange_Public(int TruthIdx, const std::vector<Gen>& TruthColl){
//Abbreviation: Get Prompt Electron Type In Same Supercluster Range (Public (shared) version of J.B.'s original gen-matching code)
//Prompt e>eee (int. conv.) case, collimated electrons can be merged in one SC & track is selected among them, and reconstructed as single electron. 
//In this case, still there will be a nearby prompt electron of LepType 1-3 within supercluster merging range.
//Supercluster merging range: |dphi|<0.3/0.3 (EB/EE), |deta|<0/0.07 (EB/EE) (additionally put 0.03 for marginal difference from size of crystals)
//Ref:JINST 10 (2015) P06005, arXiv:1502.02701
//Return Value: 1/2/3:LeptonType-1/2/3 electron found in same SC range
//             -1: No LeptonType-1/2/3 electron found in same SC range
//
//Note: Hadronic origins are not included in this algo. as I am not sure of effect of many nearby particles in jets on the ele-reco performance.
//      + it is not of my interest yet.

  if(TruthIdx<2) return false;
  if(abs(TruthColl.at(TruthIdx).PID())!=11) return false;
  if(TruthColl.at(TruthIdx).Status()!=1) return false;

  float dPhiMax=0.3, dEtaMax=0.1;
  int NearbyElType=-1;

  for(unsigned int it_gen=2; it_gen<TruthColl.size(); it_gen++){
    if(TruthColl.at(it_gen).Status()!=1) continue;
    if(abs(TruthColl.at(it_gen).PID())!=11) continue;
    if(fabs(TruthColl.at(TruthIdx).Eta()-TruthColl.at(it_gen).Eta())>dEtaMax) continue;
    if(TruthColl.at(TruthIdx).DeltaPhi(TruthColl.at(it_gen))>dPhiMax) continue;

    int LepType = GetLeptonType_Public(it_gen, TruthColl);
    if(!(LepType>=1 && LepType<=3)) continue;
    NearbyElType  = LepType;
  }

  return NearbyElType;
}


int AnalyzerCore::GetLeptonType_Public(int TruthIdx, const std::vector<Gen>& TruthColl){
//Type : 1:EW-Prompt / 2: BSM-Prompt / 3:EW/BSM-Prompt-Tau Daughter 
//       4:Internal Conversion from Soft QED Radiation (PS-level) / 5:Internal Conversion from Hard Process Photon (ME-level)
//      -1:Unmatched & not EW Conversion candidate (mis-reco. or external conversion within jets)
//      -2:Hadron Daughter / -3:Daughter of Tau from Hadron or Parton / -4:Internal Conversion Daughter having hadronic origin
//      -5:External Conversion from Hard process photon / -6:External conversion from t/EWV/EWlep
//      (-4:Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
//       0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic origin or external conversion
//Note: There is no physical distinction between type 4 vs. 5 and type -5 vs. -6. The distinction is only intended for straightforward debugging.


  //Only consider Status 1 lepton
  if(TruthIdx<2) return 0;
  if(TruthColl.at(TruthIdx).Status()!=1) return 0;
  if( !(fabs(TruthColl.at(TruthIdx).PID())==11 || fabs(TruthColl.at(TruthIdx).PID())==13) ) return 0;


  int LeptonType=0;
  int MPID=0, GrMPID=0;
  vector<int> my_history, mom_history, grmom_history;
  int LastSelfIdx=-1, MotherIdx=-1, LastSelfMIdx=-1, GrMotherIdx=-1, LastSelfGrMIdx=-1;
  int Status_orig=0, MStatus_orig=0, MStatus_last=0, GrMStatus_orig=0, GrMStatus_last=0;
  bool HadronicOrigin = false;

  my_history      = TrackGenSelfHistory(TruthColl.at(TruthIdx), TruthColl);
  LastSelfIdx     = my_history[0];
  MotherIdx       = my_history[1];
  Status_orig     = TruthColl.at(LastSelfIdx).Status();
  HadronicOrigin  = IsFromHadron(TruthColl.at(TruthIdx), TruthColl);
  
  if(   MotherIdx!=-1   ){ mom_history  = TrackGenSelfHistory(TruthColl.at(MotherIdx), TruthColl);
                           LastSelfMIdx = mom_history[0];
                           GrMotherIdx  = mom_history[1];
                           MPID         = TruthColl.at(MotherIdx).PID();
                           MStatus_orig = TruthColl.at(LastSelfMIdx).Status();
                           MStatus_last = TruthColl.at(MotherIdx).Status();
                         }
  if(  GrMotherIdx!=-1  ){ grmom_history  = TrackGenSelfHistory(TruthColl.at(GrMotherIdx), TruthColl);
                           LastSelfGrMIdx = grmom_history[0];
                           GrMPID         = TruthColl.at(GrMotherIdx).PID();
                           GrMStatus_orig = TruthColl.at(LastSelfGrMIdx).Status();
                           GrMStatus_last = TruthColl.at(GrMotherIdx).Status();
                         }

  if     ( TruthIdx==-1 )                                       LeptonType= 0;
  else if( fabs(MPID)==23 || fabs(MPID)==24 || fabs(MPID)==25 ) LeptonType= 1;
  else if( IsSignalPID(MPID) )                                  LeptonType= 2;
  else if( Status_orig>20 && Status_orig<30 )                   LeptonType= 1;//1)
  else if( fabs(MPID)>50 )                                      LeptonType=-2;
  else if( fabs(MPID)==15 && MStatus_last==2 ){
           if     ( fabs(GrMPID)==23 || fabs(GrMPID)==24 || fabs(GrMPID)==25 ) LeptonType= 3;
           else if( IsSignalPID(GrMPID) )                                      LeptonType= 3;
           else if( MStatus_orig>20  && MStatus_orig<30  )                     LeptonType= 3;//1)
           else if( HadronicOrigin )                                           LeptonType=-3;//2-a)
           else if( fabs(GrMPID)==22  && GrMStatus_orig>20 && GrMStatus_orig<30 )                     LeptonType= 5;//2-b)
           else if( fabs(GrMPID)==22 )                                                                LeptonType= 4;//2-c)
           else if( (fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15) && GrMStatus_last!=2 ) LeptonType= 4;//2-d)
           else                                                                                       LeptonType= 0;
         }
  else if( fabs(MPID)==22 ){
           if( MStatus_orig>20 && MStatus_orig<30 )                            LeptonType= 5;//3-a)
           else if( HadronicOrigin )                                           LeptonType=-4;//3-b)
           else if( fabs(GrMPID)==24 || fabs(GrMPID)==23 || fabs(GrMPID)==6  ) LeptonType= 4;//3-c)
           else if( fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15 ) LeptonType= 4;//3-d)
           else                                                                LeptonType= 0;
         }
  else if( (fabs(MPID)==11 || fabs(MPID)==13 || fabs(MPID)==15) && MStatus_last!=2 && !HadronicOrigin ) LeptonType= 4;//4-a)
  else if( ((fabs(MPID)>=1 && fabs(MPID)<=5) || fabs(MPID)==21) && MStatus_last!=2 )                    LeptonType=-4;//4-b)
  else if( fabs(MPID)==6 ) LeptonType=4;//4-c)
  else LeptonType=0;


  return LeptonType;

//**footnote
//These are based on observation in DY,ZG,TT sample(DY,ZG:amcnlo+pythia, TT:powheg+pythia) for other PS generator, convention may differ.
//1) In MG5aMC@NLO generator, output of ME level generation does not have specific guage field mother if it is not on-shell specific generation.
//   e.g. u u~ > l+ l- -> fabs(MID)=1 .
//   This line is put after checking prompt mother PID checking, in order to distinguish type 1 and 2.
//2-a) e.g. a)Had > ta+X, ta>l+2nu b) q>ta+X in jet fragmentation (ta is not from hard process (ME), as it was already checked in prev. step)
//2-b) e.g. gamma>ta(+)+ta(-)+X, ta>lnu (St=2)
//2-c) e.g. " " " " " " " " " " " " " " " " " ", but soft gamma case. this is not observed in test sample but put here just in case.
//          (hadronic origin was already checked before, so no need to worry about contamination from hadronic origins at this line)
//2-d) e.g. l>tata..+l.. , ta>l+2nu (tau pair produced from internal conversion), photon radiation is not decay: GenStatus!=2
//3-a) e.g. hard process gamma>ll
//3-b) e.g. a) Had>gamma+X, gamma>ll+X (in PS+Had stage intermediate process is omitted you see just Had>Nphoton+Mhadrons+..)
//          b) q>gamma+q, gamma>ll+X in jet fragmentation or radiations of tops.
//          c) gluon>Ngamma+Mhadrons in jet fragmentation (actually observed in samples)
//3-c) e.g. W+>W+ gamma, or t>t+gamma, gamma>ll+X, not yet observed in test sample but possible (upto radiation is observed so far)
//3-d) e.g. ta>ta+gamma, gamma>ll+X, tau not from hadron(e.g. pp>tata)
//4-a) e.g. EW lep l, l>lll... just implicit conversion. 
//4-b) e.g. q or g> Nlepton +MHadrons... in parton shower history
//4-c) e.g. t>t+ll.. implicit conversion
}



int AnalyzerCore::GetLeptonType(const Lepton& lep, const std::vector<Gen>& gens){

  //==== [Type]
  //====  1 : EWPrompt
  //====  2 : Signal Daughter
  //====  3 : EWtau daughter
  //====  4 : Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit)
  //====  5 : Internal Conversion daughter from HardScatterPhoton
  //==== -1 : Unmatched & not EW Conversion candidate
  //==== -2 : Hadron daughter
  //==== -3 : Daughter of tau from hadron or parton
  //==== -4 : Internal conversion daughter(implicit,explicit) having hadronic origin 
  //==== -5 : External conversion candidate(Hard scattered photon)
  //==== -6 : External conversion from t/EWV/EWlep
  //==== (-4: Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
  //==== -------------------------------------------------------------
  //====  0 : Error
  //==== >0 : Non-fake; Non-hadronic origin
  //==== <0 : Fakes   ; Hadronic origin or external conversion

  Gen gen_closest = GetGenMatchedLepton(lep, gens); // if gen_closest.IsEmpty(), we should then look for a near photon
  //cout << "[AnalyzerCore::GetLeptonType] [Reco] pt = " << lep.Pt() << "\t, eta = " << lep.Eta() << endl;
  //cout << "[AnalyzerCore::GetLeptonType] [Gen] Index = " << gen_closest.Index() << endl;

/*
  //==== 1) Using Simple Flags
  //==== Prompt
  if( gen_closest.isPromptFinalState() ){
    return 1;
  }
  //=== From Tau
  else if( gen_closest.isDirectPromptTauDecayProductFinalState() ){
    return 2;
  }
  else{
    //=== Fake
    return -1;
  }
*/  

  //==== 2) Use Gen history

  int LeptonType=0, MatchedTruthIdx=-1;
  //==== No matched gen lepton
  if( gen_closest.IsEmpty() ){
    Gen gen_photon_closest = GetGenMatchedPhoton(lep, gens);
    int NearPhotonType = GetGenPhotonType(gen_photon_closest,gens);
    if     ( NearPhotonType<=0 ) LeptonType=-1;//1)
    else if( NearPhotonType==1 ) LeptonType=-5;//2)
    else if( NearPhotonType==2 ) LeptonType=-6;//3)
  }
  //==== Has macthed gen lepton
  else{
    MatchedTruthIdx = gen_closest.Index();
    LeptonType = GetLeptonType_Public(MatchedTruthIdx, gens);//4)
    if(LeptonType>=4 && LeptonType<=5){//5)
      int NearbyPrElType = GetPrElType_InSameSCRange_Public(MatchedTruthIdx, gens);
      if(NearbyPrElType>0) LeptonType = NearbyPrElType;
    }
  }

  return LeptonType;

//**footnote
//1) matched to no gen-lepton nor gen-photon -> mis-reco.
//   matched to no gen-lepton, but to photon with hadronic origin -> mis-reco. (e.g. pions->e) or external conversion from photon in jets
//2) matched to no gen-lepton, but to photon with non-hadronic origin (hard process) -> external conversion from photon with non-hadronic source (ME-level)
//3) matched to no gen-lepton, but to photon with non-hadronic origin (soft QED radiation) -> external conversion from photon with non-hadronic source (PS-level)
//4) matched to gen-lepton, categorize based on the truth categorization algo.: AnalyzerCore::GetLeptonType_Public(int TruthIdx, std::vector<Gen>& TruthColl)
//5) collimated e/gm objects are merged in SC, hence if there is prompt electron within SC-merging range, reco-electron's properties actually represent pre-QED-FSR prompt-electron,
//   rather than the closest internal conversion electron. Therefore shift the type to the prompt lepton's type.
//- Note: distinction between type 4 vs. 5 and type -5 vs. -6 is unphysical. it is only for debugging.
}

int AnalyzerCore::GetGenPhotonType(const Gen& genph, const std::vector<Gen>& gens){

  //==== [Type]
  //====  0 : Invalid input or Error or HardScatter is input when hardscatter is not final state
  //====  1 : HardScatter / 2: Else prompt daughter(l,V,t)
  //==== -1 : Reserved for unmatched(Not used now) / -2: Hadronic origin

  if(genph.IsEmpty()) return 0;

  int genph_index = genph.Index();

  if( genph_index<2 ) return 0;
  if( !(gens.at(genph_index).PID()==22 && (gens.at(genph_index).Status()==1 || gens.at(genph_index).Status()==23)) ) return 0;

  if(gens.at(genph_index).Status()==23){
    if(IsFinalPhotonSt23_Public(gens)) return 1;
    else                               return 0;
  }//From this pt, only St1 Photon is treated.

  vector<int> phhist = TrackGenSelfHistory(genph, gens);

  Gen          Start = gens.at(phhist[0]);
  Gen MotherOf_Start = gens.at(phhist[1]);
  int MotherOf_Start_PID = abs(MotherOf_Start.PID()); // |PID|

  bool fromhadron = IsFromHadron(genph, gens);

  if     ( 20 < Start.Status() && Start.Status() < 30 ) return 1;//1)
  else if( MotherOf_Start_PID==23 || MotherOf_Start_PID==25) return 1;//2)
  else if( MotherOf_Start_PID==24 || MotherOf_Start_PID==6 || IsSignalPID(MotherOf_Start_PID) ) return 2;//3)
  else if( fromhadron ) return -2;//4)
  else if( MotherOf_Start_PID==11 || MotherOf_Start_PID==13 || MotherOf_Start_PID==15 ) return 2;//5)
  else return 0;

//*footnote
//1) Hard process photon's original status is 20's, unless the history is skimmed in MiniAOD.
//2) If part of intermediate history is skimmed at MinoAOD, then last photon status is 1 without any preceding history lines.
//   In that case, intermediate mother is written in history, and this line catches the case.
//3) top and charged bosons can radiate photons.
// - Note that distinction between 1 and 2 is not physical. You should use both 1 & 2 for prompt photons.
//4) this category does not include top, but photons from hadrons and quarks. Predominantly, in most of the cases they are daughter of pi0.
//   Rarer, but other mesons as eta, B, or even some quarks can also radiate energetic photons.
//5) Photons radiated from lepton FSR. Sometimes they are quite energetic.
}


bool AnalyzerCore::IsFinalPhotonSt23_Public(const std::vector<Gen>& TruthColl){
//In Some XG proc events, it is observed that some of photons' last status is 23. Presumably due to skimming of generator history between pythia and MiniAOD.
//The function returns if this is the case.
//And this is designed only for 1 hard photon case as W+G or Z+G or TT+G

  bool IsFinalGammaStatus23 = false;
  bool HasStatus23Photon    = false;
  for(unsigned int i=2; i<TruthColl.size(); i++){
    int fpid  = fabs(TruthColl.at(i).PID());
    int GenSt = TruthColl.at(i).Status();
    int MPID_direct= TruthColl.at(TruthColl.at(i).MotherIndex()).PID();
    if( !((fpid!=22 && MPID_direct==22) || (fpid==22 && (GenSt==23||GenSt==1))) ) continue;

    vector<int> my_history, mom_history;
    my_history  = TrackGenSelfHistory(TruthColl.at(i), TruthColl);
    int LastSelfIdx     = my_history[0];
    int MotherIdx       = my_history[1];
    int LastSelfSt      = TruthColl.at(LastSelfIdx).Status();
    int LastSelfMIdx=-1, MStatus_orig=-1;
    if(MotherIdx!=-1){
      mom_history  = TrackGenSelfHistory(TruthColl.at(MotherIdx), TruthColl);
      LastSelfMIdx = mom_history[0];
      MStatus_orig = TruthColl.at(LastSelfMIdx).Status();
    }

    if(fpid==22){
      if(GenSt==23) {HasStatus23Photon=true; IsFinalGammaStatus23=true;}
      else if(GenSt==1 && LastSelfSt==23) {IsFinalGammaStatus23=false; break;}//a)
    }
    else if( MPID_direct==22 && MStatus_orig==23 ){ IsFinalGammaStatus23=false; break;}//b)
  }

  if(!HasStatus23Photon) return false;

  return IsFinalGammaStatus23;

//**footnotes
//a) Status-23 photon's last is 1. Thus status-23 photon is not the last history.
//b) Daughter particle of status-23 photon is found. Thus status-23 photon is not the last history.
}


bool AnalyzerCore::IsSignalPID(int pid){

  pid = abs(pid);

  //==== HeavyNeutrino
  if(pid>=9900000) return true;
  //==== ChargedHiggs
  if(pid==36 || pid==37) return true;

  return false;

}

//==== END Gen Matching Tools
//==============================================================

TH1D* AnalyzerCore::GetHist1D(TString histname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist_TH1D.find(histname);
  if(mapit != maphist_TH1D.end()) return mapit->second;

  return h;

}

TH2D* AnalyzerCore::GetHist2D(TString histname){

  TH2D *h = NULL;
  std::map<TString, TH2D*>::iterator mapit = maphist_TH2D.find(histname);
  if(mapit != maphist_TH2D.end()) return mapit->second;

  return h;

}

TH3D* AnalyzerCore::GetHist3D(TString histname){
  
  TH3D *h = NULL;
  std::map<TString, TH3D*>::iterator mapit = maphist_TH3D.find(histname);
  if(mapit != maphist_TH3D.end()) return mapit->second;
  
  return h;
  
}


void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double *xbins){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, xbins);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
          double value_x, double value_y, double value_z,
          double weight,
          int n_binx, double x_min, double x_max,
          int n_biny, double y_min, double y_max,
          int n_binz, double z_min, double z_max){
  
  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }
  
  this_hist->Fill(value_x, value_y, value_z, weight);
  
}

void AnalyzerCore::FillHist(TString histname,
          double value_x, double value_y, double value_z,
          double weight,
          int n_binx, double *xbins,
          int n_biny, double *ybins,
          int n_binz, double *zbins){
  
  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }
  
  this_hist->Fill(value_x, value_y, value_z, weight);
  
}

TH1D* AnalyzerCore::JSGetHist1D(TString suffix, TString histname){

  TH1D *h = NULL;

  std::map< TString, std::map<TString, TH1D*> >::iterator mapit = JSmaphist_TH1D.find(suffix);
  if(mapit==JSmaphist_TH1D.end()){
    return h;
  }
  else{

    std::map<TString, TH1D*> this_maphist = mapit->second;
    std::map<TString, TH1D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = JSGetHist1D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    (JSmaphist_TH1D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value, weight);

}

TH2D* AnalyzerCore::JSGetHist2D(TString suffix, TString histname){

  TH2D *h = NULL;

  std::map< TString, std::map<TString, TH2D*> >::iterator mapit = JSmaphist_TH2D.find(suffix);
  if(mapit==JSmaphist_TH2D.end()){
    return h;
  }
  else{

    std::map<TString, TH2D*> this_maphist = mapit->second;
    std::map<TString, TH2D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double x_min, double x_max,
                  int n_biny, double y_min, double y_max){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double *xbins,
                  int n_biny, double *ybins){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::WriteHist(){

  outfile->cd();
  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }

  outfile->cd();
  for(std::map< TString, std::map<TString, TH1D*> >::iterator mapit=JSmaphist_TH1D.begin(); mapit!=JSmaphist_TH1D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH1D* > this_maphist = mapit->second;


    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH1D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }

  for(std::map< TString, std::map<TString, TH2D*> >::iterator mapit=JSmaphist_TH2D.begin(); mapit!=JSmaphist_TH2D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH2D* > this_maphist = mapit->second;

    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH2D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }

}


void AnalyzerCore::FillLeptonPlots(std::vector<Lepton *> leps, TString this_region, double weight){

  for(unsigned int i=0; i<leps.size(); i++){

    TString this_itoa = TString::Itoa(i,10);

    Lepton *lep = leps[i];

    FillHist(this_region+"/Lepton_"+this_itoa+"_Pt_"+this_region, lep->Pt(), weight, 1000, 0., 1000.);
    FillHist(this_region+"/Lepton_"+this_itoa+"_Eta_"+this_region, lep->Eta(), weight, 60, -3., 3.);
    FillHist(this_region+"/Lepton_"+this_itoa+"_RelIso_"+this_region, lep->RelIso(), weight, 100, 0., 1.);
    FillHist(this_region+"/Lepton_"+this_itoa+"_MiniRelIso_"+this_region, lep->MiniRelIso(), weight, 100, 0., 1.);

    FillHist(this_region+"/Lepton_"+this_itoa+"_dXY_"+this_region, fabs(lep->dXY()), weight, 500, 0., 0.05);
    FillHist(this_region+"/Lepton_"+this_itoa+"_dXYSig_"+this_region, fabs(lep->dXY()/lep->dXYerr()), weight, 100, 0., 10);
    FillHist(this_region+"/Lepton_"+this_itoa+"_dZ_"+this_region, fabs(lep->dZ()), weight, 500, 0., 0.5);
    FillHist(this_region+"/Lepton_"+this_itoa+"_dZSig_"+this_region, fabs(lep->dZ()/lep->dZerr()), weight, 100, 0., 10);
    FillHist(this_region+"/Lepton_"+this_itoa+"_IP3D_"+this_region, fabs(lep->IP3D()), weight, 500, 0., 0.5);
    FillHist(this_region+"/Lepton_"+this_itoa+"_IP3DSig_"+this_region, fabs(lep->IP3D()/lep->IP3Derr()), weight, 100, 0., 10);

    if(lep->LeptonFlavour()==Lepton::ELECTRON){
      Electron *el = (Electron *)lep;
      FillHist(this_region+"/Lepton_"+this_itoa+"_MVANoIso_"+this_region, el->MVANoIso(), weight, 200, -1., 1.);
    }
    else if(lep->LeptonFlavour()==Lepton::MUON){
      Muon *mu = (Muon *)lep;
      FillHist(this_region+"/Lepton_"+this_itoa+"_Chi2_"+this_region, mu->Chi2(), weight, 500, 0., 50.);
      FillHist(this_region+"/Lepton_"+this_itoa+"_TrkRelIso_"+this_region, mu->TrkIso()/mu->TuneP4().Pt(), weight, 100, 0., 1.);
    }
    else{
      cout << "[AnalyzerCore::FillLeptonPlots] lepton flavour wrong.." << endl;
      exit(EXIT_FAILURE);
    }


  }

}

void AnalyzerCore::FillJetPlots(std::vector<Jet> jets, std::vector<FatJet> fatjets, TString this_region, double weight){

  for(unsigned int i=0; i<jets.size(); i++){

    TString this_itoa = TString::Itoa(i,10);
    FillHist(this_region+"/Jet_"+this_itoa+"_Pt_"+this_region, jets.at(i).Pt(), weight, 1000, 0., 1000.);
    FillHist(this_region+"/Jet_"+this_itoa+"_Eta_"+this_region, jets.at(i).Eta(), weight, 60, -3., 3.);

  }

  for(unsigned int i=0; i<fatjets.size(); i++){

    TString this_itoa = TString::Itoa(i,10);
    FillHist(this_region+"/FatJet_"+this_itoa+"_Pt_"+this_region, fatjets.at(i).Pt(), weight, 1000, 0., 1000.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_Eta_"+this_region, fatjets.at(i).Eta(), weight, 60, -3., 3.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_Mass_"+this_region, fatjets.at(i).M(), weight, 3000, 0., 3000.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_SDMass_"+this_region, fatjets.at(i).SDMass(), weight, 3000, 0., 3000.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_LSF_"+this_region, fatjets.at(i).LSF(), weight, 100, 0., 1.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_PuppiTau21_"+this_region, fatjets.at(i).PuppiTau2()/fatjets.at(i).PuppiTau1(), weight, 100, 0., 1.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_PuppiTau31_"+this_region, fatjets.at(i).PuppiTau3()/fatjets.at(i).PuppiTau1(), weight, 100, 0., 1.);
    FillHist(this_region+"/FatJet_"+this_itoa+"_PuppiTau32_"+this_region, fatjets.at(i).PuppiTau3()/fatjets.at(i).PuppiTau2(), weight, 100, 0., 1.);
  }

}

