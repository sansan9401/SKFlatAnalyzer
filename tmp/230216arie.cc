#include"EfficiencyPlotter.cc"

void plot(){
  EfficiencyPlotter aa;
  aa.DrawPlot("mm2017/m80to100/leta","rebin:2 BMleg");
  aa.DrawPlot("mm2017/m80to100/leta_fine","BMleg");
  aa.DrawPlot("mu2017/m80to100/leta_fine","BMleg");
  aa.DrawPlot("mm2017/m80to100/leta30_fine","BMleg");
}

void plot_eta(TString filename,double eta1=2,double eta2=2.1){
  Plotter p;
  TFile f(filename);
  TH2D* data2d=(TH2D*)f.Get("muonEffi_data_probe_eta_probe_pt");
  TH2D* sim2d=(TH2D*)f.Get("muonEffi_mc_probe_eta_probe_pt");
  TH2D* sf2d=(TH2D*)f.Get("SF_probe_eta_probe_pt");
  TCanvas *c=new TCanvas;
  c->Divide(1,2);
  
  c->cd(1);
  gPad->SetPad(0,0.35,1,1);
  gPad->SetBottomMargin(0.02);
  gPad->SetTopMargin(c->GetTopMargin()/0.65);
  int ieta1=data2d->GetXaxis()->FindBin(eta1);
  int ieta2=data2d->GetXaxis()->FindBin(eta2);
  vector<TH1D*> hists;
  hists.push_back((TH1D*)data2d->ProjectionY(Form("data, %.1f<#eta<%.1f",eta1,eta1+0.1),ieta1,ieta1));
  hists.push_back((TH1D*)sim2d->ProjectionY(Form("sim, %.1f<#eta<%.1f",eta1,eta1+0.1),ieta1,ieta1));
  hists.push_back((TH1D*)data2d->ProjectionY(Form("data, %.1f<#eta<%.1f",eta2,eta2+0.1),ieta2,ieta2));
  hists.push_back((TH1D*)sim2d->ProjectionY(Form("sim, %.1f<#eta<%.1f",eta2,eta2+0.1),ieta2,ieta2));
  TLegend* leg=new TLegend(0.7,0.1,0.9,0.3);
  for(int i=0;i<hists.size();i++){
    hists[i]->SetOption("hist e");
    hists[i]->SetStats(0);
    hists[i]->SetDirectory(0);
    hists[i]->SetLineWidth(2);
    if(i%2==0){
      hists[i]->SetMarkerSize(0.7);
    }else{
      hists[i]->SetLineStyle(2);
    }
    if(i<2){
      hists[i]->SetLineColor(1);
    }else{
      hists[i]->SetLineColor(2);
    }
    leg->AddEntry(hists[i],hists[i]->GetName(),"lp");
    if(i==0){
      hists[i]->Draw();
      hists[i]->GetXaxis()->SetLabelSize(0);
      hists[i]->GetXaxis()->SetTitle("");
      hists[i]->GetYaxis()->SetTitle("eff");
      hists[i]->SetTitle(gSystem->BaseName(gSystem->DirName(filename)));
      hists[i]->GetYaxis()->SetRangeUser(0.71,1.09);
      double scale=1/TMath::Min(gPad->GetHNDC(),gPad->GetWNDC());
      hists[i]->SetTitleSize(0.04*scale,"XYZ");
      hists[i]->SetLabelSize(0.04*scale,"XYZ");
    }
    else hists[i]->Draw(hists[i]->GetOption()+TString(" same"));
  }
  leg->Draw();

  c->cd(2);
  gPad->SetPad(0,0,1,0.365);
  gPad->SetTopMargin(0.02);
  gPad->SetBottomMargin(c->GetBottomMargin()/0.35);
  gPad->SetGridx();gPad->SetGridy();

  for(int i=0;i<hists.size();i+=2){
    TH1* num=(TH1*)hists[i]->Clone();
    num->SetDirectory(0);
    num->Divide(hists[i+1]);
    if(i==0){
      num->Draw();
      num->SetTitle("");
      num->GetXaxis()->SetTitle("p_{T} [GeV]");
      num->GetYaxis()->SetTitle("SF");
      num->GetYaxis()->SetRangeUser(0.81,1.19);
      double scale=1/TMath::Min(gPad->GetHNDC(),gPad->GetWNDC());
      num->SetTitleSize(0.04*scale,"XYZ");
      num->SetLabelSize(0.04*scale,"XYZ");
      num->GetYaxis()->SetTitleOffset(1.8/scale);
    }
    else num->Draw(num->GetOption()+TString(" same"));
  }
  c->Draw();
}
void plot_eta_all(){
  plot_eta("/data6/Users/wonjun/egm_tnp_analysis/MuonEff_v10_0/UL2017_ID_all/result.root",0.4,0.2);
  plot_eta("/data6/Users/wonjun/egm_tnp_analysis/MuonEff_v10_0/UL2017_IsoMu24_all/result.root",0.4,0.2);
  plot_eta("/data6/Users/wonjun/egm_tnp_analysis/MuonEff_v10_0/UL2017_Mu17_all/result.root",0.4,0.2);
  plot_eta("/data6/Users/wonjun/egm_tnp_analysis/MuonEff_v10_0/UL2017_Mu8_all/result.root",0.4,0.2);
}
