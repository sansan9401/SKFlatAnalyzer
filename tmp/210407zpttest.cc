TString Sample="DYJetsToMuMu_MiNNLO";
//TString Sample="DYJets";
//TString Sample="DYJets_MG";

void test1(){
  Plotter p;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);
  
  vector<double> cuts={100,400,900,1600,2500,3600,4900,8100,10000};
  TGraphErrors *g=new TGraphErrors;
  g->SetName("g");
  g->SetTitle("g");
  for(int i=0;i<cuts.size();i++){
    double cut=cuts.at(i);
    TGraphErrors *gg=new TGraphErrors;
    gg->SetName(Form("gg%d",i));
    gg->SetTitle(Form("gg%d",i));
    for(int iy=0;iy<25;iy++){
      TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap_nozptweight",Form("project:x Ymin:%f Ymax:%f xmin:0 xmax:%f",0.2*iy,0.2*(iy+1),cut));
      gg->SetPoint(iy,0.2*iy+0.1,h->GetMean());
      gg->SetPointError(iy,0,h->GetMeanError());
    }
    TF1* f=new TF1(Form("f%d",i),"[0]*(1-[1]*cosh(x))");
    gg->Fit(f);
    new TCanvas;
    gg->Draw();
    g->SetPoint(i,cut,f->GetParameter(1));
    g->SetPointError(i,0,f->GetParError(1));
  }
  new TCanvas;
  g->Draw();
}
TGraphErrors* Hist2Graph(TH1* h,bool switchaxis=false){
  TGraphErrors* g=NULL;
  if(!h) return g;
  g=new TGraphErrors;
  g->SetName(h->GetName()+TString("_graph"));
  g->SetTitle(h->GetTitle()+TString("_graph"));
 
  for(int i=1;i<h->GetNcells()-1;i++){
    if(switchaxis){
      g->SetPoint(i-1,h->GetBinContent(i),h->GetBinCenter(i));
      g->SetPointError(i-1,h->GetBinError(i),0);
    }else{
      g->SetPoint(i-1,h->GetBinCenter(i),h->GetBinContent(i));
      g->SetPointError(i-1,0,h->GetBinError(i));
    }
  }
  return g;
}
TGraphErrors* Divide(TGraphErrors* g1,TGraphErrors* g2){
  TGraphErrors* g=NULL;
  if(!g1) return g;
  if(!g2) return g;
  g=new TGraphErrors;
  for(int i=0;i<g2->GetN();i++){
    double x=g2->GetX()[i];
    if(x==0) continue;
    g->SetPoint(i,x,g1->Eval(x)/g2->Eval(x));
  }
  return g;
}

double MyFunc(double* xx,double* par){
  double x=xx[0];
  double s=par[0];
  if(s<=0) s=1e-6;
  x/=s;
  vector<double> k={10,50,100,200,500,1000,2000,5000,1e100};
  int nregion=k.size();
  vector<double> a(nregion,0);
  vector<double> b(nregion,0);
  vector<double> c(nregion,0);
   
  a[0]=par[1];
  b[0]=par[2];
  c[0]=par[3];
  if(x<k[0]) return exp(a[0]*pow(x,b[0])+c[0]);
  for(int i=1;i<nregion;i++){
    a[i]=par[3+i];
    if(i==1){
      b[i]=a[i-1]*b[i-1]*pow(k[i-1],b[i-1]-1)/a[i];
      c[i]=a[i-1]*pow(k[i-1],b[i-1])+c[i-1]-a[i];
    }else{
      b[i]=a[i-1]*b[i-1]*pow(k[i-1]-k[i-2]+1,b[i-1]-1)/a[i];
      c[i]=a[i-1]*pow(k[i-1]-k[i-2]+1,b[i-1])+c[i-1]-a[i];
    }
    if(x<k[i]) return exp(a[i]*pow(x-k[i-1]+1,b[i])+c[i]);
  }
  return 0;
}    
double MyFunc2(double* xx,double* par){
  vector<double> x2={pow(xx[0],2)};
  return par[0]*xx[0]*MyFunc(&x2[0],&par[1]);
}
    
void test2(){
  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  TGraphErrors* gfinal=new TGraphErrors;
  TGraphErrors* g0=NULL;
  TCanvas* c=new TCanvas;
  for(int iy=0;iy<25;iy++){
    TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap_nozptweight",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",0.2*iy,0.2*(iy+1),400.));
    h->SetName(Form("y%d",iy));
    h->SetTitle(Form("y%d",iy));
    double scale=1/h->GetBinContent(1);
    h->Scale(scale);
    TGraphErrors* g=Hist2Graph(h,true);
    if(iy==0){
      h->Draw();
      g0=g;
      gfinal->SetPoint(iy,0,1);
    }
    else{
      h->Draw("same");
      new TCanvas;
      TGraphErrors* gg0=Divide(g,g0);
      gg0->SetTitle(Form("ratio%d",iy));
      gg0->Draw();
      gg0->GetHistogram()->GetYaxis()->SetRangeUser(0,2);
      gg0->Fit("pol0","","",0.2,0.8);
      gfinal->SetPoint(iy,0.2*iy+0.1,gg0->GetFunction("pol0")->GetParameter(0));
      gfinal->SetPointError(iy,0,gg0->GetFunction("pol0")->GetParError(0));
      c->cd();      
    }
  }
  new TCanvas;
  TF1* f=new TF1("f","(1-[0]*cosh(x))/(1-[0])");
  gfinal->Fit(f);
  gfinal->Draw();
}

void test3(){
  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  TGraphErrors* gfinal=new TGraphErrors;
  TGraphErrors* g0=NULL;
  TCanvas* c=new TCanvas;
  for(int iy=0;iy<25;iy++){
    TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap_nozptweight",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",0.2*iy,0.2*(iy+1),400.));
    //TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",0.2*iy,0.2*(iy+1),400.));
    h->SetName(Form("y%d",iy));
    h->SetTitle(Form("y%d",iy));
    h->Fit("expo","0","",0,10);
    double scale=1/h->GetFunction("expo")->Eval(0);
    h->Scale(scale);
    TGraphErrors* g=Hist2Graph(h,true);
    if(iy==0){
      h->Draw();
      g0=g;
      gfinal->SetPoint(iy,0,1);
    }
    else{
      h->Draw("same");
      new TCanvas;
      TGraphErrors* gg0=Divide(g,g0);
      gg0->SetTitle(Form("ratio%d",iy));
      gg0->Draw();
      gg0->GetHistogram()->GetYaxis()->SetRangeUser(0,2);
      gg0->Fit("pol0","","",0.1,0.8);
      gfinal->SetPoint(iy,0.2*iy+0.1,gg0->GetFunction("pol0")->GetParameter(0));
      gfinal->SetPointError(iy,0,gg0->GetFunction("pol0")->GetParError(0));
      c->cd();      
    }
  }
  new TCanvas;
  //TF1* f=new TF1("f","(1-[0]*cosh(x))/(1-[0])");
  TF1* f=new TF1("f","1-[0]*x*x");
  gfinal->Fit(f);
  //gfinal->Fit(f,"","",0,2.4);
  gfinal->Draw();
}

TGraphErrors* test_y(TString Sample_="",TString suffix="_nozptweight"){
  if(Sample_!="") Sample=Sample_;

  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap"+suffix,"project:x widthweight Ymin:-0.2 Ymax:0.2 xmin:0 xmax:10000");
  h->SetName("a");
  h->Scale(1/h->GetBinContent(1));

  int npar=12;
  TF1* f=new TF1("f",MyFunc,0,10000,npar);
  f->FixParameter(0,1);
  f->SetParameter(1,-0.03);
  f->SetParameter(2,1);
  f->SetParameter(3,0.007);
  f->SetParameter(4,-0.04);
  f->SetParameter(5,-0.02);
  f->SetParameter(6,-0.01);
  for(int i=7;i<npar;i++){
    f->SetParameter(i,-0.001);
    f->SetParLimits(i,-0.1,0);
  }
  
  h->Fit(f,"","",0,10);
  h->Fit(f,"","",0,100);
  h->Fit(f,"","",0,200);
  h->Fit(f,"","",0,500);
  h->Fit(f,"","",0,1000);
  h->Fit(f,"","",0,2000);
  h->Fit(f,"","",0,5000);
  h->Fit(f,"","",0,10000);
  f->Draw("same");


  TGraphErrors* g=new TGraphErrors;
  for(int i=0;i<25;i++){
    TH1* hh=p.GetHist(0,"mm2016a/gen_diptdirap"+suffix,Form("project:x widthweight absy Ymin:%f Ymax:%f xmin:0 xmax:%f",0.2*i,0.2*(i+1),100.));
    hh->SetStats(0);
    hh->SetTitle(Form("%.1f<|y|<%.1f",0.2*i,0.2*(i+1)));
    hh->GetXaxis()->SetTitle("p_{T}");
    TCanvas *c=new TCanvas;
    TF1* ff=new TF1(Form("ff%d",i),MyFunc2,0,100,npar+1);
    ff->SetNpx(10000);
    ff->SetParameter(0,1);
    ff->SetParameter(1,1);
    for(int i=2;i<npar+1;i++){
      ff->FixParameter(i,f->GetParameter(i-1));
    }
    //hh->Fit(ff,"","",0,20);
    hh->Fit(ff);
    ff->Draw("same");
    if(i==0){
      g->SetPoint(i,0,ff->GetParameter(1));
      g->SetPointError(i,0,ff->GetParError(1));
    }else{
      g->SetPoint(i,0.2*i+0.1,ff->GetParameter(1));
      g->SetPointError(i,0,ff->GetParError(1));
    }
  }
  TCanvas *c=new TCanvas;
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x/[1]))/(1-[0])"); fff->SetParameter(1,1);
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x))/(1-[0])");
  TF1* fff=new TF1("fff","1-[0]*x*x");
  g->Fit(fff);
  //g->Fit(fff,"","",0,2.4);
  g->Draw();
  g->GetHistogram()->GetXaxis()->SetTitle("|y|");
  g->GetHistogram()->GetYaxis()->SetTitle("p_{T}^{2} scale");
  TFile fout("output.root","update");
  if(!fout.Get(Sample)) g->Write(Sample);
  return g;
}

TGraphErrors* test_m(double ptmax=100){
  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap_nozptweight","project:x widthweight xmin:0 xmax:10000");
  h->SetName("a");
  h->Scale(1/h->GetBinContent(1));

  int npar=12;
  TF1* f=new TF1("f",MyFunc,0,10000,npar);
  f->FixParameter(0,1);
  f->SetParameter(1,-0.03);
  f->SetParameter(2,1);
  f->SetParameter(3,0.007);
  f->SetParameter(4,-0.04);
  f->SetParameter(5,-0.02);
  f->SetParameter(6,-0.01);
  for(int i=7;i<npar;i++){
    f->SetParameter(i,-0.001);
    f->SetParLimits(i,-0.1,0);
  }
  
  TCanvas* cc=new TCanvas;

  h->Fit(f,"","",0,10);
  h->Fit(f,"","",0,100);
  h->Fit(f,"","",0,200);
  h->Fit(f,"","",0,500);
  h->Fit(f,"","",0,1000);
  h->Fit(f,"","",0,2000);
  h->Fit(f,"","",0,5000);
  h->Fit(f,"","",0,10000);
  f->Draw("same");


  TGraphErrors* g=new TGraphErrors;
  vector<double> massbin={50,60,70,80,100,120,150,200,300,400,1000};
  for(int i=0;i<massbin.size()-1;i++){
    TH1* hh=p.GetHist(0,"mm2016a/gen_diptdimass_nozptweight",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",massbin[i],massbin[i+1],ptmax));
    hh->SetStats(0);
    hh->SetTitle(Form("%.0f<m<%.0f",massbin[i],massbin[i+1]));
    hh->GetXaxis()->SetTitle("p_{T}");
    TCanvas *c=new TCanvas;
    TF1* ff=new TF1(Form("ff%d",i),MyFunc2,0,100,npar+1);
    ff->SetNpx(10000);
    ff->SetParameter(0,1);
    ff->SetParameter(1,1);
    for(int i=2;i<npar+1;i++){
      ff->FixParameter(i,f->GetParameter(i-1));
    }
    //hh->Fit(ff,"","",0,20);
    hh->Fit(ff);
    ff->Draw("same");
    g->SetPoint(i,(massbin[i]+massbin[i+1])*0.5,ff->GetParameter(1));
    g->SetPointError(i,0,ff->GetParError(1));
  }
  TCanvas *c=new TCanvas;
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x/[1]))/(1-[0])"); fff->SetParameter(1,1);
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x))/(1-[0])");
  //TF1* fff=new TF1("fff","1-[0]*x*x");
  //g->Fit(fff);
  //g->Fit(fff,"","",0,2.4);
  g->Draw();
  g->GetHistogram()->GetXaxis()->SetTitle("m(Z)");
  g->GetHistogram()->GetYaxis()->SetTitle("p_{T}^{2} scale");
  return g;
}
TGraphErrors* test_isr(double ptmax=100){
  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  TGraphErrors* g=new TGraphErrors;
  vector<double> massbin={50,60,70,80,100,120,150,200,300,400,1000};
  TH1* h=p.GetHist(0,"mm2016a/gen_diptdimass_nozptweight",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",80.,100.,ptmax));
  for(int i=0;i<massbin.size()-1;i++){
    TH1* hh=p.GetHist(0,"mm2016a/gen_diptdimass_nozptweight",Form("project:x widthweight Ymin:%f Ymax:%f xmin:0 xmax:%f",massbin[i],massbin[i+1],ptmax));
    g->SetPoint(i,(massbin[i]+massbin[i+1])*0.5,hh->GetMean()/h->GetMean());
    g->SetPointError(i,0,hh->GetMeanError()/h->GetMean());
  }
  TCanvas *c=new TCanvas;
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x/[1]))/(1-[0])"); fff->SetParameter(1,1);
  //TF1* fff=new TF1("fff","(1-[0]*cosh(x))/(1-[0])");
  //TF1* fff=new TF1("fff","1-[0]*x*x");
  //g->Fit(fff);
  //g->Fit(fff,"","",0,2.4);
  g->Draw();
  g->GetHistogram()->GetXaxis()->SetTitle("m(Z)");
  g->GetHistogram()->GetYaxis()->SetTitle("<p_{T}>/<p_{T}^{Q=90}>");
  return g;
}
void plot_test_m(){
  TCanvas *c=new TCanvas;
  vector<double> ptmax={50,80,100};
  for(int i=0;i<ptmax.size();i++){
    TGraphErrors* g=test_m(ptmax[i]);
    g->SetName(Form("ptmax%d",i));
    c->cd();
    if(i==0) g->Draw("APL");
    else g->Draw("same PL");
  }  
}
void plot_test_isr(){
  TCanvas *c=new TCanvas;
  vector<double> ptmax={50,80,100};
  for(int i=0;i<ptmax.size();i++){
    TGraphErrors* g=test_isr(ptmax[i]);
    g->SetName(Form("ptmax%d",i));
    c->cd();
    if(i==0) g->Draw("APL");
    else g->Draw("same PL");
  }
}

void plot_zptweight(){
  SMPAnalyzerCore m;
  m.MCSample="DYJets";
  m.SetEra("2016preVFP");
  m.SetupZptWeight();

  TCanvas* c=new TCanvas;

  vector<TGraph*> gs;
  for(int i=0;i<6;i++){
    TGraph* g=new TGraph;
    g->SetName(Form("g%d",i));
    gs.push_back(g);
    for(int j=0;j<100;j++){
      gs[i]->SetPoint(j,j,m.GetZptWeight(j,i*0.4+0.1,Lepton::Flavour::MUON));
    }
    gs[i]->SetLineColor(i+1);
    c->cd();
    if(i==0) gs[i]->Draw("AL");
    else gs[i]->Draw("L");
  } 
}
void plot_zpt(){
  AFBPlotter aa;
  aa.DrawPlot("mm2016a/m[80,100]/y[0,0.4]/dipt","norm widthweight absy xmax:100 varibit:4 suffix:_zpty0");
  aa.DrawPlot("mm2016a/m[80,100]/y[2.0,2.4]/dipt","norm widthweight absy xmax:100 varibit:4 suffix:_zpty0");
}
void plot_pt2(){
  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);
  
  p.DrawPlot("mm2016a/gen_diptdirap_nozptweight","project:x widthweight xmin:0 xmax:100");
  p.DrawPlot("mm2016a/gen_dipt2dirap_nozptweight","project:x widthweight xmin:0 xmax:10000");
}
void plot_rapidity(){
  Plotter p;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);
  
  p.entries[0].SetHistPrefix("gen_");
  p.entries[1].SetHistPrefix("m80to100/");
  
  p.DrawPlot("mm2016a/diptdirap_nozptweight","project:y widthweight type:1");

}
void plot_test_y(){
  TGraphErrors* gmi=test_y("DYJetsToMuMu_MiNNLO");
  TGraphErrors* gmidata=test_y("DYJetsToMuMu_MiNNLO","");
  TGraphErrors* gamc=test_y("DYJets");
  TGraphErrors* gamcdata=test_y("DYJets","");
  TGraphErrors* gmg=test_y("DYJets_MG");
  TGraphErrors* gmgdata=test_y("DYJets_MG","");
  gmi->SetName("gmi");
  gmidata->SetName("gmidata");
  gamc->SetName("gamc");
  gamcdata->SetName("gamcdata");
  gmg->SetName("gmg");
  gmgdata->SetName("gmgdata");
  vector<TGraphErrors*> graphs={gmi,gamc,gmg,gmidata,gamcdata,gmgdata};
  TCanvas *c=new TCanvas;
  for(int i=0;i<graphs.size();i++){
    TGraphErrors* graph=graphs[i];
    graph->GetFunction("fff")->SetBit(TF1::EStatusBits::kNotDraw);
    graph->SetFillColor(0);
    if(TString(graph->GetName()).Contains("data")) graph->SetMarkerStyle(24);
    else graph->SetMarkerStyle(22);
    graph->SetLineColor(i%3+1);
    graph->SetMarkerColor(i%3+1);
    if(i==0) graph->Draw("APL");
    else graph->Draw("PL");
  }
  TLegend* leg=new TLegend(0.18,0.12,0.58,0.5);
  leg->AddEntry(gmi,"MiNNLO");
  leg->AddEntry(gamc,"aMC@NLO");
  leg->AddEntry(gmg,"MadGraph");
  leg->AddEntry(gmidata,"MiNNLO reweighted to data");
  leg->AddEntry(gamcdata,"aMC@NLO reweighted to data");
  leg->AddEntry(gmgdata,"MadGraph reweighted to data");
  leg->Draw();
}

TGraphErrors* test_ratio(TString Sample_="",TString suffix="_nozptweight"){
  if(Sample_!="") Sample=Sample_;

  Plotter p;
  Verbosity=0;
  p.ScanFiles("/data6/Users/hsseo/SKFlatOutput/Run2UltraLegacy_v1/ZptWeight");
  p.AddEntry(samples["2016preVFP/ZptWeight_"+Sample]);

  new TCanvas;
  TH1* h0=p.GetHist(0,"mm2016a/gen_dipt2dirap"+suffix,"project:x widthweight Ymin:-0.2 Ymax:0.2 xmin:0 xmax:10000");
  h0->SetName("h0");
  h0->Fit("expo","0","",0,10);
  double scale=1/h0->GetFunction("expo")->Eval(0);
  h0->Scale(scale);

  int npar=12;
  TF1* f0=new TF1("f0",MyFunc,0,10000,npar);
  f0->SetNpx(10000);
  f0->FixParameter(0,1);
  f0->SetParameter(1,-0.03);
  f0->SetParameter(2,1);
  f0->SetParameter(3,0.007);
  f0->SetParameter(4,-0.04);
  f0->SetParameter(5,-0.02);
  f0->SetParameter(6,-0.01);
  for(int i=7;i<npar;i++){
    f0->SetParameter(i,-0.001);
    f0->SetParLimits(i,-0.1,0);
  }
  
  h0->Fit(f0,"","",0,10);
  h0->Fit(f0,"","",0,100);
  h0->Fit(f0,"","",0,200);
  h0->Fit(f0,"","",0,500);
  h0->Fit(f0,"","",0,1000);
  h0->Fit(f0,"","",0,2000);
  h0->Fit(f0,"","",0,5000);
  h0->Fit(f0,"","",0,10000);
  f0->SetParameter(3,f0->GetParameter(3)-log(f0->Eval(0)));
  f0->Draw("same");

  for(int i=1;i<25;i+=5){
    new TCanvas;
    TH1* h=p.GetHist(0,"mm2016a/gen_dipt2dirap"+suffix,Form("project:x widthweight Ymin:%f Ymax:%f absy xmin:0 xmax:10000",i*0.2,(i+1)*0.2));
    h->SetName(Form("h%d",i));
    h->Fit("expo","0","",0,10);
    double scale=1/h->GetFunction("expo")->Eval(0);
    h->Scale(scale);

    int npar=12;
    TF1* f=new TF1(Form("f%d",i),MyFunc,0,10000,npar); 
    f->SetNpx(10000);
    f->FixParameter(0,1);
    f->SetParameter(1,-0.03);
    f->SetParameter(2,1);
    f->SetParameter(3,0.007);
    f->SetParameter(4,-0.04);
    f->SetParameter(5,-0.02);
    f->SetParameter(6,-0.01);
    for(int i=7;i<npar;i++){
      f->SetParameter(i,-0.001);
      f->SetParLimits(i,-0.1,0);
    }
    
    h->Fit(f,"","",0,10);
    h->Fit(f,"","",0,100);
    h->Fit(f,"","",0,200);
    h->Fit(f,"","",0,500);
    h->Fit(f,"","",0,1000);
    h->Fit(f,"","",0,2000);
    h->Fit(f,"","",0,5000);
    h->Fit(f,"","",0,10000);
    f->SetParameter(3,f->GetParameter(3)-log(f->Eval(0)));
    f->Draw("same");
    
    new TCanvas;
    TGraphErrors *g=new TGraphErrors;
    for(int i=0;i<99;i++){
      double x=i+1;
      g->SetPoint(i,x*x,f->GetX(f0->Eval(x*x))/x/x);
      //g->SetPoint(i,f->Eval(x*x),f->GetX(f0->Eval(x*x))/x/x);
    }
    g->GetHistogram()->GetYaxis()->SetRangeUser(0,1.1);
    g->Draw();
  }
  return NULL;
}
