// efficiency and l1 matching plots for aleko's SMP meeting
Plotter plot(int type=0,TString era="2016preVFP"){
  Plotter p;
  samples.clear();
  p.ScanFiles("data/Run2UltraLegacy_v1/"+era+"/ID/Electron/");
  if(type==0){
    p.AddEntry("MediumID_QPlus_v3_4");
    p.AddEntry("MediumID_QPlus_v3_4");
  }else if(type==1){
    if(era=="2017"||era=="2018"){
      p.AddEntry("Ele32_MediumID_QPlus_v3_4");
      p.AddEntry("Ele32_MediumID_QPlus_v3_4");      
    }else{
      p.AddEntry("Ele27_MediumID_QPlus_v3_4");
      p.AddEntry("Ele27_MediumID_QPlus_v3_4");
    }
  }else if(type==2){
    p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
    p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
  }else if(type==3){
    p.AddEntry("Ele12Leg2_MediumID_QPlus_v3_4");
    p.AddEntry("Ele12Leg2_MediumID_QPlus_v3_4");
  }
  p.systematics["altSig"]=p.MakeSystematic("altSig",Systematic::Type::ENVELOPE,(1<<Sample::Type::DATA),"_stat->_altSig");
  p.systematics["altBkg"]=p.MakeSystematic("altBkg",Systematic::Type::ENVELOPE,(1<<Sample::Type::DATA),"_stat->_altBkg");
  p.systematics["altTag"]=p.MakeSystematic("altTag",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_stat->_altTag");
  p.systematics["altMC"]=p.MakeSystematic("altMC",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_stat->_altMC");
  p.systematics["total"]=p.MakeSystematic("total",Systematic::Type::MULTI,0,"altSig altBkg altTag altMC");
  p.entries[0].title="data";
  p.entries[0].type=Sample::Type::DATA;
  p.entries[0].style.drawoption="e";
  p.entries[0].style.markerstyle=24;
  p.entries[1].title="MadGraph";
  p.entries[1].type=Sample::Type::SIGNAL;
  p.entries[1].replace["Data"]="MC";
  p.entries[1].style.markerstyle=24;
  TLatex l;
  l.SetTextSize(0.035);
  l.SetNDC();
  TCanvas* c=NULL;
  c=p.DrawPlot("EGamma_EffData2D_stat","project:y Xmin:0.0 Xmax:0.4 logx preliminary xtitle:'p_{T}(e^{+}) [GeV]' 1:ytitle:Efficiency 2:ytitle:'data / MC' era:"+era+" BRleg");
  l.DrawLatex(0.7,0.6,"0 #leq #eta(e) < 0.4");
  if(type==0){
    l.DrawLatex(0.25,0.4,"Medium ID");
  }else if(type==1){
    l.DrawLatex(0.25,0.4,"Single Electron Trigger");
  }else if(type==2){
    l.DrawLatex(0.25,0.4,"Double Electron Trigger Leg1");
  }else if(type==3){
    l.DrawLatex(0.25,0.4,"Double Electron Trigger Leg2");
  }
  c->SaveAs("fig/210507/"+era+"_pt"+Form("%d",type)+".png");
  c=p.DrawPlot("EGamma_EffData2D_stat","project:x Ymin:40 Ymax:45 preliminary xtitle:#eta(e^{+}) 1:ytitle:Efficiency 2:ytitle:'data / MC' era:"+era+" BRleg");
  l.DrawLatex(0.7,0.6,"40 #leq p_{T}(e) < 45");
  if(type==0){
    l.DrawLatex(0.25,0.4,"Medium ID");
  }else if(type==1){
    l.DrawLatex(0.25,0.4,"Single Electron Trigger");
  }else if(type==2){
    l.DrawLatex(0.25,0.4,"Double Electron Trigger Leg1");
  }else if(type==3){
    l.DrawLatex(0.25,0.4,"Double Electron Trigger Leg2");
  }
  c->SaveAs("fig/210507/"+era+"_eta"+Form("%d",type)+".png");
  return p;
}
Plotter plot_mi(int type=0,TString era="2016preVFP"){
  Plotter p;
  samples.clear();
  p.ScanFiles("data/Run2UltraLegacy_v1/"+era+"/ID/Electron/");
  if(type==0){
    p.AddEntry("MediumID_QPlus_v3_4");
    p.AddEntry("MediumID_QPlus_v3_4");
    p.AddFile("mi0","/data6/Users/hsseo/egm_tnp_analysis/ElectronEff_v3_5/"+era+"/MediumID_QPlus_v3_5.root");
    p.AddEntry("mi0");
    p.AddEntry("mi0");      
  }else if(type==1){
    p.AddEntry("Ele27_MediumID_QPlus_v3_4");
    p.AddEntry("Ele27_MediumID_QPlus_v3_4");
    p.AddFile("mi1","/data6/Users/hsseo/egm_tnp_analysis/ElectronEff_v3_5/"+era+"/Ele27_MediumID_QPlus_v3_5.root");
    p.AddEntry("mi1");
    p.AddEntry("mi1");      
  }else if(type==2){
    p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
    p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
    p.AddFile("mi2","/data6/Users/hsseo/egm_tnp_analysis/ElectronEff_v3_5/"+era+"/Ele23Leg1_MediumID_QPlus_v3_5.root");
    p.AddEntry("mi2");
    p.AddEntry("mi2");      
  }else if(type==3){
    p.AddEntry("Ele12Leg2_MediumID_QPlus_v3_4");
    p.AddEntry("Ele12Leg2_MediumID_QPlus_v3_4");
    p.AddFile("mi3","/data6/Users/hsseo/egm_tnp_analysis/ElectronEff_v3_5/"+era+"/Ele12Leg2_MediumID_QPlus_v3_5.root");
    p.AddEntry("mi3");
    p.AddEntry("mi3");      
  }
  p.systematics["altSig"]=p.MakeSystematic("altSig",Systematic::Type::ENVELOPE,(1<<Sample::Type::DATA),"_stat->_altSig");
  p.systematics["altBkg"]=p.MakeSystematic("altBkg",Systematic::Type::ENVELOPE,(1<<Sample::Type::DATA),"_stat->_altBkg");
  p.systematics["altTag"]=p.MakeSystematic("altTag",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_stat->_altTag");
  p.systematics["altMC"]=p.MakeSystematic("altMC",Systematic::Type::ENVELOPE,(1<<Sample::Type::SIGNAL),"_stat->_altMC");
  p.systematics["total"]=p.MakeSystematic("total",Systematic::Type::MULTI,0,"altSig altBkg altTag altMC");
  p.entries[0].title="data (MadGraph template)";
  p.entries[0].type=Sample::Type::DATA;
  p.entries[0].style.drawoption="e";
  p.entries[0].style.markerstyle=24;
  p.entries[1].title="MadGraph";
  p.entries[1].type=Sample::Type::SIGNAL;
  p.entries[1].replace["Data"]="MC";
  p.entries[1].style.markerstyle=24;
  p.entries[2].title="data (MiNNLO template)";
  p.entries[2].type=Sample::Type::DATA;
  p.entries[2].style.drawoption="e";
  p.entries[2].style.linecolor=1;
  p.entries[2].style.markercolor=1;
  p.entries[2].style.markerstyle=25;
  p.entries[3].title="MiNNLO";
  p.entries[3].type=Sample::Type::SIGNAL;
  p.entries[3].replace["Data"]="MC";
  p.entries[3].style.drawoption="e";
  p.entries[3].style.markerstyle=25;
  TLatex l;
  l.SetTextSize(0.035);
  l.SetNDC();
  TCanvas* c=NULL;
  c=p.DrawPlot("EGamma_EffData2D_stat","project:y Xmin:0.0 Xmax:0.4 logx preliminary xtitle:'p_{T}(e^{+}) [GeV]' 1:ytitle:Efficiency 2:ytitle:'Ratio' era:"+era+" BRleg");
  l.DrawLatex(0.65,0.55,"0.0 #leq #eta(e) < 0.4");
  if(type==0){
    l.DrawLatex(0.5,0.65,"Medium ID");
  }else if(type==1){
    l.DrawLatex(0.5,0.65,"Single Electron Trigger");
  }else if(type==2){
    l.DrawLatex(0.5,0.65,"Double Electron Trigger Leg1");
  }else if(type==3){
    l.DrawLatex(0.5,0.65,"Double Electron Trigger Leg2");
  }
  c->SaveAs("fig/210507/mi"+era+"_pt"+Form("%d",type)+".png");
  c=p.DrawPlot("EGamma_EffData2D_stat","project:x Ymin:40 Ymax:45 preliminary xtitle:#eta(e^{+}) 1:ytitle:Efficiency 2:ytitle:'Ratio' era:"+era+" BRleg");
  l.DrawLatex(0.65,0.55,"40 #leq p_{T}(e) < 45");
  if(type==0){
    l.DrawLatex(0.5,0.65,"Medium ID");
  }else if(type==1){
    l.DrawLatex(0.5,0.65,"Single Electron Trigger");
  }else if(type==2){
    l.DrawLatex(0.5,0.65,"Double Electron Trigger Leg1");
  }else if(type==3){
    l.DrawLatex(0.5,0.65,"Double Electron Trigger Leg2");
  }
  c->SaveAs("fig/210507/mi"+era+"_eta"+Form("%d",type)+".png");
  return p;
}
Plotter plot_l1(){
  TString era="2017";
  Plotter p;
  samples.clear();
  p.ScanFiles("data/Run2UltraLegacy_v1/"+era+"/ID/Electron/");
  p.AddEntry("egammaEffi.txt_EGM2D_Ele23Leg1_MediumID_QPlus_v3_2");
  p.AddEntry("egammaEffi.txt_EGM2D_Ele23Leg1_MediumID_QPlus_v3_2");
  p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
  p.AddEntry("Ele23Leg1_MediumID_QPlus_v3_4");
  p.entries[0].title="data";
  p.entries[0].type=Sample::Type::DATA;
  p.entries[0].style.drawoption="e";
  p.entries[0].style.markerstyle=24;
  p.entries[1].title="MadGraph";
  p.entries[1].type=Sample::Type::SIGNAL;
  p.entries[1].replace["Data"]="MC";
  p.entries[1].style.markerstyle=24;
  p.entries[2].title="data (L1 p_{T}-threshold)";
  p.entries[2].type=Sample::Type::DATA;
  p.entries[2].style.drawoption="e";
  p.entries[2].style.linecolor=1;
  p.entries[2].style.markercolor=1;
  p.entries[2].style.markerstyle=25;
  p.entries[3].title="MadGraph (L1 p_{T}-threshold)";
  p.entries[3].type=Sample::Type::SIGNAL;
  p.entries[3].replace["Data"]="MC";
  p.entries[3].style.drawoption="e";
  p.entries[3].style.markerstyle=25;
  TCanvas* c=p.DrawPlot("EGamma_EffData2D","project:y Xmin:0.0 Xmax:0.4 logx preliminary xtitle:'p_{T}(e^{+}) [GeV]' 1:ytitle:Efficiency 2:ytitle:'others / data' era:"+era+" BRleg");
  TLatex l;
  l.SetTextSize(0.035);
  l.SetNDC();
  l.DrawLatex(0.6,0.6,"0.0 #leq #eta(e) < 0.4");
  l.DrawLatex(0.45,0.7,"Double Electron Trigger Leg1");
  c->SaveAs("fig/210507/l1.png");
  return p;
}
void plot_all(){
  plot_mi(0,"2016postVFP");
  plot_mi(1,"2016postVFP");
  plot_mi(2,"2016postVFP");
  plot_mi(3,"2016postVFP");
  return;
  plot(0,"2018");
  plot(1,"2018");
  plot(2,"2018");
  plot(3,"2018");
  plot_l1();
}
