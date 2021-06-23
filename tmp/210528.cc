//efficiency validiation plots with systematics
{
  TString prefixes[]={"mm2016a","mm2016b","mm2017","mm2018",
		      "mu2016a","mu2016b","mu2017","mu2018",
		      "ee2016a","ee2016b","ee2017","ee2018",
		      "el2016a","el2016b","el2017","el2018",};

  EfficiencyPlotter mu("data ^amc+tau_amc+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    mu.SavePlot("amc/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'muon p_{T}' sysname:efficiencySF");
    mu.SavePlot("amc/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta' sysname:efficiencySF");
    mu.SavePlot("amc/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu) sysname:efficiencySF");
    mu.SavePlot("amc/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(#mu#mu) 1:logy sysname:efficiencySF");
    mu.SavePlot("amc/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)' sysname:efficiencySF");
  }

  EfficiencyPlotter mumg("data ^mg+tau_mg+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    mumg.SavePlot("mg/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'muon p_{T}' sysname:efficiencySF");
    mumg.SavePlot("mg/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta' sysname:efficiencySF");
    mumg.SavePlot("mg/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu) sysname:efficiencySF");
    mumg.SavePlot("mg/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(#mu#mu) 1:logy sysname:efficiencySF");
    mumg.SavePlot("mg/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)' sysname:efficiencySF");
  }

  EfficiencyPlotter mumi("data ^minnlo+tau_minnlo+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    if(prefix.Contains("2017")||prefix.Contains("2018")) continue;
    mumi.SavePlot("minnlo/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'muon p_{T}' sysname:efficiencySF");
    mumi.SavePlot("minnlo/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta' sysname:efficiencySF");
    mumi.SavePlot("minnlo/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu) sysname:efficiencySF");
    mumi.SavePlot("minnlo/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(#mu#mu) 1:logy sysname:efficiencySF");
    mumi.SavePlot("minnlo/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)' sysname:efficiencySF");
  }

  EfficiencyPlotter el("data ^amc+tau_amc+vv+wjets+tttw+ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    el.SavePlot("amc/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'electron p_{T}' sysname:efficiencySF");
    el.SavePlot("amc/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta' sysname:efficiencySF");
    el.SavePlot("amc/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ee) sysname:efficiencySF");
    el.SavePlot("amc/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(ee) 1:logy sysname:efficiencySF");
    el.SavePlot("amc/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)' sysname:efficiencySF");

    if(prefix.Contains("ee")) continue;
    el.SavePlot("amc_tight/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt_TightID_Selective xmax:80 norm xtitle:'electron p_{T}' sysname:efficiencySF");
    el.SavePlot("amc_tight/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta_TightID_Selective norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta' sysname:efficiencySF");
    el.SavePlot("amc_tight/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt_TightID_Selective norm xmax:200 rebin:2 xtitle:p_{T}(ee) sysname:efficiencySF");
    el.SavePlot("amc_tight/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass_TightID_Selective norm xtitle:m(ee) 1:logy sysname:efficiencySF");
    el.SavePlot("amc_tight/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)' sysname:efficiencySF");
  }  

  EfficiencyPlotter elmg("data ^mg+tau_mg+vv+wjets+tttw+ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    elmg.SavePlot("mg/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'electron p_{T}' sysname:efficiencySF");
    elmg.SavePlot("mg/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta' sysname:efficiencySF");
    elmg.SavePlot("mg/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ee) sysname:efficiencySF");
    elmg.SavePlot("mg/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(ee) 1:logy sysname:efficiencySF");
    elmg.SavePlot("mg/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)' sysname:efficiencySF");
  }  

  EfficiencyPlotter elmi("data ^minnlo+tau_minnlo+vv+wjets+tttw+ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    if(prefix.Contains("2017")||prefix.Contains("2018")) continue;
    elmi.SavePlot("minnlo/"+prefix+"_lpt","histname:"+prefix+"/m52to150/lpt xmax:80 norm xtitle:'electron p_{T}' sysname:efficiencySF");
    elmi.SavePlot("minnlo/"+prefix+"_leta","histname:"+prefix+"/m52to150/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta' sysname:efficiencySF");
    elmi.SavePlot("minnlo/"+prefix+"_dipt","histname:"+prefix+"/m52to150/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ee) sysname:efficiencySF");
    elmi.SavePlot("minnlo/"+prefix+"_dimass","histname:"+prefix+"/m52to150/dimass norm xtitle:m(ee) 1:logy sysname:efficiencySF");
    elmi.SavePlot("minnlo/"+prefix+"_dirap","histname:"+prefix+"/m52to150/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)' sysname:efficiencySF");
  }  
}
