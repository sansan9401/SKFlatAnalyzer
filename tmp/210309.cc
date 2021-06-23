//efficiency validiation plots
{
  TString prefixes[]={"mm2016a","mm2016b","mm2017","mm2018",
		      "mu2016a","mu2016b","mu2017","mu2018",
		      "ee2016a","ee2016b","ee2017","ee2018",
		      "el2016a","el2016b","el2017","el2018",};

  EfficiencyPlotter mu("data ^amc+tau_amc+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    mu.SavePlot("amc/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'muon p_{T}'");
    mu.SavePlot("amc/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading muon p_{T}'");
    mu.SavePlot("amc/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading muon p_{T}'");
    mu.SavePlot("amc/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta'");
    mu.SavePlot("amc/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu)");
    mu.SavePlot("amc/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(#mu#mu) 1:logy");
    mu.SavePlot("amc/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)'");

    TString suffix="_noPUweight";
    mu.SavePlot("amc/"+prefix+"_lpt"+suffix,"histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'muon p_{T}' varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_l0pt"+suffix,"histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading muon p_{T}' varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_l1pt"+suffix,"histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading muon p_{T}' varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_leta"+suffix,"histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta' varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_dipt"+suffix,"histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu) varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_dimass"+suffix,"histname:"+prefix+"/m60to120/dimass norm xtitle:m(#mu#mu) 1:logy varibit:12 suffix:_noPUweight");
    mu.SavePlot("amc/"+prefix+"_dirap"+suffix,"histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)' varibit:12 suffix:_noPUweight");
  }

  EfficiencyPlotter mumg("data ^mg+tau_mg+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    mumg.SavePlot("mg/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'muon p_{T}'");
    mumg.SavePlot("mg/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading muon p_{T}'");
    mumg.SavePlot("mg/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading muon p_{T}'");
    mumg.SavePlot("mg/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta'");
    mumg.SavePlot("mg/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu)");
    mumg.SavePlot("mg/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(#mu#mu) 1:logy");
    mumg.SavePlot("mg/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)'");

    TString suffix="_noPUweight";
    mumg.SavePlot("mg/"+prefix+"_lpt"+suffix,"histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'muon p_{T}' varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_l0pt"+suffix,"histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading muon p_{T}' varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_l1pt"+suffix,"histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading muon p_{T}' varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_leta"+suffix,"histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta' varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_dipt"+suffix,"histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu) varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_dimass"+suffix,"histname:"+prefix+"/m60to120/dimass norm xtitle:m(#mu#mu) 1:logy varibit:12 suffix:_noPUweight");
    mumg.SavePlot("mg/"+prefix+"_dirap"+suffix,"histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)' varibit:12 suffix:_noPUweight");
  }

  EfficiencyPlotter mumi("data ^minnlo+tau_minnlo+vv+wjets+tttw+1.8*ss");
  for(auto prefix:prefixes){
    if(prefix.Contains("el")||prefix.Contains("ee")) continue;
    if(prefix.Contains("2017")||prefix.Contains("2018")) continue;
    mumi.SavePlot("minnlo/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'muon p_{T}'");
    mumi.SavePlot("minnlo/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading muon p_{T}'");
    mumi.SavePlot("minnlo/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading muon p_{T}'");
    mumi.SavePlot("minnlo/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'muon #eta'");
    mumi.SavePlot("minnlo/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(#mu#mu)");
    mumi.SavePlot("minnlo/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(#mu#mu) 1:logy");
    mumi.SavePlot("minnlo/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(#mu#mu)'");
  }

  EfficiencyPlotter el("data ^amc+tau_amc+vv+wjets+tttw");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    el.SavePlot("amc/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'electron p_{T}'");
    el.SavePlot("amc/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading electron p_{T}'");
    el.SavePlot("amc/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading electron p_{T}'");
    el.SavePlot("amc/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta'");
    el.SavePlot("amc/"+prefix+"_lsceta","histname:"+prefix+"/m60to120/lsceta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta_{SC}'");
    el.SavePlot("amc/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ee)");
    el.SavePlot("amc/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(ee) 1:logy");
    el.SavePlot("amc/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)'");

    if(prefix.Contains("ee")) continue;
    el.SavePlot("amc_tight/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt_TightID_Selective xmax:80 norm xtitle:'electron p_{T}'");
    el.SavePlot("amc_tight/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt_TightID_Selective xmax:80 norm xtitle:'leading electron p_{T}'");
    el.SavePlot("amc_tight/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt_TightID_Selective xmax:80 norm xtitle:'sub-leading electron p_{T}'");
    el.SavePlot("amc_tight/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta_TightID_Selective norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta'");
    el.SavePlot("amc_tight/"+prefix+"_lsceta","histname:"+prefix+"/m60to120/lsceta_TightID_Selective norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta_{SC}'");
    el.SavePlot("amc_tight/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt_TightID_Selective norm xmax:200 rebin:2 xtitle:p_{T}(ee)");
    el.SavePlot("amc_tight/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass_TightID_Selective norm xtitle:m(ee) 1:logy");
    el.SavePlot("amc_tight/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)'");
  }  

  EfficiencyPlotter elmi("data ^minnlo+tau_minnlo+vv+wjets+tttw");
  for(auto prefix:prefixes){
    if(prefix.Contains("mu")||prefix.Contains("mm")) continue;
    if(prefix.Contains("2017")||prefix.Contains("2018")) continue;
    elmi.SavePlot("minnlo/"+prefix+"_lpt","histname:"+prefix+"/m60to120/lpt xmax:80 norm xtitle:'electron p_{T}'");
    elmi.SavePlot("minnlo/"+prefix+"_l0pt","histname:"+prefix+"/m60to120/l0pt xmax:80 norm xtitle:'leading electron p_{T}'");
    elmi.SavePlot("minnlo/"+prefix+"_l1pt","histname:"+prefix+"/m60to120/l1pt xmax:80 norm xtitle:'sub-leading electron p_{T}'");
    elmi.SavePlot("minnlo/"+prefix+"_leta","histname:"+prefix+"/m60to120/leta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta'");
    elmi.SavePlot("minnlo/"+prefix+"_lsceta","histname:"+prefix+"/m60to120/lsceta norm xmin:-2.5 xmax:2.5 noleg xtitle:'electron #eta_{SC}'");
    elmi.SavePlot("minnlo/"+prefix+"_dipt","histname:"+prefix+"/m60to120/dipt norm xmax:200 rebin:2 xtitle:p_{T}(ee)");
    elmi.SavePlot("minnlo/"+prefix+"_dimass","histname:"+prefix+"/m60to120/dimass norm xtitle:m(ee) 1:logy");
    elmi.SavePlot("minnlo/"+prefix+"_dirap","histname:"+prefix+"/m60to120/dirap norm xmin:-2.5 xmax:2.5 noleg xtitle:'y(ee)'");
  }  
}
