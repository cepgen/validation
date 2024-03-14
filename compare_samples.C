// clang-format off
#include "CepGenEnvironment.h"
// clang-format on

#include <THStack.h>

#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "CepGenAddOns/ROOTWrapper/ROOTTreeInfo.h"

void compare_samples(const vector<pair<string, string> >& samples,
                     const vector<string>& var_in,
                     const string& selection = "",
                     const string& label = "",
                     int nbins = 0,
                     double xmin = 0.,
                     double xmax = 0.,
                     bool logx = false,
                     bool logy = false,
                     bool show_legend = true,
                     bool ratio = false,
                     const string& output_file = "samples_comparison") {
  cepgen::ROOTCanvas c(output_file, label, ratio);
  //c.SetLegendX1(0.4);
  //c.SetLegendY1(0.77);
  c.SetGrid(1, 1);

  vector<TH1*> hists;
  auto* hs = c.Make<THStack>();
  size_t i = 0;
  for (const auto& smp : samples) {
    TChain chain("events");
    auto variable = var_in.at(0) + ">>htemp" + to_string(i);
    if (nbins > 0 && (xmin != 0. || xmax != 0.))
      variable += "(" + to_string(nbins) + "," + to_string(xmin) + "," + to_string(xmax) + ")";
    for (const auto& filename : cepgen::utils::split(smp.first, '+'))
      chain.Add(filename.data());
    ROOT::CepGenRun run_info;
    run_info.attach(chain.GetFile());
    const auto num_entries = chain.GetEntriesFast();
    chain.Draw(variable.data(), selection.data());
    if (auto* tmp = gPad->GetPrimitive(Form("htemp%zu", i))) {
      //auto* hist = (TH1*)tmp->Clone();
      auto* hist = cepgen::AddUnderOverflowBins((TH1D*)tmp->Clone());
      hist->Scale(run_info.xsect / run_info.num_events, "width");
      hist->SetLineColor(cepgen::ROOTCanvas::colours[i]);
      hist->SetMarkerColor(cepgen::ROOTCanvas::colours[i]);
      hist->SetLineWidth(3);
      hists.emplace_back(hist);
      //hs->Add(hist, "hist");
      hs->Add(hist, "e");
      if (show_legend)
        c.AddLegendEntry(hist, smp.second);
      ++i;
    }
  }
  hs->Draw("nostack");
  if (var_in.size() > 1) {
    string tok1 = var_in.at(1), tok2 = "";
    if (var_in.size() > 2) {  // unit specified
      tok1 += " (" + var_in.at(2) + ")";
      tok2 = " (pb/" + var_in.at(2) + ")";
    }
    hs->GetHistogram()->SetTitle((";" + tok1 + ";d#sigma/d" + var_in.at(1) + tok2).data());
  }
  c.Prettify(hs);
  if (logx)
    c.SetLogx();
  if (logy)
    c.SetLogy();
  //c.RatioPlot(hists.at(0), hists);

  c.Save("pdf,png");
}
