// clang-format off
#include "StrFunParams.h"
// clang-format on

#include <TGraph.h>
#include <TMultiGraph.h>

#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "DISsamples/DISsample.h"

void compare_datasets(bool logx = true,
                      bool logy = true,
                      bool in_w = false,
                      bool plot2d = false,
                      const string& output = "datasets_comparison",
                      const vector<string>& formats = {"pdf"}) {
  cepgen::ROOTCanvas c(output);
  c.SetLegendX1(in_w ? 0.15 : 0.4);
  c.SetLegendY1(0.77);
  c.SetGrid(1, 1);

  //--- data points

  struct DataSample {
    explicit DataSample(const DISsample& hnd, const std::string& aname, int astyle, int acolour = kBlack)
        : handler(hnd), name(aname), style(astyle), colour(acolour) {}
    DISsample handler;
    std::string name;
    int style{0};
    int colour{0};
  };
  std::vector<DataSample> samples;
  samples.emplace_back(DISsample::fromCLAS("samples/clas_0p225-4p725gev2.csv"), "CLAS", 24, kOrange);
  samples.emplace_back(DISsample::fromBCDMS("samples/bcdms.csv"), "BCDMS", 30, kCyan + 2);
  samples.emplace_back(  //DISsample::fromCLAS("samples/nmc_0p75-65gev2.csv")
                         // + DISsample::fromCLAS("samples/nmc.csv")
      DISsample::fromNMC("samples/nmc_0p8-62gev2.csv"),
      "NMC",
      25);
  samples.emplace_back(DISsample::fromCLAS("samples/e665.csv"), "E665", 27, kGreen + 2);
  samples.emplace_back(DISsample::fromZEUS("samples/zeus.csv") + DISsample::fromCLAS("samples/zeus_0p11-0p65gev2.csv") +
                           DISsample::fromCLAS("samples/zeus_1p5-15gev2.csv") +
                           DISsample::fromZEUS("samples/zeus_0p6-17gev2.csv") +
                           DISsample::fromCLAS("samples/zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26,
                       kRed + 1);
  samples.emplace_back(
      DISsample::fromCLAS("samples/h1.csv") + DISsample::fromCLAS("samples/h1_neutral_2004.csv") +
          DISsample::fromCLAS("samples/h1_4p5-1600gev2.csv") + DISsample::fromCLAS("samples/hera-h1-nc.csv") +
          DISsample::fromH11997("samples/h1_1p5-150gev2.csv") + DISsample::fromH1lowQ2("samples/h1_compton.csv") +
          DISsample::fromH1lowQ2("samples/h1_0p35-3p5gev2.csv"),
      "H1",
      22,
      kBlue - 2);
  samples.emplace_back(DISsample::fromHermes("samples/hermes.csv"), "HERMES", 28, kMagenta + 2);
  //samples.emplace_back(DISsample::fromCCFR("samples/ccfr.txt"), "CCFR", 29, kRed + 2);

  TMultiGraph mg;
  if (plot2d) {
    size_t i = 0;
    const std::string plt_type = "cont2";
    std::vector<TGraph2D*> plts;
    for (const auto& sample : samples) {
      auto* plt = new TGraph2D(in_w ? sample.handler.wQ2Values() : sample.handler.xbjQ2Values());
      plt->Draw((plt_type + (i == 0 ? "" : ",same")).data());
      plt->SetLineColorAlpha(sample.colour, 0.5);
      c.Prettify(plt->GetHistogram());
      plt->SetMinimum(0.);
      plt->SetMaximum(0.9);
      plts.emplace_back(plt);
      ++i;
    }
    if (logx) {
      if (in_w)
        plts[0]->GetHistogram()->GetXaxis()->SetLimits(1.e-1, 300.);
      else
        plts[0]->GetHistogram()->GetXaxis()->SetLimits(1.e-7, 1.);
      c.SetLogx();
    }
    if (logy) {
      plts[0]->GetHistogram()->GetYaxis()->SetLimits(1.e-2, 1.e7);
      c.SetLogy();
    }
  } else {
    for (const auto& sample : samples) {
      auto* proj = new TGraph(in_w ? sample.handler.wQ2Span() : sample.handler.xbjQ2Span());
      proj->SetMarkerStyle(sample.style);
      proj->SetMarkerSize(0.8);
      proj->SetMarkerColor(sample.colour);
      mg.Add(proj);
      if (proj->GetN() > 0)
        c.AddLegendEntry(proj, sample.name.data(), "p");
    }

    //--- general plotting

    mg.Draw("ap");
    mg.GetHistogram()->SetTitle(Form("%s\\Q^{2} (GeV^{2})", in_w ? "W (GeV)" : "x_{Bj}"));
    c.Prettify(mg.GetHistogram());
    mg.GetHistogram()->SetTitle("");

    //mg.GetHistogram()->GetXaxis()->SetLimits(xbj_lims.min(), xbj_lims.max());
    if (logx) {
      if (in_w)
        mg.GetHistogram()->GetXaxis()->SetLimits(9.e-1, 350.);
      else
        mg.GetHistogram()->GetXaxis()->SetLimits(1.e-7, 1.);
      c.SetLogx();
    }
    if (logy) {
      mg.GetHistogram()->SetMinimum(1.e-2);
      //mg.GetHistogram()->SetMaximum(TMath::Min(5., mg.GetHistogram()->GetMaximum()));
      mg.GetHistogram()->SetMaximum(1.e7);
      c.SetLogy();
    } else
      mg.GetHistogram()->GetYaxis()->SetRangeUser(1.e-6, TMath::Min(1.5, mg.GetHistogram()->GetMaximum()));
  }

  for (const auto& fmt : formats)
    c.Save(fmt);
}
