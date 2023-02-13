// clang-format off
#include "commons.h"
// clang-format on

#include <TGraph.h>
#include <TMultiGraph.h>

#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "SampleHandler.h"

void compare_datasets(bool logx = true, bool logy = true, bool in_w = false, bool plot2d = false) {
  cepgen::ROOTCanvas c("datasets_comparison");
  c.SetLegendX1(in_w ? 0.15 : 0.4);
  c.SetLegendY1(0.77);
  c.SetGrid(1, 1);

  //--- data points

  struct DataSample {
    explicit DataSample(const SampleHandler& hnd, const std::string& aname, int astyle, int acolour = kBlack)
        : handler(hnd), name(aname), style(astyle), colour(acolour) {}
    SampleHandler handler;
    std::string name;
    int style{0};
    int colour{0};
  };
  std::vector<DataSample> samples;
  samples.emplace_back(SampleHandler::fromCLAS("samples/clas_0p225-4p725gev2.csv"), "CLAS", 24, kOrange);
  samples.emplace_back(SampleHandler::fromBCDMS("samples/bcdms.csv"), "BCDMS", 30, kCyan + 2);
  samples.emplace_back(  //SampleHandler::fromCLAS("samples/nmc_0p75-65gev2.csv")
                         // + SampleHandler::fromCLAS("samples/nmc.csv")
      SampleHandler::fromNMC("samples/nmc_0p8-62gev2.csv"),
      "NMC",
      25);
  samples.emplace_back(SampleHandler::fromCLAS("samples/e665.csv"), "E665", 27, kGreen + 2);
  samples.emplace_back(SampleHandler::fromZEUS("samples/zeus.csv") +
                           SampleHandler::fromCLAS("samples/zeus_0p11-0p65gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_1p5-15gev2.csv") +
                           SampleHandler::fromZEUS("samples/zeus_0p6-17gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26,
                       kRed + 1);
  samples.emplace_back(
      SampleHandler::fromCLAS("samples/h1.csv") + SampleHandler::fromCLAS("samples/h1_neutral_2004.csv") +
          SampleHandler::fromCLAS("samples/h1_4p5-1600gev2.csv") + SampleHandler::fromCLAS("samples/hera-h1-nc.csv") +
          SampleHandler::fromH11997("samples/h1_1p5-150gev2.csv") +
          SampleHandler::fromH1lowQ2("samples/h1_compton.csv") +
          SampleHandler::fromH1lowQ2("samples/h1_0p35-3p5gev2.csv"),
      "H1",
      22,
      kBlue - 2);
  samples.emplace_back(SampleHandler::fromHermes("samples/hermes.csv"), "HERMES", 28, kMagenta + 2);

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
        plts[0]->GetHistogram()->GetXaxis()->SetLimits(1., 1.e3);
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
      proj->SetMarkerSize(0.5);
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
        mg.GetHistogram()->GetXaxis()->SetLimits(1., 1.e3);
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

  c.Save("pdf");
}
