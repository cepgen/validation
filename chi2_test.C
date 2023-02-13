// clang-format off
#include "CepGenEnvironment.h"
// clang-format on

#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "SampleHandler.h"

void chi2_dataset(const SampleHandler&, TGraphErrors&, cepgen::strfun::Parameterisation*);
void chi2_binned_dataset(const SampleHandler&,
                         const std::vector<cepgen::Limits>&,
                         TGraphErrors&,
                         cepgen::strfun::Parameterisation*);

void chi2_test(int sf = 301, bool binned = false, bool logx = true) {
  const cepgen::Limits q2_lims{1.e-1, 1.e5};
  const int num_bins = 50;

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
  samples.emplace_back(SampleHandler::fromCLAS("samples/nmc_0p75-65gev2.csv")
                           // + SampleHandler::fromCLAS("samples/nmc.csv")
                           + SampleHandler::fromNMC("samples/nmc_0p8-62gev2.csv"),
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

  pdg::MCDFileParser::parse(cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mass_width_2021.mcd");
  const auto sfs = cepgen::StructureFunctionsFactory::get().build(sf);

  cepgen::ROOTCanvas c("chi2_sf_comparison", "CepGen structure functions");
  c.SetLegendX1(0.175);
  c.SetLegendY1(0.775);  // NW
  /*c.SetLegendX1(0.45);
  c.SetLegendY1(0.25);*/ // SE

  TMultiGraph mg;
  const auto bins = q2_lims.split(num_bins, logx);
  for (const auto& sample : samples) {
    auto* gr = new TGraphErrors;
    if (binned)
      chi2_binned_dataset(sample.handler, bins, *gr, sfs.get());
    else
      chi2_dataset(sample.handler, *gr, sfs.get());
    gr->SetMarkerColor(sample.colour);
    gr->SetMarkerStyle(sample.style);
    mg.Add(gr);
    c.AddLegendEntry(gr, sample.name.data(), "ep");
  }

  mg.Draw("ap");
  mg.GetHistogram()->SetTitle(Form("%s (GeV^{2})\\Normalised #chi^{2}", binned ? "Q^{2} bin" : "Q^{2}"));
  c.Prettify(mg.GetHistogram());
  if (logx) {
    mg.GetXaxis()->SetLimits(q2_lims.min(), q2_lims.max());
    c.SetLogx();
  }
  {
    mg.GetHistogram()->SetMinimum(std::max(mg.GetHistogram()->GetMinimum(), 1.e-5));
    mg.GetHistogram()->SetMaximum(std::max(mg.GetHistogram()->GetMaximum() * 1.e2, 1.e3));
    c.SetLogy();
  }
  c.SetGrid(1, 1);
  //c.GetLegend()->SetHeader(
  //    Form("#bf{%s}", cepgen::utils::split(cepgen::StructureFunctionsFactory::get().describe(sf), ' ')[0].data()));
  c.SetTopLabel(Form("#bf{%s} F_{2}^{p}",
                     cepgen::utils::split(cepgen::StructureFunctionsFactory::get().describe(sf), ' ')[0].data()));

  c.Save("pdf,png");
}

void chi2_dataset(const SampleHandler& smp, TGraphErrors& chi2_vs_q2, cepgen::strfun::Parameterisation* sfs) {
  chi2_vs_q2.Clear();
  for (const auto& q2 : smp.q2Values()) {
    const auto gr = smp.xBjGraph(q2);
    double chi2 = 0.;
    size_t num_points = 0;
    for (int i = 0; i < gr.GetN(); ++i) {
      double xbj = gr.GetX()[i], f2 = gr.GetY()[i];
      const double f2_th = sfs->F2(xbj, q2);
      chi2 += pow(f2_th - f2, 2) / f2_th;
      num_points++;
    }
    chi2_vs_q2.AddPoint(q2, num_points > 1 ? chi2 / (num_points - 1) : 0.);
  }
}

void chi2_binned_dataset(const SampleHandler& smp,
                         const std::vector<cepgen::Limits>& bins,
                         TGraphErrors& chi2_vs_q2,
                         cepgen::strfun::Parameterisation* sfs) {
  chi2_vs_q2.Clear();
  for (const auto& bin : bins) {
    double chi2 = 0.;
    size_t num_points = 0;
    /*auto gr_th = TGraphErrors(), gr_data = TGraphErrors();
    for (const auto& q2 : smp.q2Values()) {
      if (!bin.contains(q2))
        continue;
      const auto gr = smp.xBjGraph(q2);
      for (int i = 0; i < gr.GetN(); ++i) {
        double xbj = gr.GetX()[i], f2 = gr.GetY()[i];
        gr_data.AddPoint(xbj, f2);
        gr_th.AddPoint(xbj, sfs->F2(xbj, q2));
      }
    }
    chi2_vs_q2.AddPoint(bin.x(0.5), num_points > 0 ? chi2 / num_points : 0.);
    chi2_vs_q2.SetPointError(chi2_vs_q2.GetN() - 1, 0.5 * bin.range(), 0.);*/
    for (const auto& q2 : smp.q2Values()) {
      if (!bin.contains(q2))
        continue;
      const auto gr = smp.xBjGraph(q2);
      for (int i = 0; i < gr.GetN(); ++i) {
        double xbj = gr.GetX()[i], f2 = gr.GetY()[i];
        const double f2_th = sfs->F2(xbj, q2);
        chi2 += pow(f2_th - f2, 2) / f2_th;
        num_points++;
      }
    }
    chi2_vs_q2.AddPoint(bin.x(0.5), num_points > 1 ? chi2 / (num_points - 1) : 0.);
    chi2_vs_q2.SetPointError(chi2_vs_q2.GetN() - 1, 0.5 * bin.range(), 0.);
  }
}
