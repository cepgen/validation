// clang-format off
#include "CepGenEnvironment.h"
// clang-format on

#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "DISsamples/DISsample.h"

using namespace std;

void chi2_dataset(const DISsample& smp, TGraphErrors& chi2_vs_q2, cepgen::strfun::Parameterisation* sf_mod) {
  chi2_vs_q2.Clear();
  for (const auto& q2 : smp.q2Values()) {
    const auto gr = smp.xBjGraph(q2);
    double chi2 = 0.;
    size_t num_points = 0;
    for (int i = 0; i < gr.GetN(); ++i) {
      double xbj = gr.GetX()[i], f2 = gr.GetY()[i];
      const double f2_th = sf_mod->F2(xbj, q2);
      chi2 += pow(f2_th - f2, 2) / f2_th;
      num_points++;
    }
    chi2_vs_q2.AddPoint(q2, num_points > 1 ? chi2 / (num_points - 1) : 0.);
  }
}

void chi2_binned_dataset(const DISsample& smp,
                         const vector<cepgen::Limits>& bins,
                         TGraphErrors& chi2_vs_q2,
                         cepgen::strfun::Parameterisation* sf_mod) {
  chi2_vs_q2.Clear();
  for (const auto& bin : bins) {
    double chi2 = 0.;
    size_t num_points = 0;
    for (const auto& q2 : smp.q2Values()) {
      if (!bin.contains(q2))
        continue;
      const auto gr = smp.xBjGraph(q2);
      for (int i = 0; i < gr.GetN(); ++i) {
        double xbj = gr.GetX()[i], f2 = gr.GetY()[i];
        const double f2_th = sf_mod->F2(xbj, q2);
        chi2 += pow(f2_th - f2, 2) / f2_th;
        num_points++;
      }
    }
    chi2_vs_q2.AddPoint(bin.x(0.5), num_points > 1 ? chi2 / (num_points - 1) : 0.);
    chi2_vs_q2.SetPointError(chi2_vs_q2.GetN() - 1, 0.5 * bin.range(), 0.);
  }
}

void chi2_test(int sf = 301, bool binned = false, bool logx = true, const string& output_name = "chi2_sf_comparison") {
  cepgen::initialise();

  const auto sfs = (sf == 0) ? cepgen::StructureFunctionsFactory::get().modules() : vector<int>{sf};

  const cepgen::Limits q2_lims{1.e-1, 1.e5};
  const int num_bins = 50;

  struct DataSample {
    explicit DataSample(const DISsample& hnd, const string& aname, int astyle, int acolour = kBlack)
        : handler(hnd), name(aname), style(astyle), colour(acolour) {}
    DISsample handler;
    string name;
    int style{0};
    int colour{0};
  };
  vector<DataSample> samples;
  samples.emplace_back(DISsample::fromCLAS("clas_0p225-4p725gev2.csv"), "CLAS", 24, kOrange);
  samples.emplace_back(DISsample::fromBCDMS("bcdms.csv"), "BCDMS", 30, kCyan + 2);
  samples.emplace_back(DISsample::fromCLAS("nmc_0p75-65gev2.csv")
                           // + DISsample::fromCLAS("nmc.csv")
                           + DISsample::fromNMC("nmc_0p8-62gev2.csv"),
                       "NMC",
                       25);
  samples.emplace_back(DISsample::fromCLAS("e665.csv"), "E665", 27, kGreen + 2);
  samples.emplace_back(DISsample::fromZEUS("zeus.csv") + DISsample::fromCLAS("zeus_0p11-0p65gev2.csv") +
                           DISsample::fromCLAS("zeus_1p5-15gev2.csv") + DISsample::fromZEUS("zeus_0p6-17gev2.csv") +
                           DISsample::fromCLAS("zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26,
                       kRed + 1);
  samples.emplace_back(DISsample::fromCLAS("h1.csv") + DISsample::fromCLAS("h1_neutral_2004.csv") +
                           DISsample::fromCLAS("h1_4p5-1600gev2.csv") + DISsample::fromCLAS("hera-h1-nc.csv") +
                           DISsample::fromH11997("h1_1p5-150gev2.csv") + DISsample::fromH1lowQ2("h1_compton.csv") +
                           DISsample::fromH1lowQ2("h1_0p35-3p5gev2.csv"),
                       "H1",
                       22,
                       kBlue - 2);
  samples.emplace_back(DISsample::fromHermes("hermes.csv"), "HERMES", 28, kMagenta + 2);

  for (const auto& sf_name : sfs) {
    const auto sf_mod = cepgen::StructureFunctionsFactory::get().build(sf_name);

    cepgen::ROOTCanvas c(output_name + Form("_%d", sf_name), "CepGen structure functions");
    c.SetLegendX1(0.175);
    c.SetLegendY1(0.775);  // NW
    /*c.SetLegendX1(0.45);
  c.SetLegendY1(0.25);*/ // SE

    TMultiGraph mg;
    const auto bins = q2_lims.split(num_bins, logx);
    for (const auto& sample : samples) {
      auto* gr = new TGraphErrors;
      if (binned)
        chi2_binned_dataset(sample.handler, bins, *gr, sf_mod.get());
      else
        chi2_dataset(sample.handler, *gr, sf_mod.get());
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
      mg.GetHistogram()->SetMinimum(max(mg.GetHistogram()->GetMinimum(), 1.e-5));
      mg.GetHistogram()->SetMaximum(max(mg.GetHistogram()->GetMaximum() * 10., 1.e4));
      c.SetLogy();
    }
    c.SetGrid(1, 1);
    c.SetTopLabel(
        Form("#bf{%s} F_{2}^{p}",
             cepgen::utils::split(cepgen::StructureFunctionsFactory::get().describe(sf_name), ' ')[0].data()));

    c.Save("pdf,png");
  }
}
