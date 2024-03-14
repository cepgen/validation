// clang-format off
#include "StrFunParams.h"
// clang-format on

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TSystem.h>

#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Physics/PDG.h"
#include "CepGen/Physics/Utils.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/Message.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "DISsamples/DISsample.h"

void compare_sf_vs_w(double q2 = 1.225,
                     bool logx = false,
                     bool logy = false,
                     bool plot_fl = false,
                     bool plot_legend = true,
                     const string& output_file = "sf_comparison",
                     const vector<string>& formats = {"pdf"}) {
  pdg::MCDFileParser::parse(cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mass_width_2021.mcd");

  const auto mp = cepgen::PDG::get().mass(2212);
  cepgen::Limits w_lims{1., 10.};
  const unsigned short num_points = 5000;
  const double x_pos_leg_data = 0.16, y_pos_leg_data = 0.775;
  const double x_leg_size = 0.15, y_leg_size = 0.15;

  const string mstw_grid_path = cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mstw_sf_scan_nnlo.dat";

  vector<StrFunParams> sfs;
  sfs.emplace_back(
      301,
      cepgen::ParametersList().set<cepgen::ParametersList>(
          "perturbativeSF",
          cepgen::ParametersList().setName<int>(205 /* MSTWgrid */).set<string>("gridPath", mstw_grid_path))
      //.set<vector<double> >("W2limits", {3., 16.})
  );
  sfs.emplace_back(11);
  sfs.emplace_back(101);
  sfs.emplace_back(103);
  sfs.emplace_back(102);
  sfs.emplace_back(202);
  sfs.emplace_back(203);
  sfs.emplace_back(204);
  //sfs.emplace_back(304);

  struct DataSample {
    explicit DataSample(const DISsample& hnd, const std::string& aname, int astyle)
        : handler(hnd), name(aname), style(astyle) {}
    DISsample handler;
    std::string name;
    int style{0};
  };
  std::vector<DataSample> samples;
  samples.emplace_back(DISsample::fromCLAS("samples/clas_0p225-4p725gev2.csv"), "CLAS", 24);
  samples.emplace_back(DISsample::fromBCDMS("samples/bcdms.csv"), "BCDMS", 30);
  samples.emplace_back(
      DISsample::fromCLAS("samples/nmc_0p75-65gev2.csv") + DISsample::fromNMC("samples/nmc_0p8-62gev2.csv"), "NMC", 25);
  samples.emplace_back(DISsample::fromCLAS("samples/e665.csv"), "E665", 27);
  samples.emplace_back(DISsample::fromZEUS("samples/zeus.csv") + DISsample::fromCLAS("samples/zeus_0p11-0p65gev2.csv") +
                           DISsample::fromCLAS("samples/zeus_1p5-15gev2.csv") +
                           DISsample::fromZEUS("samples/zeus_0p6-17gev2.csv") +
                           DISsample::fromCLAS("samples/zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26);
  samples.emplace_back(DISsample::fromCLAS("samples/h1.csv") + DISsample::fromCLAS("samples/h1_4p5-1600gev2.csv") +
                           //DISsample::fromCLAS("samples/h1_neutral_2004.csv") +
                           //DISsample::fromCLAS("samples/hera-h1-nc.csv") +
                           DISsample::fromH11997("samples/h1_1p5-150gev2.csv") +
                           DISsample::fromH1lowQ2("samples/h1_compton.csv") +
                           DISsample::fromH1lowQ2("samples/h1_0p35-3p5gev2.csv"),
                       "H1",
                       22);
  samples.emplace_back(DISsample::fromHermes("samples/hermes.csv"), "HERMES", 28);
  //samples.emplace_back(DISsample::fromCCFR("samples/ccfr.txt"), "CCFR", 29);

  auto leg_data =
      new TLegend(x_pos_leg_data, y_pos_leg_data, x_pos_leg_data + x_leg_size, y_pos_leg_data + y_leg_size, "", "NDC");
  cepgen::Limits w_lims_samples{1., -1.};
  const double scaling_factor = 1.1;
  if (!plot_fl)
    for (const auto& sample : samples) {
      const auto* smp_x = sample.handler.wGraph(q2).GetXaxis();
      w_lims_samples.min() = TMath::Min(w_lims_samples.min(), smp_x->GetXmin() / scaling_factor);
      w_lims_samples.max() = TMath::Max(w_lims_samples.max(), smp_x->GetXmax() * scaling_factor);
    }

  w_lims.min() = TMath::Max(w_lims.min(), w_lims_samples.min());
  w_lims.max() = TMath::Max(w_lims.max(), w_lims_samples.max());
  CG_LOG << w_lims_samples << ":" << w_lims;

  for (const auto& w : w_lims.generate(num_points, logx)) {
    const auto xbj = cepgen::utils::xBj(q2, mp * mp, w * w);
    for (auto& sf : sfs) {
      if (plot_fl)
        sf.graph().AddPoint(w, sf.sf()->FL(xbj, q2));
      else
        sf.graph().AddPoint(w, sf.sf()->F2(xbj, q2));
    }
  }

  cepgen::ROOTCanvas c(output_file);
  c.SetTopLabel(Form("Q^{2} ~ %g GeV^{2}", q2));
  c.SetLegendX1(0.4);
  c.SetLegendY1(0.77);

  TMultiGraph mg;

  //--- model curves

  for (auto& sf : sfs) {
    if (sf.name().find("LUXlike") != std::string::npos) {
      auto gr_luxlike_line = (TGraph*)sf.graph().Clone();
      gr_luxlike_line->SetLineColorAlpha(kRed + 1, 0.5);
      gr_luxlike_line->SetLineWidth(8);
      mg.Add(gr_luxlike_line);
    }
    sf.graph().SetFillStyle(0);
    sf.graph().SetFillColor(0);
    mg.Add(&sf.graph());
    if (plot_legend)
      c.AddLegendEntry(&sf.graph(), cepgen::utils::split(sf.name(), ' ')[0].data(), "l");
  }

  if (sfs.size() > 9) {
    auto leg = c.GetLegend();
    leg->SetNColumns(3);
    leg->SetX1(0.18);
    leg->SetX2(0.88);
    leg->SetY1(0.77);
    leg->SetY2(0.92);
  }

  if (!plot_fl)
    for (const auto& sample : samples) {
      auto* proj = new TGraphErrors(sample.handler.wGraph(q2));
      proj->SetMarkerStyle(sample.style);
      proj->SetMarkerSize(1.);
      mg.Add(proj, "p");
      if (proj->GetN() > 0)
        leg_data->AddEntry(proj, sample.name.data(), "ep");
    }

  //--- general plotting

  mg.Draw("al");
  mg.GetHistogram()->SetTitle(Form("W (GeV)\\%s", (plot_fl ? "F_{L}^{p}" : "F_{2}^{p}")));
  c.Prettify(mg.GetHistogram());
  mg.GetHistogram()->SetTitle("");

  mg.GetHistogram()->GetXaxis()->SetLimits(w_lims.min(), w_lims.max());
  if (logx)
    c.SetLogx();
  if (logy) {
    mg.GetHistogram()->SetMinimum(1.e-3);
    c.SetLogy();
  } else {
    const double y_max = (plot_fl ? 5.0 : 1.5);
    mg.GetHistogram()->GetYaxis()->SetRangeUser(1.e-4, TMath::Min(y_max, mg.GetHistogram()->GetMaximum()));
  }
  if (leg_data->GetNRows() > 0) {
    leg_data->SetTextFont(132);
    leg_data->SetTextSize(0.035);
    leg_data->SetLineWidth(0);
    leg_data->SetFillStyle(0);
    c.Place(leg_data);
  }
  c.SetGrid(1, 1);
  for (const auto& fmt : formats)
    c.Save(fmt);
}
