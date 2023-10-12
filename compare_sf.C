// clang-format off
#include "StrFunParams.h"
// clang-format on

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TSystem.h>

#include "CepGen/Generator.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/StructureFunctions/Parameterisation.h"
#include "CepGen/Utils/Environment.h"
#include "CepGen/Utils/String.h"
#include "CepGenAddOns/ROOTWrapper/ROOTCanvas.h"
#include "SampleHandler.h"

void compare_sf(double q2 = 1.225,
                bool logx = false,
                bool logy = false,
                bool right = false,
                bool show_legend = true,
                bool plot_fl = false,
                const string& output_file = "sf_comparison",
                const vector<string>& formats = {"pdf"}) {
  //const cepgen::Limits xbj_lims(2.5e-6, 0.99);
  const cepgen::Limits xbj_lims(1.1e-5, 0.99);
  //const cepgen::Limits xbj_lims(1.e-3, 0.99);
  //const cepgen::Limits xbj_lims(1.e-3, 0.9);
  const unsigned short num_points = 1000;
  const double x_pos_leg_data = (right) ? 0.56 : 0.16;
  //const double x_pos_leg_data = (right) ? 0.5 : 0.16;
  const double y_pos_leg_data = 0.35;
  const double x_leg_size = 0.15, y_leg_size = 0.15;
  //const double y_pos_leg_data = 0.7;
  const string mstw_grid_path = cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mstw_sf_scan_nnlo.dat";

  cepgen::initialise();
  pdg::MCDFileParser::parse(cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/mass_width_2021.mcd");

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
  //sfs.emplace_back(402);
  //sfs.emplace_back(205);
  //sfs.emplace_back(206);
  //sfs.emplace_back(401, cepgen::ParametersList().set<string>("pdfSet", "cteq6l1"), "CTEQ6l1", kGray + 1);
  /*sfs.emplace_back(401, cepgen::ParametersList().set<string>("pdfSet", "cteq6l1"), "CTEQ6l1:full", kBlue);
  sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "cteq6l1").set<int>("mode", 1), "CTEQ6l1:valence", kBlue, 2);
  sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "cteq6l1").set<int>("mode", 2), "CTEQ6l1:sea", kBlue, 3);*/
  /*sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "LUXqed17_plus_PDF4LHC15_nnlo_100"), "LUXqed17", kBlue - 2, 3);*/
  /*sfs.emplace_back(205, cepgen::ParametersList().set<string>("gridPath", mstw_grid_path), "MSTW:grid-NNLO", kOrange);
  sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "MSTW2008lo90cl"), "MSTW:partonic-LO", kOrange, 2);
  sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "MSTW2008nlo90cl"), "MSTW:partonic-NLO", kOrange, 3);
  sfs.emplace_back(
      401, cepgen::ParametersList().set<string>("pdfSet", "MSTW2008nnlo90cl"), "MSTW:partonic-NNLO", kOrange, 5);*/
  sfs.emplace_back(303,
                   cepgen::ParametersList().set<string>(
                       "gridFile", cepgen::utils::env::get("HOME") + "/work/dev/cepgen/External/a08tmc.dat"));

  for (const auto& xbj : xbj_lims.generate(num_points, logx))
    for (auto& sf : sfs) {
      if (plot_fl)
        sf.graph().AddPoint(xbj, sf.sf()->FL(xbj, q2));
      else
        sf.graph().AddPoint(xbj, sf.sf()->F2(xbj, q2));
    }

  cepgen::ROOTCanvas c(output_file);
  c.SetTopLabel(Form("Q^{2} ~ %g GeV^{2}", q2));
  c.SetLegendX1(0.4);
  c.SetLegendY1(0.77);
  c.SetGrid(1, 1);

  TMultiGraph mg;

  //--- model curves

  for (auto& sf : sfs) {
    if (sf.name().find("LUXlike") != std::string::npos) {
      auto gr_luxlike_line = (TGraph*)sf.graph().Clone();
      gr_luxlike_line->SetLineColorAlpha(kRed + 1, 0.5);
      gr_luxlike_line->SetLineWidth(8);
      mg.Add(gr_luxlike_line);
    }
    mg.Add(&sf.graph());
    if (show_legend)
      c.AddLegendEntry(&sf.graph(), cepgen::utils::split(sf.name(), ' ')[0].data(), "l");
  }

  //--- data points

  //auto leg_data = new TLegend( x_pos_leg_data, 0.3, x_pos_leg_data+0.15, 0.55, "" );
  auto leg_data =
      new TLegend(x_pos_leg_data, y_pos_leg_data, x_pos_leg_data + x_leg_size, y_pos_leg_data + y_leg_size, "");

  struct DataSample {
    explicit DataSample(const SampleHandler& hnd, const std::string& aname, int astyle)
        : handler(hnd), name(aname), style(astyle) {}
    SampleHandler handler;
    std::string name;
    int style{0};
  };
  std::vector<DataSample> samples;
  samples.emplace_back(SampleHandler::fromCLAS("samples/clas_0p225-4p725gev2.csv"), "CLAS", 24);
  samples.emplace_back(SampleHandler::fromBCDMS("samples/bcdms.csv"), "BCDMS", 30);
  samples.emplace_back(
      SampleHandler::fromCLAS("samples/nmc_0p75-65gev2.csv") + SampleHandler::fromNMC("samples/nmc_0p8-62gev2.csv"),
      "NMC",
      25);
  samples.emplace_back(SampleHandler::fromCLAS("samples/e665.csv"), "E665", 27);
  samples.emplace_back(SampleHandler::fromZEUS("samples/zeus.csv") +
                           SampleHandler::fromCLAS("samples/zeus_0p11-0p65gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_1p5-15gev2.csv") +
                           SampleHandler::fromZEUS("samples/zeus_0p6-17gev2.csv") +
                           SampleHandler::fromCLAS("samples/zeus_8p5-5000gev2.csv"),
                       "ZEUS",
                       26);
  samples.emplace_back(SampleHandler::fromCLAS("samples/h1.csv") +
                           SampleHandler::fromCLAS("samples/h1_4p5-1600gev2.csv") +
                           //SampleHandler::fromCLAS("samples/h1_neutral_2004.csv") +
                           //SampleHandler::fromCLAS("samples/hera-h1-nc.csv") +
                           SampleHandler::fromH11997("samples/h1_1p5-150gev2.csv") +
                           SampleHandler::fromH1lowQ2("samples/h1_compton.csv") +
                           SampleHandler::fromH1lowQ2("samples/h1_0p35-3p5gev2.csv"),
                       "H1",
                       22);
  samples.emplace_back(SampleHandler::fromHermes("samples/hermes.csv"), "HERMES", 28);

  if (!plot_fl)
    for (const auto& sample : samples) {
      auto* proj = new TGraphErrors(sample.handler.xBjGraph(q2));
      proj->SetMarkerStyle(sample.style);
      proj->SetMarkerSize(1.);
      mg.Add(proj, "p");
      if (proj->GetN() > 0)
        leg_data->AddEntry(proj, sample.name.data(), "ep");
    }

  //--- general plotting

  mg.Draw("al");
  mg.GetHistogram()->SetTitle(Form("x_{Bj}\\%s", (plot_fl ? "F_{L}^{p}" : "F_{2}^{p}")));
  c.Prettify(mg.GetHistogram());
  mg.GetHistogram()->SetTitle("");

  mg.GetHistogram()->GetXaxis()->SetLimits(xbj_lims.min(), xbj_lims.max());
  if (logx)
    c.SetLogx();
  if (logy) {
    mg.GetHistogram()->SetMinimum(1.e-3);
    mg.GetHistogram()->SetMaximum(TMath::Min(5., mg.GetHistogram()->GetMaximum()));
    c.SetLogy();
  }
  //else mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, 1.5 );
  else {
    //mg.GetHistogram()->SetMinimum( 1.e-4 );
    mg.GetHistogram()->GetYaxis()->SetRangeUser(1.e-4, TMath::Min(1.5, mg.GetHistogram()->GetMaximum()));
    //mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, TMath::Min( 1.9, mg.GetHistogram()->GetMaximum() ) );
  }

  if (auto* leg = c.GetLegend()) {
    if (sfs.size() > 6) {
      leg->SetNColumns(2);
      leg->SetX1(0.18);
      leg->SetX2(0.88);
      leg->SetY1(0.62);
      leg->SetY2(0.9);
    } else if (sfs.size() > 9) {
      leg->SetNColumns(3);
      leg->SetX1(0.18);
      leg->SetX2(0.88);
      leg->SetY1(0.77);
      leg->SetY2(0.92);
    }
  }
  if (leg_data->GetNRows() > 0) {
    leg_data->SetTextFont(132);
    leg_data->SetTextSize(0.035);
    leg_data->SetLineWidth(0);
    leg_data->SetFillStyle(0);
    leg_data->Draw();
  }

  for (const auto& fmt : formats)
    c.Save(fmt);
}
