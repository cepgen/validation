#include "Canvas.h"

R__ADD_INCLUDE_PATH( /home/forthomme/work/dev/cepgen-private )
R__ADD_LIBRARY_PATH( /home/forthomme/work/dev/cepgen-private/build )
R__LOAD_LIBRARY( libgslcblas.so )
R__LOAD_LIBRARY( libgsl.so )
R__LOAD_LIBRARY( libmuparser.so )
R__LOAD_LIBRARY( libCepGenCore.so )

#include "CepGen/Generator.h"
#include "CepGen/Core/Exception.h"
#include "CepGen/Physics/MCDFileParser.h"
#include "CepGen/Modules/StructureFunctionsFactory.h"
#include "CepGen/FormFactors/Parameterisation.h"

#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystem.h"

void compare_ff( bool logx = false, bool logy = false, bool right = false, bool show_legend = true )
{
  const double q2_min = 1., q2_max = 1.e5;
  const double lq2min = log10( q2_min ), lq2max = log10( q2_max );
  const unsigned short num_points = 1000;

  pdg::MCDFileParser::parse( "/home/forthomme/work/dev/cepgen-private/External/mass_width_2020.mcd" );
  //cepgen::initialise();
  //cepgen::loadLibrary( "/home/forthomme/work/dev/cepgen-private/build/libCepGenCore.so" );
  //cepgen::dumpModules();

  cepgen::Canvas c( "ff_comparison", "CepGen form factors" );
  c.SetLegendX1( 0.4 );
  c.SetLegendY1( 0.77 );

  TMultiGraph mg;

  vector<int> colours = { kBlack, kRed+1, kBlue+1, kGreen+2, kOrange+1, kMagenta+1 };
  vector<TGraph*> v_g_fe, v_g_fm;
  vector<pair<string,string> > v_p_ffnames = {
    { "StandardDipole", "Standard dipole" },
    { "Mergell", "Mergell et al." },
    { "Brash", "Brash et al." },
    { "Arrington", "Arrington et al." }
  };
  size_t j = 0;
  CG_INFO("")<<v_p_ffnames;
  CG_INFO("")<<cepgen::formfac::FormFactorsFactory::get().modules();
  cout << "haha" << endl;
  for ( const auto& p_ffnames : v_p_ffnames ) {
    const auto& ffmod = cepgen::formfac::FormFactorsFactory::get().build( p_ffnames.first );
    v_g_fe.emplace_back();
    auto& gr_fe = *v_g_fe.rbegin();
    gr_fe->SetLineColor( colours[j] );
    gr_fe->SetLineWidth( 2 );
    gr_fe->SetLineStyle( 1 );
    mg.Add( gr_fe );
    v_g_fm.emplace_back();
    auto& gr_fm = *v_g_fm.rbegin();
    gr_fm->SetLineColor( colours[j] );
    gr_fm->SetLineWidth( 2 );
    gr_fm->SetLineStyle( 2 );
    mg.Add( gr_fm );
    if ( show_legend )
      c.AddLegendEntry( gr_fe, p_ffnames.second.c_str(), "l" );
    for ( unsigned int i = 0; i < num_points; ++i ) {
      double q2 = q2_min + ( q2_max-q2_min ) * i / ( num_points-1 );
      if ( logx )
        q2 = pow( 10, lq2min + i*( lq2max-lq2min )/( num_points-1 ) );
      auto fn = (*ffmod)( cepgen::mode::Beam::ProtonElastic, q2 );
      //cout << q2 << "|" << xbj << "|" << fn << endl;
      gr_fe->SetPoint( i, q2, fn.FE );
      gr_fm->SetPoint( i, q2, fn.FM );
    }
    ++j;
  }
  cout << "haha" << endl;

  //--- general plotting

  mg.Draw( "al" );
  mg.GetHistogram()->SetTitle( "Q^{2}\\F_{E,M}^{p}" );
  c.Prettify( mg.GetHistogram() );
  mg.GetHistogram()->SetTitle( "" );

  mg.GetHistogram()->GetXaxis()->SetLimits( q2_min, q2_max );
  if ( logx ) c.SetLogx();
  if ( logy ) {
    mg.GetHistogram()->SetMinimum( 1.e-3 );
    mg.GetHistogram()->SetMaximum( TMath::Min( 5., mg.GetHistogram()->GetMaximum() ) );
    c.SetLogy();
  }
  //else mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, 1.5 );
  else {
    //mg.GetHistogram()->SetMinimum( 1.e-4 );
    mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, TMath::Min( 1.5, mg.GetHistogram()->GetMaximum() ) );
    //mg.GetHistogram()->GetYaxis()->SetRangeUser( 1.e-4, TMath::Min( 1.9, mg.GetHistogram()->GetMaximum() ) );
  }

  c.Save( "pdf" );
}
