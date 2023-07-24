
#include "TApplication.h"

#include "iostream"
#include "iomanip"

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TRandom.h"

using namespace std;

Double_t _gaussexp( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];

  Double_t Norm = par[0];
  Double_t Q = par[1];
  Double_t s = par[2];
  Double_t alpha = par[3];
  Double_t w = par[4];
  
  Double_t arg = 0.0; 
  
  if ( s!=0.0 ) arg = ( xx - Q )/s;    
  else cout << "Error: The code tries to divide by zero." << endl;

  Double_t gn = 0.5*TMath::Erfc( -Q/( sqrt(2.0)*s ) );
  
  Double_t result = 0.0;
  if ( xx>=0.0 ) result = w*alpha*TMath::Exp( -xx*alpha ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );

  result *= Norm;
  
  return result;
  
}


Int_t example()
{  
  time_t start;  
  time( &start );

  cout << "" << endl;
  cout << " the macro now begins ... ( example.cc )" << endl;
  cout << "" << endl;
  
  gROOT->Reset();
  gStyle->SetCanvasBorderMode( 0 );
  gStyle->SetCanvasColor( 0 );
  gStyle->SetPadColor( 0 );
  gStyle->SetPadBorderMode( 0 );
  gStyle->SetFrameBorderMode( 0 );
  gStyle->SetTitleColor( 0 );   
  gStyle->SetTitleFillColor( 0 );  
  gStyle->SetTitleBorderSize( 0 );
  gStyle->SetTitleX( 0.10 );
  gStyle->SetTitleY( 0.98 );
  gStyle->SetTitleFont( 22, "" );
  gStyle->SetTitleSize( 0.055, ""  );
  gStyle->SetStatColor( 0 );
  gStyle->SetStatFont( 22 );
  gStyle->SetStatBorderSize( 1 );
  gStyle->SetStatX( 0.90 );
  gStyle->SetStatY( 0.90 );
  gStyle->SetStatFontSize( 0.04 );
  gStyle->SetOptStat( 1110 );
  gStyle->SetTitleFont( 22, "XYZ"  );
  gStyle->SetTitleSize( 0.05, "XYZ"  );
  gStyle->SetTitleColor( kBlack, "XYZ"  );
  gStyle->SetTitleAlign(13);
  gStyle->SetLabelFont( 22, "XYZ"  );
  gStyle->SetLabelSize( 0.04, "XYZ"  );
  gStyle->SetOptStat( 1110 ); 
  gStyle->SetOptFit( 0 ); 
  
  gStyle->SetPalette( 1 ); 
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable( NRGBs, stops, red, green, blue, NCont );
  gStyle->SetNumberContours( NCont );

  gROOT->ForceStyle();
    
  
  // Input <--------------------- !!!

  Int_t ntot = 2.5e+6;
  
  Float_t Q0 = 0.0;
  Float_t s0 = 20.0;
    
  Float_t mu = 1.0;
  
  Float_t Q = 100.0;
  Float_t s = 25.0;
  Float_t alpha = 1.0/30.0;
  Float_t w = 0.2;
  
    
  // <--------------------------- !!!
  
  TCanvas *c1 = new TCanvas( "c1", "" );
  c1->cd();
  //c1->SetLogy();
  
  TF1 *spe = new TF1( "spe" , _gaussexp, 0.0, 2000.0, 5 );
  spe->SetParNames( "Norm", "Q", "#sigma", "#alpha", "w" );
  spe->SetNpx( 2000 );
  spe->SetLineColor( kRed );
  spe->SetLineWidth( 2 );
  
  spe->SetParameters( 1.0, Q, s, alpha, w );

  
  TH1D *h1 = new TH1D( "h1", "; Charge [AU]; Entries", 250, -200.0, 800.0 );
  gRandom->SetSeed( 0 );
  
  for ( Int_t i=0; i<ntot; i++ )
    {
      Double_t QQ = gRandom->Gaus( Q0, s0 );
      
      Int_t n = gRandom->Poisson( mu );
      for ( Int_t j=0; j<n; j++ ) QQ += spe->GetRandom();

      h1->Fill( QQ );
      
    }

  h1->Draw( "" );
  c1->Update();
  c1->WaitPrimitive();
  
  cout << "" << endl;
  cout << "" << endl;
  
  time_t end;      
  time( &end );
      
  Int_t dura = difftime( end, start );      
  Int_t min = dura / 60; Int_t sec = dura % 60;
    
  cout << " ---> "<< Form( "%02d:%02d", min, sec ) << endl;  
  cout << "" << endl;
  
  return 0;
  
}

int main()
{
  TApplication theApp( "App", 0, 0 );
    
  example();

}
