#ifdef __CLING__
#pragma cling optimize(0)
#endif
void first_LeapFrog__test_FFT_SURF_0_fl()
{
//=========Macro generated from canvas: canvas2/Canvas2
//=========  (Fri Jun 17 10:03:07 2022) by ROOT version 6.26/00
   TCanvas *canvas2 = new TCanvas("canvas2", "Canvas2",0,0,600,400);
   canvas2->SetHighLightColor(2);
   canvas2->Range(-0.353858,-41.65123,2.727147,7.990372);
   canvas2->SetFillColor(0);
   canvas2->SetBorderMode(0);
   canvas2->SetBorderSize(2);
   canvas2->SetLogx();
   canvas2->SetLogy();
   canvas2->SetGridx();
   canvas2->SetGridy();
   canvas2->SetFrameBorderMode(0);
   canvas2->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("first_LeapFrog__test_FFT_SURF_0;k;|u(k,t)|^2");
   
   Double_t Graph_fx1[250] = {
   1,
   2,
   3,
   4,
   5,
   6,
   7,
   8,
   9,
   10,
   11,
   12,
   13,
   14,
   15,
   16,
   17,
   18,
   19,
   20,
   21,
   22,
   23,
   24,
   25,
   26,
   27,
   28,
   29,
   30,
   31,
   32,
   33,
   34,
   35,
   36,
   37,
   38,
   39,
   40,
   41,
   42,
   43,
   44,
   45,
   46,
   47,
   48,
   49,
   50,
   51,
   52,
   53,
   54,
   55,
   56,
   57,
   58,
   59,
   60,
   61,
   62,
   63,
   64,
   65,
   66,
   67,
   68,
   69,
   70,
   71,
   72,
   73,
   74,
   75,
   76,
   77,
   78,
   79,
   80,
   81,
   82,
   83,
   84,
   85,
   86,
   87,
   88,
   89,
   90,
   91,
   92,
   93,
   94,
   95,
   96,
   97,
   98,
   99,
   100,
   101,
   102,
   103,
   104,
   105,
   106,
   107,
   108,
   109,
   110,
   111,
   112,
   113,
   114,
   115,
   116,
   117,
   118,
   119,
   120,
   121,
   122,
   123,
   124,
   125,
   126,
   127,
   128,
   129,
   130,
   131,
   132,
   133,
   134,
   135,
   136,
   137,
   138,
   139,
   140,
   141,
   142,
   143,
   144,
   145,
   146,
   147,
   148,
   149,
   150,
   151,
   152,
   153,
   154,
   155,
   156,
   157,
   158,
   159,
   160,
   161,
   162,
   163,
   164,
   165,
   166,
   167,
   168,
   169,
   170,
   171,
   172,
   173,
   174,
   175,
   176,
   177,
   178,
   179,
   180,
   181,
   182,
   183,
   184,
   185,
   186,
   187,
   188,
   189,
   190,
   191,
   192,
   193,
   194,
   195,
   196,
   197,
   198,
   199,
   200,
   201,
   202,
   203,
   204,
   205,
   206,
   207,
   208,
   209,
   210,
   211,
   212,
   213,
   214,
   215,
   216,
   217,
   218,
   219,
   220,
   221,
   222,
   223,
   224,
   225,
   226,
   227,
   228,
   229,
   230,
   231,
   232,
   233,
   234,
   235,
   236,
   237,
   238,
   239,
   240,
   241,
   242,
   243,
   244,
   245,
   246,
   247,
   248,
   249,
   250};
   Double_t Graph_fy1[250] = {
   125,
   7.64304e-32,
   1.543575e-32,
   1.733074e-32,
   2.738148e-32,
   1.387432e-31,
   3.275623e-31,
   3.493625e-32,
   9.921242e-31,
   2.540719e-31,
   7.78999e-32,
   1.863226e-31,
   2.721406e-33,
   3.684179e-31,
   5.505719e-32,
   7.433157e-32,
   1.536481e-31,
   5.292467e-32,
   1.606185e-31,
   1.028188e-32,
   3.443872e-31,
   3.861012e-32,
   3.557274e-32,
   4.109615e-32,
   5.730403e-32,
   7.795222e-33,
   1.310926e-32,
   8.002545e-33,
   2.178493e-31,
   2.778473e-32,
   2.26781e-31,
   6.626651e-32,
   1.206654e-32,
   4.994256e-33,
   2.134222e-32,
   4.546896e-32,
   6.008766e-33,
   5.592484e-32,
   3.535335e-33,
   4.717113e-32,
   7.209756e-32,
   2.154205e-33,
   6.862362e-34,
   3.750054e-33,
   1.081471e-32,
   4.819468e-34,
   2.321836e-35,
   4.021729e-33,
   3.566843e-31,
   1.705901e-33,
   3.274045e-31,
   5.314645e-33,
   7.155917e-33,
   8.136628e-33,
   1.164685e-32,
   8.685834e-33,
   3.23452e-32,
   2.373808e-34,
   2.041745e-32,
   4.958777e-33,
   1.542553e-32,
   1.943525e-32,
   2.663346e-32,
   1.889711e-32,
   1.519592e-32,
   4.890768e-33,
   1.017049e-32,
   1.881115e-32,
   3.647599e-32,
   2.315012e-34,
   2.8517e-33,
   2.610509e-32,
   1.395747e-34,
   9.895069e-34,
   1.275706e-32,
   3.324223e-33,
   6.095191e-33,
   7.549952e-33,
   1.362406e-32,
   1.350589e-32,
   2.121249e-32,
   1.682837e-33,
   7.600084e-33,
   1.745769e-33,
   8.881419e-33,
   3.709781e-33,
   1.833296e-32,
   1.810276e-32,
   6.552732e-32,
   1.305473e-33,
   8.301869e-32,
   2.161978e-34,
   4.391079e-33,
   3.031182e-33,
   2.020589e-32,
   8.709236e-33,
   1.273388e-32,
   2.815971e-34,
   6.076192e-32,
   8.848344e-34,
   1.658366e-31,
   4.270255e-32,
   3.004398e-32,
   1.557747e-31,
   5.542349e-32,
   2.137085e-32,
   1.603465e-32,
   1.025433e-32,
   1.9467e-32,
   3.534986e-33,
   9.024764e-32,
   2.172572e-32,
   6.50695e-32,
   9.197939e-33,
   2.097422e-33,
   8.263836e-33,
   6.568201e-32,
   6.838956e-33,
   7.731074e-32,
   9.263869e-33,
   1.130898e-31,
   1.342975e-33,
   3.868292e-32,
   3.64475e-33,
   1.421584e-33,
   2.616301e-33,
   7.158457e-32,
   2.705894e-33,
   9.825838e-32,
   2.580359e-33,
   1.564978e-32,
   3.391947e-33,
   1.118953e-32,
   3.180191e-33,
   5.723102e-33,
   1.814893e-33,
   2.117562e-33,
   2.519595e-33,
   1.405021e-31,
   5.414397e-33,
   1.526218e-32,
   3.230285e-33,
   3.753225e-34,
   3.010592e-33,
   6.25045e-33,
   7.3579e-34,
   1.55813e-32,
   5.109298e-33,
   1.139826e-31,
   7.550141e-33,
   1.173232e-31,
   2.825645e-34,
   1.047949e-32,
   1.183309e-33,
   4.920382e-33,
   6.274973e-34,
   3.760556e-32,
   1.588146e-33,
   2.586723e-32,
   1.408301e-33,
   3.672416e-32,
   5.861254e-33,
   8.56658e-33,
   4.498658e-33,
   4.73295e-33,
   1.999644e-33,
   2.727163e-33,
   3.884341e-34,
   2.625571e-32,
   4.05855e-36,
   1.061454e-31,
   1.346151e-32,
   1.283625e-33,
   7.605728e-33,
   8.094471e-33,
   3.69978e-33,
   8.154331e-33,
   8.992776e-35,
   3.468943e-32,
   5.104694e-33,
   4.665746e-33,
   7.447639e-33,
   6.583472e-33,
   4.936618e-33,
   3.131064e-35,
   3.093481e-33,
   2.21827e-34,
   8.224351e-34,
   7.702852e-33,
   4.776215e-33,
   8.60905e-32,
   8.269004e-33,
   7.645579e-33,
   1.582513e-33,
   1.469746e-33,
   4.868502e-33,
   1.776983e-33,
   2.628448e-33,
   1.541266e-31,
   5.183146e-33,
   7.805851e-32,
   2.086252e-33,
   2.523528e-32,
   4.825318e-33,
   2.727826e-32,
   2.552289e-33,
   1.169774e-32,
   1.111769e-33,
   5.756277e-33,
   3.014845e-33,
   2.800687e-32,
   9.526337e-33,
   5.267031e-32,
   1.074979e-32,
   4.014709e-32,
   3.511514e-32,
   7.46102e-33,
   9.723369e-33,
   2.123752e-33,
   6.716067e-33,
   2.032039e-32,
   4.058063e-33,
   2.903007e-33,
   4.859126e-33,
   3.988122e-33,
   2.847968e-33,
   3.025179e-32,
   1.768031e-33,
   2.549316e-32,
   2.847131e-33,
   6.088894e-32,
   1.191218e-32,
   1.198317e-33,
   1.256005e-34,
   4.402954e-33,
   3.030313e-34,
   1.031442e-32,
   3.292339e-34,
   5.286714e-32,
   5.383027e-33,
   2.829921e-32,
   3.191087e-33,
   7.373301e-33,
   4.514959e-33,
   4.34175e-33,
   5.117749e-34,
   6.310607e-33,
   4.417781e-33,
   1.014255e-31,
   1.843211e-33};
   TGraph *graph = new TGraph(250,Graph_fx1,Graph_fy1);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",250,0.9,274.9);
   Graph_Graph1->SetMinimum(3.652695e-36);
   Graph_Graph1->SetMaximum(137.5);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph1);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx2[250] = {
   1,
   2,
   3,
   4,
   5,
   6,
   7,
   8,
   9,
   10,
   11,
   12,
   13,
   14,
   15,
   16,
   17,
   18,
   19,
   20,
   21,
   22,
   23,
   24,
   25,
   26,
   27,
   28,
   29,
   30,
   31,
   32,
   33,
   34,
   35,
   36,
   37,
   38,
   39,
   40,
   41,
   42,
   43,
   44,
   45,
   46,
   47,
   48,
   49,
   50,
   51,
   52,
   53,
   54,
   55,
   56,
   57,
   58,
   59,
   60,
   61,
   62,
   63,
   64,
   65,
   66,
   67,
   68,
   69,
   70,
   71,
   72,
   73,
   74,
   75,
   76,
   77,
   78,
   79,
   80,
   81,
   82,
   83,
   84,
   85,
   86,
   87,
   88,
   89,
   90,
   91,
   92,
   93,
   94,
   95,
   96,
   97,
   98,
   99,
   100,
   101,
   102,
   103,
   104,
   105,
   106,
   107,
   108,
   109,
   110,
   111,
   112,
   113,
   114,
   115,
   116,
   117,
   118,
   119,
   120,
   121,
   122,
   123,
   124,
   125,
   126,
   127,
   128,
   129,
   130,
   131,
   132,
   133,
   134,
   135,
   136,
   137,
   138,
   139,
   140,
   141,
   142,
   143,
   144,
   145,
   146,
   147,
   148,
   149,
   150,
   151,
   152,
   153,
   154,
   155,
   156,
   157,
   158,
   159,
   160,
   161,
   162,
   163,
   164,
   165,
   166,
   167,
   168,
   169,
   170,
   171,
   172,
   173,
   174,
   175,
   176,
   177,
   178,
   179,
   180,
   181,
   182,
   183,
   184,
   185,
   186,
   187,
   188,
   189,
   190,
   191,
   192,
   193,
   194,
   195,
   196,
   197,
   198,
   199,
   200,
   201,
   202,
   203,
   204,
   205,
   206,
   207,
   208,
   209,
   210,
   211,
   212,
   213,
   214,
   215,
   216,
   217,
   218,
   219,
   220,
   221,
   222,
   223,
   224,
   225,
   226,
   227,
   228,
   229,
   230,
   231,
   232,
   233,
   234,
   235,
   236,
   237,
   238,
   239,
   240,
   241,
   242,
   243,
   244,
   245,
   246,
   247,
   248,
   249,
   250};
   Double_t Graph_fy2[250] = {
   -nan,
   -nan,
   -nan,
   -nan,
   nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan,
   -nan};
   graph = new TGraph(250,Graph_fx2,Graph_fy2);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   graph->SetLineColor(2);
   graph->SetMarkerColor(2);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",250,0.9,274.9);
   Graph_Graph2->SetMinimum(-nan);
   Graph_Graph2->SetMaximum(-nan);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetTitleOffset(1);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph2);
   
   multigraph->Add(graph,"");
   multigraph->Draw("APL");
   multigraph->GetXaxis()->SetLimits(0.9, 262.45);
   multigraph->GetXaxis()->SetTitle("k");
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetTitleOffset(1);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("|u(k,t)|^2");
   multigraph->GetYaxis()->SetLabelFont(42);
   multigraph->GetYaxis()->SetTitleFont(42);
   
   TLegend *leg = new TLegend(0.75,0.75,0.85,0.85,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","t = 0.","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph","t = 0.299","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.2248993,0.9304839,0.7751007,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("first_LeapFrog__test_FFT_SURF_0");
   pt->Draw();
   canvas2->Modified();
   canvas2->cd();
   canvas2->SetSelected(canvas2);
}
