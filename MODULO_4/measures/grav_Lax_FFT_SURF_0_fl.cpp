#ifdef __CLING__
#pragma cling optimize(0)
#endif
void grav_Lax_FFT_SURF_0_fl()
{
//=========Macro generated from canvas: canvas2/Canvas2
//=========  (Fri Jun 17 17:25:42 2022) by ROOT version 6.26/00
   TCanvas *canvas2 = new TCanvas("canvas2", "Canvas2",0,0,600,400);
   canvas2->SetHighLightColor(2);
   canvas2->Range(-0.353858,-36.68524,2.727147,5.989447);
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
   multigraph->SetTitle("grav_Lax_FFT_SURF_0;k;|u(k,t)|^2");
   
   Double_t Graph_fx13[250] = {
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
   Double_t Graph_fy13[250] = {
   7.115849,
   5.292213,
   3.230885,
   1.619122,
   0.6660567,
   0.224914,
   0.06234414,
   0.01418562,
   0.002649567,
   0.000406233,
   5.112685e-05,
   5.281979e-06,
   4.479381e-07,
   3.118265e-08,
   1.781893e-09,
   8.35842e-11,
   3.218402e-12,
   1.017256e-13,
   2.63933e-15,
   5.621223e-17,
   9.827506e-19,
   1.41032e-20,
   1.661553e-22,
   1.608102e-24,
   1.251065e-26,
   1.149288e-28,
   8.737373e-31,
   1.874128e-30,
   2.074575e-30,
   1.505644e-30,
   1.73608e-30,
   1.548416e-30,
   1.913432e-30,
   1.434191e-30,
   1.457521e-30,
   1.739781e-30,
   1.137339e-30,
   1.440771e-30,
   1.205388e-30,
   1.268723e-30,
   1.104455e-30,
   1.13895e-30,
   1.443313e-30,
   9.204016e-31,
   1.268281e-30,
   1.122704e-30,
   1.056605e-30,
   8.566467e-31,
   9.462538e-31,
   9.892346e-31,
   7.61195e-31,
   1.022959e-30,
   7.115618e-31,
   8.422116e-31,
   8.367423e-31,
   7.139508e-31,
   8.945759e-31,
   6.275514e-31,
   8.478385e-31,
   7.1039e-31,
   5.883026e-31,
   7.504838e-31,
   6.538887e-31,
   6.591634e-31,
   5.808641e-31,
   6.461461e-31,
   5.837001e-31,
   4.834751e-31,
   6.100731e-31,
   4.840565e-31,
   4.672996e-31,
   4.150736e-31,
   5.174933e-31,
   4.963187e-31,
   4.760522e-31,
   5.428517e-31,
   4.208854e-31,
   3.969034e-31,
   4.224831e-31,
   4.941285e-31,
   4.140012e-31,
   3.327391e-31,
   5.359954e-31,
   3.818302e-31,
   3.759495e-31,
   4.314699e-31,
   4.09706e-31,
   3.752147e-31,
   3.155975e-31,
   3.868466e-31,
   3.020632e-31,
   3.393898e-31,
   3.837489e-31,
   3.320852e-31,
   3.83941e-31,
   2.768323e-31,
   3.154302e-31,
   2.319974e-31,
   2.738292e-31,
   3.836917e-31,
   2.144781e-31,
   2.319974e-31,
   2.799715e-31,
   3.031648e-31,
   2.356986e-31,
   2.944035e-31,
   2.686414e-31,
   2.165012e-31,
   2.425209e-31,
   2.17822e-31,
   3.002528e-31,
   2.230166e-31,
   3.10356e-31,
   2.88793e-31,
   1.801608e-31,
   2.756348e-31,
   2.182474e-31,
   2.090861e-31,
   2.424413e-31,
   2.727784e-31,
   2.569913e-31,
   2.064799e-31,
   2.227146e-31,
   2.248754e-31,
   1.814627e-31,
   2.317741e-31,
   2.202075e-31,
   1.702632e-31,
   1.868017e-31,
   2.456688e-31,
   1.933407e-31,
   1.764606e-31,
   2.053416e-31,
   1.669414e-31,
   1.739593e-31,
   1.801748e-31,
   1.872756e-31,
   2.215012e-31,
   1.420517e-31,
   1.9472e-31,
   1.215907e-31,
   1.976161e-31,
   1.562488e-31,
   2.08669e-31,
   1.871677e-31,
   1.960577e-31,
   1.897062e-31,
   7.033726e-32,
   2.090157e-31,
   1.672355e-31,
   1.243186e-31,
   2.018722e-31,
   1.095341e-31,
   1.817466e-31,
   1.378401e-31,
   1.809985e-31,
   1.839838e-31,
   1.527099e-31,
   1.490591e-31,
   1.307661e-31,
   1.742236e-31,
   1.213821e-31,
   1.482825e-31,
   1.865334e-31,
   1.18036e-31,
   1.937038e-31,
   1.22179e-31,
   1.849581e-31,
   1.349239e-31,
   1.52047e-31,
   1.645462e-31,
   7.576123e-32,
   1.581921e-31,
   1.304821e-31,
   1.416731e-31,
   1.199193e-31,
   1.277134e-31,
   1.509989e-31,
   1.477917e-31,
   1.327922e-31,
   1.449109e-31,
   1.352011e-31,
   1.171886e-31,
   1.246825e-31,
   1.641671e-31,
   7.58336e-32,
   1.35891e-31,
   1.050439e-31,
   1.363479e-31,
   1.571234e-31,
   1.07692e-31,
   1.626935e-31,
   1.138145e-31,
   1.389821e-31,
   9.171456e-32,
   1.589238e-31,
   9.575426e-32,
   1.325549e-31,
   3.231417e-31,
   2.284644e-31,
   4.191739e-31,
   1.325549e-31,
   1.464171e-31,
   2.270393e-31,
   8.547018e-32,
   1.286538e-31,
   1.131132e-31,
   1.463776e-31,
   9.59061e-32,
   1.095082e-31,
   9.970448e-32,
   9.154539e-32,
   1.124647e-31,
   1.03086e-31,
   1.170627e-31,
   1.261193e-31,
   9.532574e-32,
   1.359525e-31,
   1.166289e-31,
   1.058395e-31,
   1.531699e-31,
   8.083475e-32,
   1.460267e-31,
   9.296653e-32,
   1.262018e-31,
   1.016879e-31,
   1.271803e-31,
   1.137932e-31,
   1.388899e-31,
   1.144033e-31,
   1.305811e-31,
   1.45352e-31,
   6.503909e-32,
   1.086264e-31,
   1.281954e-31,
   8.49799e-32,
   1.318847e-31,
   1.093413e-31,
   1.243512e-31,
   8.087137e-32,
   1.282951e-31,
   1.415528e-31,
   1.017892e-31,
   1.227911e-31,
   1.074583e-31,
   1.155681e-31,
   1.009742e-31,
   1.072851e-31,
   1.009742e-31,
   1.009742e-31};
   TGraph *graph = new TGraph(250,Graph_fx13,Graph_fy13);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   
   TH1F *Graph_Graph13 = new TH1F("Graph_Graph13","Graph",250,0.9,274.9);
   Graph_Graph13->SetMinimum(5.853518e-32);
   Graph_Graph13->SetMaximum(7.827434);
   Graph_Graph13->SetDirectory(0);
   Graph_Graph13->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph13->SetLineColor(ci);
   Graph_Graph13->GetXaxis()->SetLabelFont(42);
   Graph_Graph13->GetXaxis()->SetTitleOffset(1);
   Graph_Graph13->GetXaxis()->SetTitleFont(42);
   Graph_Graph13->GetYaxis()->SetLabelFont(42);
   Graph_Graph13->GetYaxis()->SetTitleFont(42);
   Graph_Graph13->GetZaxis()->SetLabelFont(42);
   Graph_Graph13->GetZaxis()->SetTitleOffset(1);
   Graph_Graph13->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph13);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx14[250] = {
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
   Double_t Graph_fy14[250] = {
   7.116834,
   5.295032,
   3.230133,
   1.617008,
   0.6642067,
   0.2239075,
   0.06193726,
   0.01406024,
   0.002619067,
   0.0004003167,
   5.02251e-05,
   5.167128e-06,
   4.369131e-07,
   3.026035e-08,
   1.711282e-09,
   8.353318e-11,
   2.555812e-12,
   1.717918e-13,
   1.388346e-14,
   8.931814e-15,
   7.254992e-15,
   5.362721e-15,
   4.122724e-15,
   3.19474e-15,
   2.50286e-15,
   1.980137e-15,
   1.580847e-15,
   1.272724e-15,
   1.032697e-15,
   8.440636e-16,
   6.945996e-16,
   5.752596e-16,
   4.792846e-16,
   4.015773e-16,
   3.382596e-16,
   2.863569e-16,
   2.4357e-16,
   2.081084e-16,
   1.785688e-16,
   1.538435e-16,
   1.330529e-16,
   1.154949e-16,
   1.006051e-16,
   8.792819e-17,
   7.709441e-17,
   6.780236e-17,
   5.980506e-17,
   5.289934e-17,
   4.691727e-17,
   4.171951e-17,
   3.719004e-17,
   3.323183e-17,
   2.976348e-17,
   2.671644e-17,
   2.403283e-17,
   2.166357e-17,
   1.956698e-17,
   1.770748e-17,
   1.605468e-17,
   1.45825e-17,
   1.326855e-17,
   1.20935e-17,
   1.104067e-17,
   1.00956e-17,
   9.245738e-18,
   8.480174e-18,
   7.789375e-18,
   7.165028e-18,
   6.599838e-18,
   6.087408e-18,
   5.622119e-18,
   5.199017e-18,
   4.813728e-18,
   4.462389e-18,
   4.141578e-18,
   3.848258e-18,
   3.579727e-18,
   3.33359e-18,
   3.107699e-18,
   2.900147e-18,
   2.709224e-18,
   2.5334e-18,
   2.371301e-18,
   2.221695e-18,
   2.083475e-18,
   1.955644e-18,
   1.837301e-18,
   1.727635e-18,
   1.625911e-18,
   1.531467e-18,
   1.4437e-18,
   1.362066e-18,
   1.28607e-18,
   1.215259e-18,
   1.149226e-18,
   1.087596e-18,
   1.030029e-18,
   9.762166e-19,
   9.258732e-19,
   8.787385e-19,
   8.345761e-19,
   7.931676e-19,
   7.543135e-19,
   7.178315e-19,
   6.835512e-19,
   6.513197e-19,
   6.209931e-19,
   5.924405e-19,
   5.655413e-19,
   5.401819e-19,
   5.162622e-19,
   4.936845e-19,
   4.723617e-19,
   4.522122e-19,
   4.3316e-19,
   4.151354e-19,
   3.98074e-19,
   3.81914e-19,
   3.666009e-19,
   3.520825e-19,
   3.383102e-19,
   3.252383e-19,
   3.128266e-19,
   3.010335e-19,
   2.898249e-19,
   2.791662e-19,
   2.690244e-19,
   2.593714e-19,
   2.501795e-19,
   2.414214e-19,
   2.330743e-19,
   2.251151e-19,
   2.175221e-19,
   2.102759e-19,
   2.033581e-19,
   1.967499e-19,
   1.904368e-19,
   1.844011e-19,
   1.786301e-19,
   1.731094e-19,
   1.678264e-19,
   1.627687e-19,
   1.579255e-19,
   1.532853e-19,
   1.488389e-19,
   1.445759e-19,
   1.40488e-19,
   1.365666e-19,
   1.328034e-19,
   1.291914e-19,
   1.257227e-19,
   1.223918e-19,
   1.191907e-19,
   1.161145e-19,
   1.131575e-19,
   1.103135e-19,
   1.075786e-19,
   1.049465e-19,
   1.024139e-19,
   9.997576e-20,
   9.762785e-20,
   9.536727e-20,
   9.31887e-20,
   9.108958e-20,
   8.906694e-20,
   8.711637e-20,
   8.523647e-20,
   8.342257e-20,
   8.167317e-20,
   7.998543e-20,
   7.835651e-20,
   7.678436e-20,
   7.526662e-20,
   7.380131e-20,
   7.238602e-20,
   7.101955e-20,
   6.969904e-20,
   6.842354e-20,
   6.719075e-20,
   6.599991e-20,
   6.48484e-20,
   6.37356e-20,
   6.265981e-20,
   6.161935e-20,
   6.061391e-20,
   5.964113e-20,
   5.87006e-20,
   5.77911e-20,
   5.691113e-20,
   5.606047e-20,
   5.523706e-20,
   5.444102e-20,
   5.367092e-20,
   5.292579e-20,
   5.220532e-20,
   5.150825e-20,
   5.083416e-20,
   5.018243e-20,
   4.955132e-20,
   4.894223e-20,
   4.835241e-20,
   4.778225e-20,
   4.72314e-20,
   4.66987e-20,
   4.618408e-20,
   4.568663e-20,
   4.520638e-20,
   4.474241e-20,
   4.429462e-20,
   4.386233e-20,
   4.344533e-20,
   4.304291e-20,
   4.265517e-20,
   4.228112e-20,
   4.192147e-20,
   4.157439e-20,
   4.124122e-20,
   4.092025e-20,
   4.061202e-20,
   4.031617e-20,
   4.003202e-20,
   3.975978e-20,
   3.949889e-20,
   3.924956e-20,
   3.901111e-20,
   3.878348e-20,
   3.856686e-20,
   3.836011e-20,
   3.816426e-20,
   3.797818e-20,
   3.780237e-20,
   3.763614e-20,
   3.747995e-20,
   3.733316e-20,
   3.719596e-20,
   3.706805e-20,
   3.69495e-20,
   3.683992e-20,
   3.673936e-20,
   3.664818e-20,
   3.656576e-20,
   3.649192e-20,
   3.642733e-20,
   3.637094e-20,
   3.632359e-20,
   3.62853e-20,
   3.625469e-20,
   3.623369e-20,
   3.622042e-20,
   3.619457e-20};
   graph = new TGraph(250,Graph_fx14,Graph_fy14);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   graph->SetLineColor(2);
   graph->SetMarkerColor(2);
   
   TH1F *Graph_Graph14 = new TH1F("Graph_Graph14","Graph",250,0.9,274.9);
   Graph_Graph14->SetMinimum(3.257511e-20);
   Graph_Graph14->SetMaximum(7.828517);
   Graph_Graph14->SetDirectory(0);
   Graph_Graph14->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph14->SetLineColor(ci);
   Graph_Graph14->GetXaxis()->SetLabelFont(42);
   Graph_Graph14->GetXaxis()->SetTitleOffset(1);
   Graph_Graph14->GetXaxis()->SetTitleFont(42);
   Graph_Graph14->GetYaxis()->SetLabelFont(42);
   Graph_Graph14->GetYaxis()->SetTitleFont(42);
   Graph_Graph14->GetZaxis()->SetLabelFont(42);
   Graph_Graph14->GetZaxis()->SetTitleOffset(1);
   Graph_Graph14->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph14);
   
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
   entry=leg->AddEntry("Graph","t = 0.998","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(2);
   entry->SetMarkerStyle(1);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.3121477,0.9304839,0.6878523,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("grav_Lax_FFT_SURF_0");
   pt->Draw();
   canvas2->Modified();
   canvas2->cd();
   canvas2->SetSelected(canvas2);
}
