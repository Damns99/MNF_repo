#ifdef __CLING__
#pragma cling optimize(0)
#endif
void grav_LaxWendroff_SURF_0_fl()
{
//=========Macro generated from canvas: canvas2/Canvas2
//=========  (Fri Jun 17 17:25:45 2022) by ROOT version 6.26/00
   TCanvas *canvas2 = new TCanvas("canvas2", "Canvas2",0,0,600,400);
   canvas2->SetHighLightColor(2);
   canvas2->Range(-0.787125,-0.2099107,0.585125,1.191039);
   canvas2->SetFillColor(0);
   canvas2->SetBorderMode(0);
   canvas2->SetBorderSize(2);
   canvas2->SetGridx();
   canvas2->SetGridy();
   canvas2->SetFrameBorderMode(0);
   canvas2->SetFrameBorderMode(0);
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("");
   multigraph->SetTitle("grav_LaxWendroff_SURF_0;x;u(x,t)");
   
   Double_t Graph_fx17[500] = {
   -0.6,
   -0.598,
   -0.596,
   -0.594,
   -0.592,
   -0.59,
   -0.588,
   -0.586,
   -0.584,
   -0.582,
   -0.58,
   -0.578,
   -0.576,
   -0.574,
   -0.572,
   -0.57,
   -0.568,
   -0.566,
   -0.564,
   -0.562,
   -0.56,
   -0.558,
   -0.556,
   -0.554,
   -0.552,
   -0.55,
   -0.548,
   -0.546,
   -0.544,
   -0.542,
   -0.54,
   -0.538,
   -0.536,
   -0.534,
   -0.532,
   -0.53,
   -0.528,
   -0.526,
   -0.524,
   -0.522,
   -0.52,
   -0.518,
   -0.516,
   -0.514,
   -0.512,
   -0.51,
   -0.508,
   -0.506,
   -0.504,
   -0.502,
   -0.5,
   -0.498,
   -0.496,
   -0.494,
   -0.492,
   -0.49,
   -0.488,
   -0.486,
   -0.484,
   -0.482,
   -0.48,
   -0.478,
   -0.476,
   -0.474,
   -0.472,
   -0.47,
   -0.468,
   -0.466,
   -0.464,
   -0.462,
   -0.46,
   -0.458,
   -0.456,
   -0.454,
   -0.452,
   -0.45,
   -0.448,
   -0.446,
   -0.444,
   -0.442,
   -0.44,
   -0.438,
   -0.436,
   -0.434,
   -0.432,
   -0.43,
   -0.428,
   -0.426,
   -0.424,
   -0.422,
   -0.42,
   -0.418,
   -0.416,
   -0.414,
   -0.412,
   -0.41,
   -0.408,
   -0.406,
   -0.404,
   -0.402,
   -0.4,
   -0.398,
   -0.396,
   -0.394,
   -0.392,
   -0.39,
   -0.388,
   -0.386,
   -0.384,
   -0.382,
   -0.38,
   -0.378,
   -0.376,
   -0.374,
   -0.372,
   -0.37,
   -0.368,
   -0.366,
   -0.364,
   -0.362,
   -0.36,
   -0.358,
   -0.356,
   -0.354,
   -0.352,
   -0.35,
   -0.348,
   -0.346,
   -0.344,
   -0.342,
   -0.34,
   -0.338,
   -0.336,
   -0.334,
   -0.332,
   -0.33,
   -0.328,
   -0.326,
   -0.324,
   -0.322,
   -0.32,
   -0.318,
   -0.316,
   -0.314,
   -0.312,
   -0.31,
   -0.308,
   -0.306,
   -0.304,
   -0.302,
   -0.3,
   -0.298,
   -0.296,
   -0.294,
   -0.292,
   -0.29,
   -0.288,
   -0.286,
   -0.284,
   -0.282,
   -0.28,
   -0.278,
   -0.276,
   -0.274,
   -0.272,
   -0.27,
   -0.268,
   -0.266,
   -0.264,
   -0.262,
   -0.26,
   -0.258,
   -0.256,
   -0.254,
   -0.252,
   -0.25,
   -0.248,
   -0.246,
   -0.244,
   -0.242,
   -0.24,
   -0.238,
   -0.236,
   -0.234,
   -0.232,
   -0.23,
   -0.228,
   -0.226,
   -0.224,
   -0.222,
   -0.22,
   -0.218,
   -0.216,
   -0.214,
   -0.212,
   -0.21,
   -0.208,
   -0.206,
   -0.204,
   -0.202,
   -0.2,
   -0.198,
   -0.196,
   -0.194,
   -0.192,
   -0.19,
   -0.188,
   -0.186,
   -0.184,
   -0.182,
   -0.18,
   -0.178,
   -0.176,
   -0.174,
   -0.172,
   -0.17,
   -0.168,
   -0.166,
   -0.164,
   -0.162,
   -0.16,
   -0.158,
   -0.156,
   -0.154,
   -0.152,
   -0.15,
   -0.148,
   -0.146,
   -0.144,
   -0.142,
   -0.14,
   -0.138,
   -0.136,
   -0.134,
   -0.132,
   -0.13,
   -0.128,
   -0.126,
   -0.124,
   -0.122,
   -0.12,
   -0.118,
   -0.116,
   -0.114,
   -0.112,
   -0.11,
   -0.108,
   -0.106,
   -0.104,
   -0.102,
   -0.1,
   -0.098,
   -0.096,
   -0.094,
   -0.092,
   -0.09,
   -0.088,
   -0.086,
   -0.084,
   -0.082,
   -0.08,
   -0.078,
   -0.076,
   -0.074,
   -0.072,
   -0.07,
   -0.068,
   -0.066,
   -0.064,
   -0.062,
   -0.06,
   -0.058,
   -0.056,
   -0.054,
   -0.052,
   -0.05,
   -0.048,
   -0.046,
   -0.044,
   -0.042,
   -0.04,
   -0.038,
   -0.036,
   -0.034,
   -0.032,
   -0.03,
   -0.028,
   -0.026,
   -0.024,
   -0.022,
   -0.02,
   -0.018,
   -0.016,
   -0.014,
   -0.012,
   -0.01,
   -0.008,
   -0.006,
   -0.004,
   -0.002,
   0,
   0.002,
   0.004,
   0.006,
   0.008,
   0.01,
   0.012,
   0.014,
   0.016,
   0.018,
   0.02,
   0.022,
   0.024,
   0.026,
   0.028,
   0.03,
   0.032,
   0.034,
   0.036,
   0.038,
   0.04,
   0.042,
   0.044,
   0.046,
   0.048,
   0.05,
   0.052,
   0.054,
   0.056,
   0.058,
   0.06,
   0.062,
   0.064,
   0.066,
   0.068,
   0.07,
   0.072,
   0.074,
   0.076,
   0.078,
   0.08,
   0.082,
   0.084,
   0.086,
   0.088,
   0.09,
   0.092,
   0.094,
   0.096,
   0.098,
   0.1,
   0.102,
   0.104,
   0.106,
   0.108,
   0.11,
   0.112,
   0.114,
   0.116,
   0.118,
   0.12,
   0.122,
   0.124,
   0.126,
   0.128,
   0.13,
   0.132,
   0.134,
   0.136,
   0.138,
   0.14,
   0.142,
   0.144,
   0.146,
   0.148,
   0.15,
   0.152,
   0.154,
   0.156,
   0.158,
   0.16,
   0.162,
   0.164,
   0.166,
   0.168,
   0.17,
   0.172,
   0.174,
   0.176,
   0.178,
   0.18,
   0.182,
   0.184,
   0.186,
   0.188,
   0.19,
   0.192,
   0.194,
   0.196,
   0.198,
   0.2,
   0.202,
   0.204,
   0.206,
   0.208,
   0.21,
   0.212,
   0.214,
   0.216,
   0.218,
   0.22,
   0.222,
   0.224,
   0.226,
   0.228,
   0.23,
   0.232,
   0.234,
   0.236,
   0.238,
   0.24,
   0.242,
   0.244,
   0.246,
   0.248,
   0.25,
   0.252,
   0.254,
   0.256,
   0.258,
   0.26,
   0.262,
   0.264,
   0.266,
   0.268,
   0.27,
   0.272,
   0.274,
   0.276,
   0.278,
   0.28,
   0.282,
   0.284,
   0.286,
   0.288,
   0.29,
   0.292,
   0.294,
   0.296,
   0.298,
   0.3,
   0.302,
   0.304,
   0.306,
   0.308,
   0.31,
   0.312,
   0.314,
   0.316,
   0.318,
   0.32,
   0.322,
   0.324,
   0.326,
   0.328,
   0.33,
   0.332,
   0.334,
   0.336,
   0.338,
   0.34,
   0.342,
   0.344,
   0.346,
   0.348,
   0.35,
   0.352,
   0.354,
   0.356,
   0.358,
   0.36,
   0.362,
   0.364,
   0.366,
   0.368,
   0.37,
   0.372,
   0.374,
   0.376,
   0.378,
   0.38,
   0.382,
   0.384,
   0.386,
   0.388,
   0.39,
   0.392,
   0.394,
   0.396,
   0.398};
   Double_t Graph_fy17[500] = {
   5.380186e-32,
   8.687828e-32,
   1.400652e-31,
   2.254522e-31,
   3.623129e-31,
   5.813239e-31,
   9.312317e-31,
   1.489369e-30,
   2.378221e-30,
   3.791466e-30,
   6.034861e-30,
   9.590306e-30,
   1.521608e-29,
   2.41034e-29,
   3.812052e-29,
   6.01928e-29,
   9.489327e-29,
   1.49359e-28,
   2.347105e-28,
   3.682466e-28,
   5.76833e-28,
   9.021247e-28,
   1.408602e-27,
   2.195912e-27,
   3.417802e-27,
   5.311092e-27,
   8.239977e-27,
   1.27636e-26,
   1.973902e-26,
   3.047777e-26,
   4.698355e-26,
   7.231254e-26,
   1.111185e-25,
   1.704765e-25,
   2.611246e-25,
   3.993337e-25,
   6.097186e-25,
   9.294542e-25,
   1.414594e-24,
   2.149515e-24,
   3.261027e-24,
   4.939392e-24,
   7.469606e-24,
   1.127787e-23,
   1.700049e-23,
   2.558592e-23,
   3.844552e-23,
   5.767606e-23,
   8.638743e-23,
   1.291846e-22,
   1.92875e-22,
   2.875056e-22,
   4.278798e-22,
   6.357734e-22,
   9.431659e-22,
   1.396944e-21,
   2.065737e-21,
   3.049833e-21,
   4.495545e-21,
   6.615974e-21,
   9.720985e-21,
   1.426041e-20,
   2.088616e-20,
   3.054152e-20,
   4.458899e-20,
   6.499348e-20,
   9.458387e-20,
   1.374262e-19,
   1.99355e-19,
   2.887285e-19,
   4.17501e-19,
   6.027406e-19,
   8.687774e-19,
   1.250235e-18,
   1.796305e-18,
   2.576757e-18,
   3.690388e-18,
   5.276862e-18,
   7.533289e-18,
   1.073739e-17,
   1.52798e-17,
   2.17091e-17,
   3.079436e-17,
   4.361197e-17,
   6.166594e-17,
   8.705427e-17,
   1.226987e-16,
   1.726613e-16,
   2.4258e-16,
   3.402674e-16,
   4.765305e-16,
   6.662944e-16,
   9.301368e-16,
   1.296381e-15,
   1.803946e-15,
   2.506222e-15,
   3.476327e-15,
   4.814231e-15,
   6.656384e-15,
   9.188718e-15,
   1.266417e-14,
   1.742623e-14,
   2.394062e-14,
   3.283767e-14,
   4.496914e-14,
   6.148396e-14,
   8.392944e-14,
   1.143857e-13,
   1.556448e-13,
   2.114474e-13,
   2.867975e-13,
   3.88377e-13,
   5.250935e-13,
   7.088021e-13,
   9.552532e-13,
   1.285337e-12,
   1.726716e-12,
   2.315953e-12,
   3.101299e-12,
   4.14632e-12,
   5.53461e-12,
   7.375924e-12,
   9.81411e-12,
   1.303739e-11,
   1.72916e-11,
   2.289735e-11,
   3.027194e-11,
   3.995769e-11,
   5.265815e-11,
   6.928449e-11,
   9.101471e-11,
   1.193692e-10,
   1.563069e-10,
   2.043473e-10,
   2.667258e-10,
   3.475891e-10,
   4.522437e-10,
   5.874676e-10,
   7.619044e-10,
   9.86557e-10,
   1.275408e-09,
   1.646194e-09,
   2.121378e-09,
   2.729356e-09,
   3.505965e-09,
   4.496349e-09,
   5.757284e-09,
   7.360043e-09,
   9.39395e-09,
   1.197075e-08,
   1.522998e-08,
   1.934562e-08,
   2.453415e-08,
   3.10645e-08,
   3.92702e-08,
   4.956405e-08,
   6.245622e-08,
   7.857596e-08,
   9.869811e-08,
   1.237751e-07,
   1.549753e-07,
   1.937301e-07,
   2.417891e-07,
   3.012878e-07,
   3.748275e-07,
   4.655716e-07,
   5.773599e-07,
   7.148451e-07,
   8.836542e-07,
   1.090581e-06,
   1.343812e-06,
   1.653196e-06,
   2.030558e-06,
   2.490069e-06,
   3.048686e-06,
   3.726653e-06,
   4.548105e-06,
   5.541751e-06,
   6.741689e-06,
   8.188335e-06,
   9.929504e-06,
   1.202167e-05,
   1.453138e-05,
   1.753696e-05,
   2.113035e-05,
   2.541935e-05,
   3.053002e-05,
   3.66096e-05,
   4.382965e-05,
   5.238973e-05,
   6.25215e-05,
   7.44934e-05,
   8.861583e-05,
   0.0001052471,
   0.0001247998,
   0.0001477484,
   0.0001746372,
   0.0002060895,
   0.0002428176,
   0.0002856338,
   0.0003354626,
   0.0003933542,
   0.0004604989,
   0.0005382432,
   0.000628107,
   0.0007318024,
   0.000851254,
   0.0009886205,
   0.001146318,
   0.001327046,
   0.001533811,
   0.001769957,
   0.002039195,
   0.002345633,
   0.002693807,
   0.003088715,
   0.003535856,
   0.004041255,
   0.00461151,
   0.005253819,
   0.005976023,
   0.006786635,
   0.007694881,
   0.008710727,
   0.009844917,
   0.011109,
   0.01251534,
   0.01407718,
   0.01580862,
   0.01772463,
   0.01984109,
   0.02217477,
   0.02474331,
   0.02756523,
   0.03065989,
   0.03404745,
   0.03774886,
   0.04178575,
   0.04618039,
   0.05095563,
   0.05613476,
   0.06174144,
   0.06779953,
   0.07433302,
   0.08136582,
   0.08892162,
   0.0970237,
   0.1056948,
   0.1149567,
   0.1248303,
   0.1353353,
   0.1464897,
   0.15831,
   0.1708106,
   0.1840036,
   0.1978987,
   0.2125028,
   0.2278199,
   0.2438505,
   0.2605918,
   0.2780373,
   0.2961764,
   0.3149945,
   0.3344727,
   0.3545875,
   0.3753111,
   0.3966107,
   0.4184491,
   0.4407841,
   0.463569,
   0.4867523,
   0.5102778,
   0.5340851,
   0.5581096,
   0.5822822,
   0.6065307,
   0.6307788,
   0.6549476,
   0.6789553,
   0.7027177,
   0.726149,
   0.749162,
   0.7716687,
   0.7935807,
   0.8148103,
   0.8352702,
   0.854875,
   0.8735412,
   0.8911879,
   0.9077375,
   0.9231163,
   0.9372549,
   0.9500886,
   0.9615584,
   0.9716108,
   0.9801987,
   0.9872816,
   0.9928259,
   0.9968051,
   0.9992003,
   1,
   0.9992003,
   0.9968051,
   0.9928259,
   0.9872816,
   0.9801987,
   0.9716108,
   0.9615584,
   0.9500886,
   0.9372549,
   0.9231163,
   0.9077375,
   0.8911879,
   0.8735412,
   0.854875,
   0.8352702,
   0.8148103,
   0.7935807,
   0.7716687,
   0.749162,
   0.726149,
   0.7027177,
   0.6789553,
   0.6549476,
   0.6307788,
   0.6065307,
   0.5822822,
   0.5581096,
   0.5340851,
   0.5102778,
   0.4867523,
   0.463569,
   0.4407841,
   0.4184491,
   0.3966107,
   0.3753111,
   0.3545875,
   0.3344727,
   0.3149945,
   0.2961764,
   0.2780373,
   0.2605918,
   0.2438505,
   0.2278199,
   0.2125028,
   0.1978987,
   0.1840036,
   0.1708106,
   0.15831,
   0.1464897,
   0.1353353,
   0.1248303,
   0.1149567,
   0.1056948,
   0.0970237,
   0.08892162,
   0.08136582,
   0.07433302,
   0.06779953,
   0.06174144,
   0.05613476,
   0.05095563,
   0.04618039,
   0.04178575,
   0.03774886,
   0.03404745,
   0.03065989,
   0.02756523,
   0.02474331,
   0.02217477,
   0.01984109,
   0.01772463,
   0.01580862,
   0.01407718,
   0.01251534,
   0.011109,
   0.009844917,
   0.008710727,
   0.007694881,
   0.006786635,
   0.005976023,
   0.005253819,
   0.00461151,
   0.004041255,
   0.003535856,
   0.003088715,
   0.002693807,
   0.002345633,
   0.002039195,
   0.001769957,
   0.001533811,
   0.001327046,
   0.001146318,
   0.0009886205,
   0.000851254,
   0.0007318024,
   0.000628107,
   0.0005382432,
   0.0004604989,
   0.0003933542,
   0.0003354626,
   0.0002856338,
   0.0002428176,
   0.0002060895,
   0.0001746372,
   0.0001477484,
   0.0001247998,
   0.0001052471,
   8.861583e-05,
   7.44934e-05,
   6.25215e-05,
   5.238973e-05,
   4.382965e-05,
   3.66096e-05,
   3.053002e-05,
   2.541935e-05,
   2.113035e-05,
   1.753696e-05,
   1.453138e-05,
   1.202167e-05,
   9.929504e-06,
   8.188335e-06,
   6.741689e-06,
   5.541751e-06,
   4.548105e-06,
   3.726653e-06,
   3.048686e-06,
   2.490069e-06,
   2.030558e-06,
   1.653196e-06,
   1.343812e-06,
   1.090581e-06,
   8.836542e-07,
   7.148451e-07,
   5.773599e-07,
   4.655716e-07,
   3.748275e-07,
   3.012878e-07,
   2.417891e-07,
   1.937301e-07,
   1.549753e-07,
   1.237751e-07,
   9.869811e-08,
   7.857596e-08,
   6.245622e-08,
   4.956405e-08,
   3.92702e-08,
   3.10645e-08,
   2.453415e-08,
   1.934562e-08,
   1.522998e-08,
   1.197075e-08,
   9.39395e-09,
   7.360043e-09,
   5.757284e-09,
   4.496349e-09,
   3.505965e-09,
   2.729356e-09,
   2.121378e-09,
   1.646194e-09,
   1.275408e-09,
   9.86557e-10,
   7.619044e-10,
   5.874676e-10,
   4.522437e-10,
   3.475891e-10,
   2.667258e-10,
   2.043473e-10,
   1.563069e-10,
   1.193692e-10,
   9.101471e-11,
   6.928449e-11,
   5.265815e-11,
   3.995769e-11,
   3.027194e-11,
   2.289735e-11,
   1.72916e-11,
   1.303739e-11,
   9.81411e-12,
   7.375924e-12,
   5.53461e-12,
   4.14632e-12,
   3.101299e-12,
   2.315953e-12,
   1.726716e-12,
   1.285337e-12,
   9.552532e-13,
   7.088021e-13,
   5.250935e-13,
   3.88377e-13,
   2.867975e-13,
   2.114474e-13,
   1.556448e-13,
   1.143857e-13,
   8.392944e-14,
   6.148396e-14,
   4.496914e-14,
   3.283767e-14,
   2.394062e-14,
   1.742623e-14};
   TGraph *graph = new TGraph(500,Graph_fx17,Graph_fy17);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   
   TH1F *Graph_Graph17 = new TH1F("Graph_Graph17","Graph",500,-0.6998,0.4978);
   Graph_Graph17->SetMinimum(4.842168e-32);
   Graph_Graph17->SetMaximum(1.1);
   Graph_Graph17->SetDirectory(0);
   Graph_Graph17->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph17->SetLineColor(ci);
   Graph_Graph17->GetXaxis()->SetLabelFont(42);
   Graph_Graph17->GetXaxis()->SetTitleOffset(1);
   Graph_Graph17->GetXaxis()->SetTitleFont(42);
   Graph_Graph17->GetYaxis()->SetLabelFont(42);
   Graph_Graph17->GetYaxis()->SetTitleFont(42);
   Graph_Graph17->GetZaxis()->SetLabelFont(42);
   Graph_Graph17->GetZaxis()->SetTitleOffset(1);
   Graph_Graph17->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph17);
   
   multigraph->Add(graph,"");
   
   Double_t Graph_fx18[500] = {
   -0.6,
   -0.598,
   -0.596,
   -0.594,
   -0.592,
   -0.59,
   -0.588,
   -0.586,
   -0.584,
   -0.582,
   -0.58,
   -0.578,
   -0.576,
   -0.574,
   -0.572,
   -0.57,
   -0.568,
   -0.566,
   -0.564,
   -0.562,
   -0.56,
   -0.558,
   -0.556,
   -0.554,
   -0.552,
   -0.55,
   -0.548,
   -0.546,
   -0.544,
   -0.542,
   -0.54,
   -0.538,
   -0.536,
   -0.534,
   -0.532,
   -0.53,
   -0.528,
   -0.526,
   -0.524,
   -0.522,
   -0.52,
   -0.518,
   -0.516,
   -0.514,
   -0.512,
   -0.51,
   -0.508,
   -0.506,
   -0.504,
   -0.502,
   -0.5,
   -0.498,
   -0.496,
   -0.494,
   -0.492,
   -0.49,
   -0.488,
   -0.486,
   -0.484,
   -0.482,
   -0.48,
   -0.478,
   -0.476,
   -0.474,
   -0.472,
   -0.47,
   -0.468,
   -0.466,
   -0.464,
   -0.462,
   -0.46,
   -0.458,
   -0.456,
   -0.454,
   -0.452,
   -0.45,
   -0.448,
   -0.446,
   -0.444,
   -0.442,
   -0.44,
   -0.438,
   -0.436,
   -0.434,
   -0.432,
   -0.43,
   -0.428,
   -0.426,
   -0.424,
   -0.422,
   -0.42,
   -0.418,
   -0.416,
   -0.414,
   -0.412,
   -0.41,
   -0.408,
   -0.406,
   -0.404,
   -0.402,
   -0.4,
   -0.398,
   -0.396,
   -0.394,
   -0.392,
   -0.39,
   -0.388,
   -0.386,
   -0.384,
   -0.382,
   -0.38,
   -0.378,
   -0.376,
   -0.374,
   -0.372,
   -0.37,
   -0.368,
   -0.366,
   -0.364,
   -0.362,
   -0.36,
   -0.358,
   -0.356,
   -0.354,
   -0.352,
   -0.35,
   -0.348,
   -0.346,
   -0.344,
   -0.342,
   -0.34,
   -0.338,
   -0.336,
   -0.334,
   -0.332,
   -0.33,
   -0.328,
   -0.326,
   -0.324,
   -0.322,
   -0.32,
   -0.318,
   -0.316,
   -0.314,
   -0.312,
   -0.31,
   -0.308,
   -0.306,
   -0.304,
   -0.302,
   -0.3,
   -0.298,
   -0.296,
   -0.294,
   -0.292,
   -0.29,
   -0.288,
   -0.286,
   -0.284,
   -0.282,
   -0.28,
   -0.278,
   -0.276,
   -0.274,
   -0.272,
   -0.27,
   -0.268,
   -0.266,
   -0.264,
   -0.262,
   -0.26,
   -0.258,
   -0.256,
   -0.254,
   -0.252,
   -0.25,
   -0.248,
   -0.246,
   -0.244,
   -0.242,
   -0.24,
   -0.238,
   -0.236,
   -0.234,
   -0.232,
   -0.23,
   -0.228,
   -0.226,
   -0.224,
   -0.222,
   -0.22,
   -0.218,
   -0.216,
   -0.214,
   -0.212,
   -0.21,
   -0.208,
   -0.206,
   -0.204,
   -0.202,
   -0.2,
   -0.198,
   -0.196,
   -0.194,
   -0.192,
   -0.19,
   -0.188,
   -0.186,
   -0.184,
   -0.182,
   -0.18,
   -0.178,
   -0.176,
   -0.174,
   -0.172,
   -0.17,
   -0.168,
   -0.166,
   -0.164,
   -0.162,
   -0.16,
   -0.158,
   -0.156,
   -0.154,
   -0.152,
   -0.15,
   -0.148,
   -0.146,
   -0.144,
   -0.142,
   -0.14,
   -0.138,
   -0.136,
   -0.134,
   -0.132,
   -0.13,
   -0.128,
   -0.126,
   -0.124,
   -0.122,
   -0.12,
   -0.118,
   -0.116,
   -0.114,
   -0.112,
   -0.11,
   -0.108,
   -0.106,
   -0.104,
   -0.102,
   -0.1,
   -0.098,
   -0.096,
   -0.094,
   -0.092,
   -0.09,
   -0.088,
   -0.086,
   -0.084,
   -0.082,
   -0.08,
   -0.078,
   -0.076,
   -0.074,
   -0.072,
   -0.07,
   -0.068,
   -0.066,
   -0.064,
   -0.062,
   -0.06,
   -0.058,
   -0.056,
   -0.054,
   -0.052,
   -0.05,
   -0.048,
   -0.046,
   -0.044,
   -0.042,
   -0.04,
   -0.038,
   -0.036,
   -0.034,
   -0.032,
   -0.03,
   -0.028,
   -0.026,
   -0.024,
   -0.022,
   -0.02,
   -0.018,
   -0.016,
   -0.014,
   -0.012,
   -0.01,
   -0.008,
   -0.006,
   -0.004,
   -0.002,
   0,
   0.002,
   0.004,
   0.006,
   0.008,
   0.01,
   0.012,
   0.014,
   0.016,
   0.018,
   0.02,
   0.022,
   0.024,
   0.026,
   0.028,
   0.03,
   0.032,
   0.034,
   0.036,
   0.038,
   0.04,
   0.042,
   0.044,
   0.046,
   0.048,
   0.05,
   0.052,
   0.054,
   0.056,
   0.058,
   0.06,
   0.062,
   0.064,
   0.066,
   0.068,
   0.07,
   0.072,
   0.074,
   0.076,
   0.078,
   0.08,
   0.082,
   0.084,
   0.086,
   0.088,
   0.09,
   0.092,
   0.094,
   0.096,
   0.098,
   0.1,
   0.102,
   0.104,
   0.106,
   0.108,
   0.11,
   0.112,
   0.114,
   0.116,
   0.118,
   0.12,
   0.122,
   0.124,
   0.126,
   0.128,
   0.13,
   0.132,
   0.134,
   0.136,
   0.138,
   0.14,
   0.142,
   0.144,
   0.146,
   0.148,
   0.15,
   0.152,
   0.154,
   0.156,
   0.158,
   0.16,
   0.162,
   0.164,
   0.166,
   0.168,
   0.17,
   0.172,
   0.174,
   0.176,
   0.178,
   0.18,
   0.182,
   0.184,
   0.186,
   0.188,
   0.19,
   0.192,
   0.194,
   0.196,
   0.198,
   0.2,
   0.202,
   0.204,
   0.206,
   0.208,
   0.21,
   0.212,
   0.214,
   0.216,
   0.218,
   0.22,
   0.222,
   0.224,
   0.226,
   0.228,
   0.23,
   0.232,
   0.234,
   0.236,
   0.238,
   0.24,
   0.242,
   0.244,
   0.246,
   0.248,
   0.25,
   0.252,
   0.254,
   0.256,
   0.258,
   0.26,
   0.262,
   0.264,
   0.266,
   0.268,
   0.27,
   0.272,
   0.274,
   0.276,
   0.278,
   0.28,
   0.282,
   0.284,
   0.286,
   0.288,
   0.29,
   0.292,
   0.294,
   0.296,
   0.298,
   0.3,
   0.302,
   0.304,
   0.306,
   0.308,
   0.31,
   0.312,
   0.314,
   0.316,
   0.318,
   0.32,
   0.322,
   0.324,
   0.326,
   0.328,
   0.33,
   0.332,
   0.334,
   0.336,
   0.338,
   0.34,
   0.342,
   0.344,
   0.346,
   0.348,
   0.35,
   0.352,
   0.354,
   0.356,
   0.358,
   0.36,
   0.362,
   0.364,
   0.366,
   0.368,
   0.37,
   0.372,
   0.374,
   0.376,
   0.378,
   0.38,
   0.382,
   0.384,
   0.386,
   0.388,
   0.39,
   0.392,
   0.394,
   0.396,
   0.398};
   Double_t Graph_fy18[500] = {
   -0.01877715,
   -0.01877548,
   -0.01877379,
   -0.01877205,
   -0.01877028,
   -0.01876847,
   -0.01876663,
   -0.01876476,
   -0.01876287,
   -0.01876095,
   -0.01875902,
   -0.01875706,
   -0.01875508,
   -0.01875309,
   -0.01875109,
   -0.01874907,
   -0.01874705,
   -0.01874501,
   -0.01874297,
   -0.01874093,
   -0.01873888,
   -0.01873683,
   -0.01873478,
   -0.01873274,
   -0.01873069,
   -0.01872865,
   -0.01872661,
   -0.01872458,
   -0.01872255,
   -0.01872054,
   -0.01871853,
   -0.01871653,
   -0.01871455,
   -0.01871258,
   -0.01871062,
   -0.01870867,
   -0.01870674,
   -0.01870482,
   -0.01870292,
   -0.01870104,
   -0.01869918,
   -0.01869733,
   -0.0186955,
   -0.0186937,
   -0.01869191,
   -0.01869014,
   -0.01868839,
   -0.01868667,
   -0.01868497,
   -0.01868329,
   -0.01868163,
   -0.01868,
   -0.01867839,
   -0.0186768,
   -0.01867524,
   -0.01867371,
   -0.01867219,
   -0.01867071,
   -0.01866925,
   -0.01866782,
   -0.01866641,
   -0.01866503,
   -0.01866368,
   -0.01866235,
   -0.01866105,
   -0.01865978,
   -0.01865853,
   -0.01865732,
   -0.01865613,
   -0.01865497,
   -0.01865384,
   -0.01865273,
   -0.01865165,
   -0.01865061,
   -0.01864959,
   -0.0186486,
   -0.01864764,
   -0.0186467,
   -0.0186458,
   -0.01864492,
   -0.01864408,
   -0.01864326,
   -0.01864247,
   -0.01864171,
   -0.01864098,
   -0.01864028,
   -0.01863961,
   -0.01863897,
   -0.01863835,
   -0.01863777,
   -0.01863721,
   -0.01863668,
   -0.01863618,
   -0.01863571,
   -0.01863527,
   -0.01863486,
   -0.01863448,
   -0.01863412,
   -0.01863379,
   -0.0186335,
   -0.01863323,
   -0.01863298,
   -0.01863277,
   -0.01863259,
   -0.01863243,
   -0.0186323,
   -0.0186322,
   -0.01863213,
   -0.01863208,
   -0.01863207,
   -0.01863208,
   -0.01863212,
   -0.01863218,
   -0.01863228,
   -0.0186324,
   -0.01863254,
   -0.01863272,
   -0.01863292,
   -0.01863315,
   -0.01863341,
   -0.01863369,
   -0.018634,
   -0.01863433,
   -0.0186347,
   -0.01863508,
   -0.0186355,
   -0.01863594,
   -0.01863641,
   -0.0186369,
   -0.01863742,
   -0.01863797,
   -0.01863854,
   -0.01863914,
   -0.01863976,
   -0.01864041,
   -0.01864108,
   -0.01864178,
   -0.01864251,
   -0.01864326,
   -0.01864404,
   -0.01864484,
   -0.01864567,
   -0.01864653,
   -0.01864741,
   -0.01864831,
   -0.01864924,
   -0.0186502,
   -0.01865118,
   -0.01865219,
   -0.01865322,
   -0.01865428,
   -0.01865536,
   -0.01865647,
   -0.0186576,
   -0.01865876,
   -0.01865994,
   -0.01866114,
   -0.01866237,
   -0.01866362,
   -0.01866489,
   -0.01866618,
   -0.01866748,
   -0.01866881,
   -0.01867015,
   -0.01867149,
   -0.01867285,
   -0.01867421,
   -0.01867557,
   -0.01867692,
   -0.01867826,
   -0.01867957,
   -0.01868085,
   -0.01868209,
   -0.01868326,
   -0.01868435,
   -0.01868535,
   -0.01868621,
   -0.01868692,
   -0.01868744,
   -0.01868773,
   -0.01868773,
   -0.01868739,
   -0.01868664,
   -0.0186854,
   -0.01868357,
   -0.01868104,
   -0.01867768,
   -0.01867335,
   -0.01866785,
   -0.018661,
   -0.01865255,
   -0.01864222,
   -0.01862971,
   -0.01861464,
   -0.01859659,
   -0.01857508,
   -0.01854955,
   -0.01851937,
   -0.01848382,
   -0.01844206,
   -0.01839316,
   -0.01833605,
   -0.01826953,
   -0.01819223,
   -0.01810261,
   -0.01799893,
   -0.01787924,
   -0.01774138,
   -0.01758288,
   -0.01740102,
   -0.01719277,
   -0.01695475,
   -0.01668321,
   -0.01637402,
   -0.01602261,
   -0.01562392,
   -0.01517245,
   -0.01466211,
   -0.01408627,
   -0.01343769,
   -0.0127085,
   -0.01189015,
   -0.01097336,
   -0.00994815,
   -0.00880373,
   -0.007528527,
   -0.00611014,
   -0.004535315,
   -0.002789938,
   -0.0008590117,
   0.001273341,
   0.003623883,
   0.006210253,
   0.009050947,
   0.0121653,
   0.01557343,
   0.01929624,
   0.02335529,
   0.0277728,
   0.03257151,
   0.03777461,
   0.04340564,
   0.04948835,
   0.05604658,
   0.06310407,
   0.07068436,
   0.07881053,
   0.08750508,
   0.09678968,
   0.106685,
   0.1172104,
   0.1283838,
   0.1402212,
   0.152737,
   0.165943,
   0.1798486,
   0.1944605,
   0.2097823,
   0.2258145,
   0.242554,
   0.259994,
   0.2781238,
   0.2969286,
   0.3163893,
   0.3364824,
   0.3571798,
   0.3784489,
   0.4002523,
   0.4225478,
   0.4452888,
   0.4684238,
   0.4918969,
   0.5156477,
   0.5396118,
   0.5637206,
   0.5879018,
   0.6120797,
   0.6361757,
   0.6601084,
   0.683794,
   0.7071472,
   0.7300813,
   0.7525087,
   0.7743417,
   0.7954929,
   0.8158759,
   0.8354055,
   0.853999,
   0.8715759,
   0.8880592,
   0.9033757,
   0.9174563,
   0.9302371,
   0.9416592,
   0.9516696,
   0.9602217,
   0.967275,
   0.9727963,
   0.9767593,
   0.9791452,
   0.9799425,
   0.9791474,
   0.9767637,
   0.9728028,
   0.9672835,
   0.9602321,
   0.9516819,
   0.9416731,
   0.9302524,
   0.9174729,
   0.9033933,
   0.8880776,
   0.8715948,
   0.8540182,
   0.8354247,
   0.8158948,
   0.7955113,
   0.7743592,
   0.7525251,
   0.7300963,
   0.7071606,
   0.6838055,
   0.6601176,
   0.6361826,
   0.6120839,
   0.5879031,
   0.5637188,
   0.5396068,
   0.5156393,
   0.4918849,
   0.4684082,
   0.4452694,
   0.4225246,
   0.4002251,
   0.3784178,
   0.3571448,
   0.3364433,
   0.3163463,
   0.2968816,
   0.2780729,
   0.2599392,
   0.2424954,
   0.2257522,
   0.2097164,
   0.194391,
   0.1797757,
   0.1658667,
   0.1526575,
   0.1401386,
   0.1282981,
   0.1171219,
   0.1065937,
   0.09669579,
   0.08740866,
   0.0787117,
   0.07058322,
   0.06300074,
   0.05594115,
   0.04938093,
   0.04329631,
   0.03766345,
   0.03245861,
   0.02765824,
   0.02323914,
   0.01917855,
   0.01545427,
   0.01204471,
   0.008928984,
   0.006086956,
   0.003499289,
   0.001147481,
   -0.0009861102,
   -0.002918251,
   -0.004664825,
   -0.006240829,
   -0.007660385,
   -0.008936746,
   -0.01008232,
   -0.01110868,
   -0.01202661,
   -0.01284611,
   -0.01357645,
   -0.01422617,
   -0.01480316,
   -0.01531466,
   -0.0157673,
   -0.01616715,
   -0.01651974,
   -0.0168301,
   -0.01710282,
   -0.01734203,
   -0.01755148,
   -0.01773453,
   -0.01789424,
   -0.01803331,
   -0.0181542,
   -0.01825909,
   -0.01834992,
   -0.01842843,
   -0.01849615,
   -0.01855446,
   -0.01860456,
   -0.01864751,
   -0.01868424,
   -0.0187156,
   -0.01874229,
   -0.01876496,
   -0.01878415,
   -0.01880035,
   -0.01881398,
   -0.01882541,
   -0.01883494,
   -0.01884286,
   -0.0188494,
   -0.01885476,
   -0.01885912,
   -0.01886263,
   -0.01886541,
   -0.01886759,
   -0.01886925,
   -0.01887047,
   -0.01887132,
   -0.01887186,
   -0.01887214,
   -0.01887219,
   -0.01887206,
   -0.01887178,
   -0.01887136,
   -0.01887083,
   -0.01887021,
   -0.01886951,
   -0.01886875,
   -0.01886792,
   -0.01886705,
   -0.01886614,
   -0.01886519,
   -0.01886421,
   -0.01886321,
   -0.01886218,
   -0.01886113,
   -0.01886007,
   -0.01885899,
   -0.01885789,
   -0.01885677,
   -0.01885565,
   -0.01885451,
   -0.01885336,
   -0.0188522,
   -0.01885103,
   -0.01884985,
   -0.01884866,
   -0.01884746,
   -0.01884625,
   -0.01884503,
   -0.0188438,
   -0.01884256,
   -0.01884131,
   -0.01884005,
   -0.01883878,
   -0.01883751,
   -0.01883622,
   -0.01883493,
   -0.01883362,
   -0.01883231,
   -0.01883099,
   -0.01882965,
   -0.01882831,
   -0.01882696,
   -0.0188256,
   -0.01882423,
   -0.01882286,
   -0.01882147,
   -0.01882007,
   -0.01881867,
   -0.01881725,
   -0.01881583,
   -0.01881439,
   -0.01881295,
   -0.0188115,
   -0.01881003,
   -0.01880856,
   -0.01880708,
   -0.01880559,
   -0.01880409,
   -0.01880258,
   -0.01880106,
   -0.01879953,
   -0.018798,
   -0.01879645,
   -0.01879489,
   -0.01879333,
   -0.01879175,
   -0.01879016,
   -0.01878857,
   -0.01878697,
   -0.01878535,
   -0.01878373,
   -0.0187821,
   -0.01878045,
   -0.0187788};
   graph = new TGraph(500,Graph_fx18,Graph_fy18);
   graph->SetName("Graph");
   graph->SetTitle("Graph");
   graph->SetFillStyle(1000);
   graph->SetLineColor(2);
   graph->SetMarkerColor(2);
   
   TH1F *Graph_Graph18 = new TH1F("Graph_Graph18","Graph",500,-0.6998,0.4978);
   Graph_Graph18->SetMinimum(-0.1187537);
   Graph_Graph18->SetMaximum(1.079824);
   Graph_Graph18->SetDirectory(0);
   Graph_Graph18->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph18->SetLineColor(ci);
   Graph_Graph18->GetXaxis()->SetLabelFont(42);
   Graph_Graph18->GetXaxis()->SetTitleOffset(1);
   Graph_Graph18->GetXaxis()->SetTitleFont(42);
   Graph_Graph18->GetYaxis()->SetLabelFont(42);
   Graph_Graph18->GetYaxis()->SetTitleFont(42);
   Graph_Graph18->GetZaxis()->SetLabelFont(42);
   Graph_Graph18->GetZaxis()->SetTitleOffset(1);
   Graph_Graph18->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph18);
   
   multigraph->Add(graph,"");
   multigraph->Draw("APL");
   multigraph->GetXaxis()->SetLimits(-0.6499, 0.4479);
   multigraph->GetXaxis()->SetTitle("x");
   multigraph->GetXaxis()->SetLabelFont(42);
   multigraph->GetXaxis()->SetTitleOffset(1);
   multigraph->GetXaxis()->SetTitleFont(42);
   multigraph->GetYaxis()->SetTitle("u(x,t)");
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
   
   TPaveText *pt = new TPaveText(0.2802685,0.9304839,0.7197315,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("grav_LaxWendroff_SURF_0");
   pt->Draw();
   canvas2->Modified();
   canvas2->cd();
   canvas2->SetSelected(canvas2);
}
