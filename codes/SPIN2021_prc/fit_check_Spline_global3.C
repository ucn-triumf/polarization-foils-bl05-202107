#include "../bin/MakeNiki.h"
//#include "/home/nop/data/Tools/GetMRcut.C"
#include <TH3.h>
TString path_R = "results/";

// Double_t Distance = 20.000; //tentative
// Double_t Conversion = 395.6;
// Double_t dist_det   = 337.; //sample to detector [mm]
// Double_t xdirect    = 64.71;
// Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
// Bool_t useThinout = 0; //thinning out the event <1e4.

// definition of shared parameter
// background function
int iparB[3] = { 0,      // exp amplitude in B histo
                 1,
                 2    // exp common parameter
};
 
// signal + background function
int iparSB[3] = { 0,
                  1, // exp amplitude in S+B histo
                  2 // exp common parameter
                  
};
 
// Create the GlobalCHi2 structure
 
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f11,
                ROOT::Math::IMultiGenFunction & f22) :
      fChi2_1(&f11), fChi2_2(&f22) {}
 
   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[3];
      for (int i = 0; i < 3; ++i) p1[i] = par[iparB[i] ];
 
      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparSB[i] ];
 
      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }
 
   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};
//8/29 add
const double c=299792458; //[m/s]
const double c_nm=299792458.e9; //[nm/s]
const double m1=939.5654133; //[MeV/c^2]
const double m2=939.5654133e15; //[neV/c^2]
//Double_t m=939.5654133e9/(pow(c,2)); //[meV/(m/s)^2]
double m_nc2=939.5654133e15/(pow(c,2)); //[neV/(m/s)^2]
double m_nc2nm=m2/(pow(c_nm,2)); //[neV/(nm/s)^2]


const double h1=4.135667696e-15; //[eV.s]
const double h=4.135667696e-12; //[meV.s]
const double hbar=4.135667696e-6/(2.*TMath::Pi());//neV.s

const double V_Fe1=209.0602;//neV
const double V_Si=54.0078;//neV
const double mu_n=60.3;//[neV T^-1] 
const double mu_n2=-60.3;//[neV T^-1] 

const double in_T=2.;//T(tesla)
const double V_Fe_p1=V_Fe1+mu_n*in_T;//neV
const double V_Fe_m1=V_Fe1-mu_n*in_T;//neV
  /*double k1b=sqrt(2*m_nc2*(E1-(V_Fe-mu_n*in_T)))/hbar;//m^-1
  double k2=sqrt(2*m_nc2*(E1-(V_Si))/hbar;//m^-1
  double alpha1=k1b=sqrt(2*m_nc2*((V_Fe-mu_n*in_T)-E1))/hbar;
  double alpha2=k1b=sqrt(2*m_nc2*((V_Fe+mu_n*in_T)-E1))/hbar;
*/
  //E>V1(+V)
//double qc=0.13;//par[2];
//double ww=7.16389E-02;//par[3];
//double R0=1.;//par[4];

//double mm=par[3];
//double alpha=-8.03288E-02;//par[5];

double uprate=0.5;//par[3];
double downrate=0.5;//par[4];

double mm=5.;
double qc_up=0.22;
//double mm2=6.24502;

//double qc=1.88418e-1;
double qc=0.154971;

//double qc=1.58418e-1;
double mm2=6.72226;
//double mm2=16.72226;

double q_coff_samp=2*TMath::Pi()* ((60.4483-45.1074)/666.);
double q_coff=2*TMath::Pi()* ((89.3295-65.4023)/1439.);//((89.3295-60.4483)/1439.);
double q_coff_true=2*TMath::Pi()* ((89.3295-60.4583)/1439.);

double Spl1(double x1) {

    double x2=q_coff_samp/x1;
    
    double x=q_coff/q_coff_true*x2;
    //double x=q_coff_true/x1;
    
    const int fNp = 300, fKstep = 0;
    const double fDelta = -1, fXmin = 0.002, fXmax = 1.198;
    const double fX[300] = { 0.002, 0.006, 0.01, 0.014, 0.018,
                          0.022, 0.026, 0.03, 0.034, 0.038,
                          0.042, 0.046, 0.05, 0.054, 0.058,
                          0.062, 0.066, 0.07, 0.074, 0.078,
                          0.082, 0.086, 0.09, 0.094, 0.098,
                          0.102, 0.106, 0.11, 0.114, 0.118,
                          0.122, 0.126, 0.13, 0.134, 0.138,
                          0.142, 0.146, 0.15, 0.154, 0.158,
                          0.162, 0.166, 0.17, 0.174, 0.178,
                          0.182, 0.186, 0.19, 0.194, 0.198,
                          0.202, 0.206, 0.21, 0.214, 0.218,
                          0.222, 0.226, 0.23, 0.234, 0.238,
                          0.242, 0.246, 0.25, 0.254, 0.258,
                          0.262, 0.266, 0.27, 0.274, 0.278,
                          0.282, 0.286, 0.29, 0.294, 0.298,
                          0.302, 0.306, 0.31, 0.314, 0.318,
                          0.322, 0.326, 0.33, 0.334, 0.338,
                          0.342, 0.346, 0.35, 0.354, 0.358,
                          0.362, 0.366, 0.37, 0.374, 0.378,
                          0.382, 0.386, 0.39, 0.394, 0.398,
                          0.402, 0.406, 0.41, 0.414, 0.418,
                          0.422, 0.426, 0.43, 0.434, 0.438,
                          0.442, 0.446, 0.45, 0.454, 0.458,
                          0.462, 0.466, 0.47, 0.474, 0.478,
                          0.482, 0.486, 0.49, 0.494, 0.498,
                          0.502, 0.506, 0.51, 0.514, 0.518,
                          0.522, 0.526, 0.53, 0.534, 0.538,
                          0.542, 0.546, 0.55, 0.554, 0.558,
                          0.562, 0.566, 0.57, 0.574, 0.578,
                          0.582, 0.586, 0.59, 0.594, 0.598,
                          0.602, 0.606, 0.61, 0.614, 0.618,
                          0.622, 0.626, 0.63, 0.634, 0.638,
                          0.642, 0.646, 0.65, 0.654, 0.658,
                          0.662, 0.666, 0.67, 0.674, 0.678,
                          0.682, 0.686, 0.69, 0.694, 0.698,
                          0.702, 0.706, 0.71, 0.714, 0.718,
                          0.722, 0.726, 0.73, 0.734, 0.738,
                          0.742, 0.746, 0.75, 0.754, 0.758,
                          0.762, 0.766, 0.77, 0.774, 0.778,
                          0.782, 0.786, 0.79, 0.794, 0.798,
                          0.802, 0.806, 0.81, 0.814, 0.818,
                          0.822, 0.826, 0.83, 0.834, 0.838,
                          0.842, 0.846, 0.85, 0.854, 0.858,
                          0.862, 0.866, 0.87, 0.874, 0.878,
                          0.882, 0.886, 0.89, 0.894, 0.898,
                          0.902, 0.906, 0.91, 0.914, 0.918,
                          0.922, 0.926, 0.93, 0.934, 0.938,
                          0.942, 0.946, 0.95, 0.954, 0.958,
                          0.962, 0.966, 0.97, 0.974, 0.978,
                          0.982, 0.986, 0.99, 0.994, 0.998,
                          1.002, 1.006, 1.01, 1.014, 1.018,
                          1.022, 1.026, 1.03, 1.034, 1.038,
                          1.042, 1.046, 1.05, 1.054, 1.058,
                          1.062, 1.066, 1.07, 1.074, 1.078,
                          1.082, 1.086, 1.09, 1.094, 1.098,
                          1.102, 1.106, 1.11, 1.114, 1.118,
                          1.122, 1.126, 1.13, 1.134, 1.138,
                          1.142, 1.146, 1.15, 1.154, 1.158,
                          1.162, 1.166, 1.17, 1.174, 1.178,
                          1.182, 1.186, 1.19, 1.194, 1.198 };
    const double fY[300] = { 1, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0,
                          0, 0.826772, 0.697917, 0.650862, 0.537525,
                          0.527491, 0.502817, 0.477046, 0.477663, 0.471612,
                          0.468058, 0.490673, 0.49523, 0.498246, 0.489787,
                          0.491766, 0.496617, 0.498853, 0.493772, 0.511259,
                          0.531188, 0.529615, 0.52608, 0.503353, 0.517948,
                          0.511434, 0.516854, 0.531556, 0.528384, 0.533591,
                          0.533623, 0.536679, 0.530889, 0.526723, 0.544358,
                          0.550417, 0.541349, 0.53252, 0.537405, 0.533762,
                          0.545372, 0.541603, 0.538317, 0.538556, 0.544574,
                          0.553247, 0.546539, 0.53217, 0.550793, 0.554797,
                          0.545386, 0.540198, 0.550931, 0.545729, 0.557383,
                          0.534095, 0.558067, 0.532161, 0.554636, 0.563272,
                          0.566608, 0.547038, 0.563977, 0.549819, 0.555607,
                          0.545972, 0.554741, 0.555081, 0.555807, 0.59528,
                          0.599909, 0.60842, 0.608734, 0.602171, 0.574354,
                          0.601932, 0.622413, 0.6149, 0.649133, 0.648208,
                          0.646136, 0.648936, 0.628589, 0.641852, 0.623592,
                          0.649156, 0.63189, 0.617886, 0.64418, 0.638153,
                          0.632425, 0.65449, 0.661854, 0.653293, 0.641437,
                          0.685784, 0.676735, 0.651846, 0.670251, 0.675627,
                          0.686973, 0.665005, 0.664969, 0.713855, 0.688728,
                          0.66383, 0.679787, 0.693635, 0.699666, 0.686275,
                          0.735499, 0.689655, 0.736211, 0.70725, 0.740741,
                          0.71408, 0.723705, 0.741456, 0.759513, 0.768049,
                          0.772955, 0.82554, 0.828096, 0.836224, 0.857143,
                          0.879592, 0.855131, 0.867713, 0.910515, 0.892929,
                          0.924312, 0.932018, 0.914414, 0.955457, 0.944215,
                          0.961098, 0.936364, 0.973333, 0.952703, 0.980561,
                          0.968326, 0.973799, 0.985714, 0.985949, 0.99095,
                          0.986014, 0.990741, 0.987562, 0.981609, 0.974684,
                          0.989556, 0.989717, 0.987562, 0.989333, 0.995049,
                          0.982353, 0.994203, 0.99359, 0.996815, 0.991304,
                          0.985207, 0.990323, 0.98615, 0.996324, 0.983165,
                          0.993103, 0.972789, 0.982079, 0.977987, 0.984848,
                          0.834921, 0.916981, 0.939502, 0.965649, 0.952174,
                          0.965517, 0.973684, 0.97541, 0.964444, 0.982833,
                          0.986301, 0.970149, 0.990741, 0.975728, 0.995098,
                          0.983425, 0.969388, 0.994565, 0.979695, 0.97861,
                          0.98324, 0.97561, 0.988304, 0.974359, 0.980519,
                          0.987578, 1, 0.988304, 0.986014, 0.96875,
                          0.984496, 0.986301, 0.975806, 0.992187, 1,
                          0.957265, 0.97561, 0.96875, 0.98374, 0.966667,
                          0.992, 0.975207, 0.96748, 0.973684, 0.990385,
                          0.974576, 0.948718, 0.973684, 0.97561, 0.929204,
                          0.888, 0.857143, 0, 1, 1,
                          0, 0, 0, 0, 0,
                          1, 1, 0, 0, 1,
                          1, 1, 1, 0, 0,
                          0, 0, 0, 0, 1,
                          0, 0, 1, 0, 0 };
    const double fB[300] = { -480.662, -72.1688, 19.3376, -5.18149, 1.38837,
                          -0.372014, 0.0996808, -0.0267094, 0.00715676, -0.00191765,
                          0.000513832, -0.000137681, 3.68915e-05, -9.88505e-06, 2.64869e-06,
                          -7.09715e-07, 1.90167e-07, -5.09552e-08, 1.36534e-08, -3.65842e-09,
                          9.80272e-10, -2.62666e-10, 7.03934e-11, -1.89075e-11, 5.23669e-12,
                          -2.03923e-12, 2.92024e-12, -9.64172e-12, 3.56466e-11, -1.32945e-10,
                          4.96133e-10, -1.85159e-09, 6.91021e-09, -2.57893e-08, 9.62468e-08,
                          -3.59198e-07, 1.34055e-06, -5.00298e-06, 1.86714e-05, -6.96826e-05,
                          0.000260059, -0.000970553, 0.00362215, -0.0135181, 0.0504501,
                          -0.188282, 0.702679, -2.62243, 9.78706, -36.5258,
                          136.316, 111.34, -58.2392, -10.3154, -20.7924,
                          0.957295, -9.06811, -2.51904, 0.279, -2.67207,
                          3.20561, 4.14503, 0.592716, -0.836383, -1.32879,
                          1.29157, 1.28449, -1.1137, 1.03702, 6.26963,
                          1.9465, -0.288441, -4.624, -0.912341, 2.17394,
                          -1.72273, 3.89676, 1.22757, -0.159235, 0.935638,
                          0.345361, -0.00152861, -2.38944, 2.09295, 4.11948,
                          -0.800859, -3.17263, 0.0690546, -0.0614875, 1.10816,
                          1.60365, -1.64225, -0.326104, 0.66165, 2.37234,
                          0.867, -4.36607, 0.789538, 4.39778, -1.41013,
                          -2.81192, 1.70848, 0.136837, 1.89269, -2.86882,
                          0.857101, -0.046699, -2.12084, 5.95663, 1.62749,
                          -3.48723, 0.146049, 0.929469, -1.7786, -0.0924275,
                          -0.737123, 2.39155, -1.99724, 6.39644, 6.56075,
                          0.437325, 1.54519, 0.000722776, -6.23502, -0.845676,
                          9.43855, -0.864313, 3.74498, 5.92417, -2.46049,
                          1.67033, -3.67509, -0.130526, -1.11572, 0.846099,
                          3.20893, -7.45857, 3.17312, 3.98361, -3.9071,
                          2.82849, 4.84566, -0.139178, -5.18677, 5.57372,
                          7.25978, -8.13983, -0.153848, 3.89238, 2.42052,
                          -1.03305, -6.25503, 9.55066, 4.69022, -10.4926,
                          -0.23907, 4.74325, 3.62023, -4.31483, 8.11847,
                          -1.2847, -0.444237, 3.59585, -0.742789, 2.77255,
                          -5.22478, 5.3497, 4.3578, 4.07514, -0.713617,
                          8.86086, 8.38798, -1.05692, 3.85266, 7.43137,
                          -1.05193, -4.73269, 11.0736, 1.97628, -0.066475,
                          8.63762, -5.16779, 4.61043, 4.30531, 0.518682,
                          -2.14867, 2.18754, 2.57476, -0.232228, 3.77526,
                          -3.15148, 3.75888, 1.15729, 0.723992, -0.126328,
                          -0.169596, 0.647669, -1.25991, -2.45675, 1.42793,
                          2.70527, -0.973711, -0.30589, 1.90934, -1.71599,
                          -0.280706, 2.20384, -0.107027, 0.183621, -2.34152,
                          0.476281, -0.29996, 1.43041, -0.920934, 0.0148896,
                          -1.55368, -1.58204, -0.386544, 7.02695, -25.644,
                          -11.7509, 21.747, 3.1987, 1.95899, -1.53056,
                          4.06456, 1.40511, -2.26557, 0.727333, 4.9233,
                          -4.0278, 1.67546, 0.655427, -0.113014, 3.06463,
                          -6.37255, 3.14281, 2.15611, -4.03644, 2.02302,
                          -1.39703, 1.31519, -0.065847, -1.98986, 2.18683,
                          3.15653, -0.20255, -1.80148, -3.08105, -0.539922,
                          4.10241, -2.70618, 0.205052, 6.30054, -7.26205,
                          -3.44424, 2.74631, 1.07283, -0.939994, 1.12462,
                          2.63656, -5.26582, 0.0364696, 3.97812, 1.22978,
                          -8.22819, 0.432938, 5.82737, -3.57353, -24.8938,
                          37.4413, -178.917, 12.2259, 237.156, -210.85,
                          -143.757, 35.8781, 0.244646, -36.8567, 147.182,
                          198.128, -189.695, -189.349, 197.092, 150.98,
                          -51.0124, 53.0697, -161.267, -158.004, 43.2808,
                          -15.1196, 17.1976, -53.671, 197.486, 13.7261,
                          -252.39, 245.836, 19.0469, -322.023, 519.047 };
    const double fC[300] = { 70873.4, 31250, -8373.41, 2243.65, -601.184,
                          161.087, -43.1631, 11.5655, -3.09897, 0.830366,
                          -0.222496, 0.0596176, -0.0159745, 0.00428035, -0.00114692,
                          0.000307315, -8.23449e-05, 2.20643e-05, -5.9121e-06, 1.58414e-06,
                          -4.24469e-07, 1.13735e-07, -3.04698e-08, 8.14459e-09, -2.10854e-09,
                          2.89557e-10, 9.50311e-10, -4.0908e-09, 1.54129e-08, -5.75608e-08,
                          2.1483e-07, -8.0176e-07, 2.99221e-06, -1.11671e-05, 4.16761e-05,
                          -0.000155537, 0.000580473, -0.00216636, 0.00808495, -0.0301734,
                          0.112609, -0.420262, 1.56844, -5.85349, 21.8455,
                          -81.5286, 304.269, -1135.55, 4237.92, -15816.1,
                          59026.6, -65270.6, 22875.7, -10894.8, 8275.53,
                          -2838.1, 331.746, 1305.52, -606.013, -131.754,
                          1601.17, -1366.32, 478.242, -835.517, 712.415,
                          -57.3251, 55.5572, -655.105, 1192.79, 115.367,
                          -1196.15, 637.414, -1721.31, 2649.22, -1877.65,
                          903.485, 501.388, -1168.69, 821.987, -548.269,
                          400.7, -487.422, -109.555, 1230.15, -723.518,
                          -506.565, -86.3763, 896.797, -929.432, 1221.84,
                          -1097.97, 286.493, 42.5442, 204.394, 223.279,
                          -599.615, -708.652, 1997.55, -1095.49, -356.482,
                          6.03415, 1124.07, -1516.98, 1955.94, -3146.32,
                          4077.8, -4303.75, 3785.22, -1765.85, 683.564,
                          -1962.25, 2870.57, -2674.71, 1997.7, -1576.15,
                          1414.98, -632.809, -464.39, 2562.81, -2521.73,
                          990.877, -713.911, 327.795, -1886.73, 3234.06,
                          -663.009, -1912.71, 3065.03, -2520.23, 424.067,
                          608.639, -1945, 2831.14, -3077.44, 3567.89,
                          -2977.18, 310.309, 2347.61, -2144.99, 172.314,
                          1511.58, -1007.29, -238.92, -1022.98, 3713.1,
                          -3291.59, -558.318, 2554.81, -1543.26, 1175.29,
                          -2038.69, 733.19, 3218.23, -4433.34, 637.641,
                          1925.74, -680.158, 399.402, -2383.16, 5491.49,
                          -7842.28, 8052.39, -7042.37, 5957.71, -5078.88,
                          3079.54, -435.924, 187.947, -258.612, -938.577,
                          3332.2, -3450.42, 1089.19, 138.202, 756.477,
                          -2877.3, 1957.12, 1994.45, -4268.76, 3758.07,
                          -1582.05, -1869.31, 4313.86, -4390.14, 3443.48,
                          -4110.32, 5194.37, -5097.57, 4395.82, -3393.95,
                          1662.27, 65.3248, -715.723, 607.398, -819.978,
                          809.161, -604.845, 127.951, -427.163, 1398.33,
                          -1079, 159.253, 7.70185, 546.105, -1452.44,
                          1811.26, -1190.12, 612.402, -539.74, -91.5455,
                          795.996, -990.056, 1422.65, -2010.49, 2244.44,
                          -2636.58, 2629.5, -2330.62, 4183.99, -12351.7,
                          15825, -7450.56, 2813.48, -3123.41, 2251.02,
                          -852.238, 187.378, -1105.05, 1853.28, -804.285,
                          -1433.49, 2859.31, -3114.32, 2922.21, -2127.79,
                          -231.502, 2610.34, -2857.02, 1308.88, 205.98,
                          -1060.99, 1739.05, -2084.31, 1603.3, -559.133,
                          801.558, -1641.33, 1241.59, -1561.49, 2196.77,
                          -1036.19, -665.963, 1393.77, 130.1, -3520.75,
                          4475.2, -2927.56, 2509.19, -3012.4, 3528.55,
                          -3150.57, 1174.97, 150.597, 834.816, -1521.9,
                          -842.588, 3007.87, -1659.26, -690.96, -4639.1,
                          20222.9, -74312.4, 122098, -65865.5, -46135.9,
                          62909, -18000.2, 9091.85, -18367.2, 64376.9,
                          -51640.4, -45315.3, 45401.6, 51208.8, -62736.9,
                          12238.8, 13781.8, -67365.8, 68181.6, -17860.5,
                          3260.39, 4818.92, -22536.1, 85325.4, -131265,
                          64736.3, 59820.3, -116518, 31250, 0.004 };
    const double fD[300] = { -3.30195e+06, -3.30195e+06, 884755, -237069, 63522.6,
                          -17020.8, 4560.71, -1222.04, 327.445, -87.7385,
                          23.5095, -6.29934, 1.6879, -0.452272, 0.121186,
                          -0.0324717, 0.00870076, -0.00233136, 0.000624687, -0.000167384,
                          4.48503e-05, -1.20171e-05, 3.21787e-06, -8.54428e-07, 1.99841e-07,
                          5.50628e-08, -4.20093e-07, 1.62531e-06, -6.08114e-06, 2.26992e-05,
                          -8.47158e-05, 0.000316164, -0.00117994, 0.0044036, -0.0164344,
                          0.0613342, -0.228902, 0.854275, -3.1882, 11.8985,
                          -44.4059, 165.725, -618.494, 2308.25, -8614.51,
                          32149.8, -119985, 447789, -1.67117e+06, 6.23689e+06,
                          -1.03581e+07, 7.34553e+06, -2.81421e+06, 1.59753e+06, -926136,
                          264153, 81148, -159295, 39521.6, 144411,
                          -247291, 153714, -109480, 128994, -64145,
                          9406.86, -59221.9, 153991, -89784.8, -109293,
                          152797, -196560, 364211, -377239, 231761,
                          -33508.1, -139173, 165890, -114188, 79080.8,
                          -74010.2, 31488.9, 111642, -162806, 18079.4,
                          35015.8, 81931.1, -152186, 179273, -193318,
                          115372, -20329.1, 13487.5, 1573.76, -68574.5,
                          -9086.41, 225517, -257754, 61584.4, 30209.7,
                          93169.3, -220087, 289410, -425189, 602010,
                          -698463, 674081, -462589, 204118, -220484,
                          402734, -462107, 389367, -297821, 249261,
                          -170649, 14034.9, 252267, -423712, 292718,
                          -142066, 86808.8, -184544, 426733, -324756,
                          -104141, 414811, -465438, 245358, 15381,
                          -212803, 398011, -492381, 553777, -545423,
                          273958, 169775, -374384, 193109, 111606,
                          -209906, 64030.8, -65338.3, 394673, -583724,
                          227772, 259428, -341506, 226546, -267832,
                          230990, 207087, -637631, 422582, 107341,
                          -217158, 89963.3, -231881, 656221, -1.11115e+06,
                          1.32456e+06, -1.2579e+06, 1.08334e+06, -919716, 679869,
                          -292956, 51989.3, -37213.3, -56663.7, 355898,
                          -565218, 378301, -79249.2, 51522.9, -302815,
                          402868, 3110.81, -521934, 668903, -445010,
                          -23938.1, 515264, -725333, 652802, -629484,
                          775391, -857662, 791116, -649148, 421352,
                          -133079, -65087.3, 110260, -118948, 135762,
                          -117834, 61066.3, -46259.5, 152125, -206444,
                          103188, -12629.3, 44866.9, -166545, 271974,
                          -250115, 150210, -96011.8, 37349.5, 73961.8,
                          -148838, 201059, -286095, 354577, -406752,
                          438840, -413343, 542884, -1.37798e+06, 2.34806e+06,
                          -1.93963e+06, 855336, -494740, 447869, -258605,
                          86634.7, -107702, 246527, -221463, -52433.8,
                          357733, -497802, 503043, -420833, 158024,
                          236820, -455614, 347159, -91908.8, -105581,
                          233337, -318613, 307301, -180203, 113391,
                          -203574, 240244, -233590, 313188, -269413,
                          30851.9, 171644, -105306, -304237, 666329,
                          -616896, 453063, -460132, 545079, -556593,
                          360462, -85364.8, 57018.2, -196393, 56609.6,
                          320872, -388928, 80692, -329011, 2.07183e+06,
                          -7.87793e+06, 1.63675e+07, -1.56636e+07, 1.64414e+06, 9.08707e+06,
                          -6.74244e+06, 2.25767e+06, -2.28825e+06, 6.89534e+06, -9.66811e+06,
                          527091, 7.55974e+06, 483933, -9.49548e+06, 6.24797e+06,
                          128582, -6.7623e+06, 1.12956e+07, -7.17017e+06, 1.76007e+06,
                          129877, -2.27958e+06, 8.98845e+06, -1.80492e+07, 1.63335e+07,
                          -409662, -1.46948e+07, 1.2314e+07, 1.2314e+07, 0.00185641 };
    int klow=0;
    if(x<=fXmin) klow=0;
    else if(x>=fXmax) klow=fNp-1;
    else {
      if(fKstep) {
        // Equidistant knots, use histogramming
        klow = int((x-fXmin)/fDelta);
        if (klow < fNp-1) klow = fNp-1;
      } else {
        int khig=fNp-1, khalf;
        // Non equidistant knots, binary search
        while(khig-klow>1)
          if(x>fX[khalf=(klow+khig)/2]) klow=khalf;
          else khig=khalf;
      }
    }
    // Evaluate now
    double dx=x-fX[klow];
    //double ff;
    double ff=(fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
    /*
    if(x>0.135){
      if(x<0.145){
        ff=0.99;
      }
      else{
        ff=(fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
      }
    }
    else{
      ff=(fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
    }*/
    //*/

    /*
    if(x>0.135){
      if(x<0.145){
        
        return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
      }
    }
    else{
      return (fY[klow]+dx*(fB[klow]+dx*(fC[klow]+dx*fD[klow])));
    }*/
    return ff;
    
}



double func_R0(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double V_Fe=par[2];
  double V_Fe_p=V_Fe+mu_n*in_mag;
  double V_Fe_m=V_Fe-mu_n*in_mag;
  double R1,R11;
  //double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  //double qc=1.36239E-01;
  //double qc=1.30030e-01*sin(0.0231098)/sin(0.0159505);
  //double qc=1.24590e-01*sin(0.0231098/2.)/sin(0.0159505/2.);
  //double mm2=5.88427e+00 ;
  //double mm2=6.25737;
  //double qc=1.70382e-01;//
  
  

  double ww=2.5E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.;//par[3];
  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];

  //model parameter
  double qc_up=0.217;
  double mm=5.2;


  // 5.88427e+00 
  //double mm2=5.83252E+00;
  //double R1;


  double funcdown;
  if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe-mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe-mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        funcdown=1.;
      }
    }
    

    double Rup;
    if(q1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(q1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }

    /*
    double Rdown;
    if(q1<qc){
       Rdown=downrate*R0;  
    }
    
    else{
      if(q1>=qc){
        //double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        double down_R=downrate*R0/pow((1.+mm2*(q1-qc)),4);
        Rdown=down_R;//+down_R;
      }
    }
    */
    double Rdown;
    if(q1>0.13){
      if(q1<0.145){
        Rdown=1.-Rup;
      }
      else{
        //RR=Spl1(x1);
        Rdown=Spl1(q1)-Rup;
      }
    }
    else{
      //RR=Spl1(x1);
      Rdown=Spl1(q1)-Rup;
    }
    //Rdown=Spl1(q1)-Rup;
    
    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);

      R11=(R1_nume/R1_deno);//*(1.-Rdown)+funcdown*Rdown;
      //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        //R1=(R1_nume/R1_deno)*(funcP-0.5);
        R11=(R1_nume/R1_deno);//*(1.-Rdown)+funcdown*Rdown;
        //R1=(R1_nume/R1_deno)*(1.5-funcP)+funcdown*(funcP-0.5);
      }
      else{
        //R1=1.*(1.5-funcP);
        //R1=1*(funcP-0.5);
        R11=1.;//*(1.-Rdown)+funcdown*Rdown;
        //R1=1.*(1.5-funcP)+funcdown*(funcP-0.5);   
      }
    }

    R1=R11*(1.-Rdown)+funcdown*Rdown;

    return R1;

  }
  TF1 *f0 = new TF1("",func_R0,0.01,0.5,3);

double Spl2(double *x, double *par){
  double x1=x[0];
  double RR;
  //double R0=0.99;
  double ww=2.5E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;
  double mm=5.2;

  double Rup;
    if(x1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(x1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((x1-mm*qc_up)/ww))*(1.-alpha*(x1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }
  
  if(x1>0.13){
      if(x1<0.145){
        RR=1.-Rup;
      }
      else{
        RR=Spl1(x1)-Rup;
      }
  }
  else{
    RR=Spl1(x1)-Rup;
  }
  //RR=Spl1(x1);
  return RR;
}
double Spl3(double *x, double *par){
  double x1=x[0];
  double RR;
  //double R0=0.99;
  double ww=2.5E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;
  double mm=5.2;

  double Rup;
    if(x1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(x1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((x1-mm*qc_up)/ww))*(1.-alpha*(x1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }
  
  if(x1>0.13){
      if(x1<0.145){
        RR=1.-Rup;
      }
      else{
        RR=Spl1(x1)-Rup;
      }
  }
  else{
    RR=Spl1(x1)-Rup;
  }
  //RR=Spl1(x1);
  double RR1=1-RR;
  return RR1;
}

TF1 *sp2 = new TF1("",Spl2,0.1,0.5,0);
TF1 *sp3 = new TF1("",Spl3,0.1,0.5,0);
  
double func_R00(double *qqq,double *par){
  double q1=qqq[0];
  double E1=pow(hbar*q1,2)/8./m_nc2nm;
  double d1=par[0];
  double in_mag=par[1];
  double V_Fe=par[2];
  double R1;
  double V_Fe_p=V_Fe+mu_n*in_mag;
  double V_Fe_m=V_Fe-mu_n*in_mag;
  //double q1=qqq[0];
  //double E1=pow(hbar*q1,2)/8./m_nc2nm;
  //double qc=1.36239E-01;
  //double qc=1.30030e-01;
  
  //double qc=1.24590e-01*sin(0.0231098/2.)/sin(0.0159505/2.);
  //double mm2=5.88427e+00 ;
  //double mm2=6.24502;
  

  double ww=2.5E-03;//par[1];
  double R0=0.99;
  
  //double mm=par[3];
  double alpha=0.28;//par[3];

  double uprate=0.5;//par[3];
  double downrate=0.5;//par[4];
  double qc_up=0.217;

  double mm=5.2;
  //double mm2=5.83252E+00;
  
  //double R1;

  double funcdown;
  if(E1>V_Fe_m){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe-mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      
    }
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe-mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        funcdown=(R1_nume/R1_deno);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      }
      else{
        funcdown=1.;
      }
    }
    

  double Rup;
    if(q1<qc_up){
       Rup=uprate*R0;  
    }
    
    else{
      if(q1>=qc_up){
        double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        //double down_R=downrate*R0*(1.-tanh((q1-qc)/ww))*(1.-alpha*(q1-qc));
        Rup=up_R;//+down_R;
      }
    }

/*
    double Rdown;
    if(q1<qc){
       Rdown=downrate*R0;  
    }
    
    else{
      if(q1>=qc){
        //double up_R=uprate*0.5*R0*(1.-tanh((q1-mm*qc_up)/ww))*(1.-alpha*(q1-qc_up));
        double down_R=downrate*R0/pow((1.+mm2*(q1-qc)),4);
        Rdown=down_R;//+down_R;
      }
    }
    */
   double Rdown;
    if(q1>0.13){
      if(q1<0.145){
        Rdown=1.-Rup;
      }
      else{
        //RR=Spl1(x1);
        Rdown=Spl1(q1)-Rup;
      }
    }
    else{
      //RR=Spl1(x1);
      Rdown=Spl1(q1)-Rup;
    }

    if(E1>V_Fe_p){
      double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
      //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k1a=sqrt(2.*m_nc2nm*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
      double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

      double R1_nume=pow((k0*k1a-k1a*k2)*cos(k1a*d1),2)+pow((k0*k2-k1a*k1a)*sin(k1a*d1),2);
      double R1_deno=pow((k0*k1a+k1a*k2)*cos(k1a*d1),2)+pow((k1a*k1a+k0*k2)*sin(k1a*d1),2);
      //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
      //R1=(R1_nume/R1_deno)*(funcP-0.5);


      //R1=(R1_nume/R1_deno)*(funcP-0.5)+funcdown*(1.5-funcP);
      R1=(R1_nume/R1_deno)*Rdown+funcdown*(1.-Rdown);
      
    }
    
    else{
      if(E1>V_Si){
        double k0=sqrt(2.*m_nc2nm*E1)/hbar;//nm^-1
        //double k1a=sqrt(2.*m_nc2*(E1-(V_Fe+mu_n*in_mag)))/hbar;//nm^-1
        double alpha=sqrt(2.*m_nc2nm*((V_Fe+mu_n*in_mag)-E1))/hbar;//nm^-1
        double k2=sqrt(2.*m_nc2nm*(E1-V_Si))/hbar;//nm^-1

        double R1_nume=pow((k0*alpha-alpha*k2)*cosh(alpha*d1),2)+pow((k0*k2+alpha*alpha)*sinh(alpha*d1),2);
        double R1_deno=pow((k0*alpha+alpha*k2)*cosh(alpha*d1),2)+pow((alpha*alpha-k0*k2)*sinh(alpha*d1),2);
        //R1=(R1_nume/R1_deno)*(1.5-funcP);//k1a*k1a-k0*k2 k0*k2-k1a*k1a
        //R1=(R1_nume/R1_deno)*(funcP-0.5);
        //R1=(R1_nume/R1_deno)*(funcP-0.5)+funcdown*(1.5-funcP);
        R1=(R1_nume/R1_deno)*Rdown+funcdown*(1.-Rdown);
      }
      else{
        //R1=1.*(1.5-funcP);
        //R1=1*(funcP-0.5);


        //R1=1.*(funcP-0.5)+funcdown*(1.5-funcP);
        R1=1.*Rdown+funcdown*(1.-Rdown);
        
      
      }
    }
    //R1=Rdown;

    return R1;

  }
  TF1 *f1 = new TF1("",func_R00,0.01,0.5,3);



Double_t Distance = 18.101;//[m]
Double_t Conversion = 395.6;
Double_t dist_det   = 666.; //sample to detector [mm]
//Double_t xdirect    = 63.29;//62.9442
Double_t xdirect    = 60.4483;//62.9442;
Bool_t useMRfirst = 0; //use only MR events to avoid frame overlap
Bool_t useThinout = 0; //thinning out the event <1e4.

void InitColor(){
  //set default color
  gROOT->GetColor(2)->SetRGB(220./255.,  50./255.,  47./255.); // Red
  gROOT->GetColor(3)->SetRGB(135./255., 194./255.,  63./255.); // Green
  gROOT->GetColor(4)->SetRGB( 38./255., 139./255., 210./255.); // Blue
  gROOT->GetColor(5)->SetRGB(250./255., 202./255.,  18./255.); // Yellow
  gROOT->GetColor(6)->SetRGB(236./255.,   0./255., 140./255.); // Magenta
  gROOT->GetColor(7)->SetRGB(135./255., 206./255., 250./255.); // Cyan
  gROOT->GetColor(8)->SetRGB(102./255., 205./255., 170./255.); // Lightgreen
  return;
}

TTree* GetTree(TString filestr){

  TString ROOTstr = filestr(0,19);//filestr　ファイルを番号によって読むものを変える
  ROOTstr += ".root";
  TString path = "data/210713_SiFe/";
  TString ROOTstr_path = path+ROOTstr;
  TFile *file = TFile::Open(ROOTstr_path.Data());

  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  //if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}

//9/8
TTree* GetTree1(TString ROOTstr_path1){

  TFile *file = TFile::Open(ROOTstr_path1.Data());
  //TFile *file = TFile::Open(ROOTstr.Data());
  if ( file->IsOpen() ) printf("ROOT file opened successfully\n");
  TTree* tup=(TTree*)file->Get("T");
  if(useThinout==1)tup->SetMaxEntryLoop(10000);
  //  tup->SetDirectory(NULL);
  //  file->Close();
  return tup;
}
////////

const string scan_id = "90nm_scan_fine_3";
const string run_id="20210717002421";
const Int_t num1 = 5; // this should be the half of the number of the files obtained by the scan 


Int_t fit_check_Spline_global3(){

  


  
  // load CSV file with magnetic field
  vector<Int_t> vec_index;
  vector<Double_t> vec_I, vec_H; 
  // vec_index.clear(); vec_I.clear(); vec_H.clear();
  ifstream fcsv(Form("data_scans/magnetic_%s.csv", run_id.c_str())); 
  if (!fcsv.is_open())
  {
      exit(EXIT_FAILURE);
  }
  string str;
  getline(fcsv, str);
  cout << str << endl; // print out the first row
  while (getline(fcsv, str)){
      Int_t csv_i;     
      Double_t current;
      Double_t magfield;
      sscanf(str.c_str(), "%d,%lf,%lf", &csv_i, &current, &magfield);
      // str >> index >> ",'" >> current >>",'">> magfield;
      vec_index.push_back(csv_i); 
      vec_I.push_back(current);
      vec_H.push_back(magfield);
      // cout << " " << index
    //  << " " << current
      //   << " " << magfield << endl;
      // }
  // i_csv++;
  }

  InitColor();
  TH1::SetDefaultSumw2();

  const Int_t num = 4;
  Int_t kp[num];
  Int_t kp2[num];
  TTree* tup[num];
  TTree* tup2[num];
  TH1F* hx[num];
  TH1F* hlambda[num];
  TH1F* hlm0[num];
  TH1F* hlm[num];
  TH1F* hratio[num];
  TH1F* hq[num];
  TH1F* hq0[num];
  TH3F* hxylambda[num];

  TH1F* hpolratio[num];
  TH1F* hpolratio2[num];
  TH1F* hpolratio3[num];
  TH1F* hq2[num];
  TH1F* hq02[num];

  TH1F* hq0BG[num];
  TH1F* hqBG[num];
  TH1F* hq02BG[num];
  TH1F* hq2BG[num];

  Int_t kp11[num1];
  Int_t kp211[num1];
  TTree* tup11[num1];
  TTree* tup211[num1];
  //TTree* tup02[num1];

  TH1F* hx11[num1];
  TH1F* hlambda11[num1];
  TH1F* hratio11[num1];
  TH1F* hq11[num1];
  TH1F* hq011[num1];
  TH3F* hxylambda11[num1];

  TH1F* hpolratio11[num1];
  TH1F* hpolratio211[num1];
  TH1F* hpolratio311[num1];
  TH1F* hq211[num1];
  TH1F* hq0211[num1];

  TH1F* hq0BG11[num1];
  TH1F* hqBG11[num1];
  TH1F* hq02BG11[num1];
  TH1F* hq2BG11[num1];

  TTree* tup0;
  Int_t kp0;


  Double_t angle11[num1];
  Double_t angle211[num1];
  for(Int_t i=0; i< num; i++){
    angle11[i]=TMath::Abs(45.0572 - xdirect)/dist_det;
    angle211[i]=TMath::Abs(45.0572 - xdirect)/dist_det;
  }  

  Double_t angledeg11[num];
  Double_t angledeg211[num];
  
  angledeg11[0]=angle11[0]*180./TMath::Pi()/2.;
  angledeg11[1]=angle11[1]*180./TMath::Pi()/2.;
  angledeg11[2]=angle11[2]*180./TMath::Pi()/2.;

  
  
  


  TString namestr[num];
  TString namestr2[num];
  //off
  namestr[0]="20210714193654_list.root";
  namestr[1]="20210716232122_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  //namestr[2]="20210717004515_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr[2]="20210716233530_list.root"; //AFP ON
  
  namestr[3]="20210717022140_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  //namestr[4]="20210715085349_list.root"; //0.264A Fe 30 nm, theta = 0.69 deg., x = 0.0 mm, B = 1.5 mT from -8 mT
  //namestr[5]="20210715082606_list.root"; //0.378A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 0 mT  -> 2 mT
  //namestr[6]="20210715083711_list.root"; //0.6A Fe 30 nm, theta = 0.69 deg. x = 0.0 mm , B = 2 mT  -> 3 mT

  //on
  namestr2[0]="20210714193654_list.root";
  namestr2[1]="20210716233530_list.root"; //2A Fe 90 nm,  x = 0.0 mm, B = -8.13 mT
  //namestr2[2]="20210717005023_list.root"; //0A Fe 90 nm, theta = 0.69 deg. x = 0.0 mm , B = -0.32198 mT
  namestr2[2]="20210716233530_list.root"; //AFP ON
  namestr2[3]="20210717022920_list.root"; //0.265A Fe 90 nm, theta = 0.69 deg., x = 0.0 mm, B = -1.35656  mT 
  
  TString degstr[num];
  TString degstr2[num];
  //off
  degstr[0]="Direct(M1 reflect)";
  degstr[1]="B = 8.13 mT, SF OFF";
  degstr[2]="B = 8.13 mT, SF ON";
  //degstr[3]="B = 0.908 mT";
  degstr[3]="B = 1.35 mT";
  //degstr[5]="B = 1.80 mT";
  //degstr[6]="B = 2.66 mT";

  //on
  degstr2[0]="Direct(M1 reflect)";
  degstr2[1]="B = 8.13 mT";
  degstr2[2]="B = 0.322 mT";
  //degstr2[3]="B = 0.908 mT";
  degstr2[3]="B = 1.35 mT";
  //degstr2[5]="B = 1.80 mT";
  //degstr2[6]="B = 2.66 mT";
  

  Double_t angle[num];
  Double_t angle2[num];
  angle[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  //angle[1] = TMath::Abs(47.4493 - xdirect)/dist_det; //rad 47.1868(30nm)47.2685
  angle[1] = TMath::Abs(45.1074 - xdirect)/dist_det; 
  
  //angle[2] = TMath::Abs(47.09 - xdirect)/dist_det; //rad
  angle[2] = TMath::Abs(45.0572 - xdirect)/dist_det; //rad 47.04
  angle[3] = TMath::Abs(47.4013 - xdirect)/dist_det; //rad 47.2
  //angle[4] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[5] = TMath::Abs(47.2 - xdirect)/dist_det; //rad
  //angle[6] = TMath::Abs(47.2 - xdirect)/dist_det; //rad

  angle2[0] = TMath::Abs(70.5 - xdirect)/dist_det; //rad
  angle2[1] = TMath::Abs(45.1074 - xdirect)/dist_det; //rad 47.07
  angle2[2] = TMath::Abs(45.0572- xdirect)/dist_det; //rad 47.1 
  angle2[3] = TMath::Abs(47.4845 - xdirect)/dist_det; //rad 47.19
  //angle2[4] = TMath::Abs(47.24 - xdirect)/dist_det; //rad
  //angle2[5] = TMath::Abs(47.21 - xdirect)/dist_det; //rad
  //angle2[6] = TMath::Abs(47.11 - xdirect)/dist_det; //rad

  Double_t angledeg[num1];
  Double_t angledeg2[num1];
  
  angledeg[0]=angle[0]*180./TMath::Pi()/2.;
  angledeg[1]=angle[1]*180./TMath::Pi()/2.;
  angledeg[2]=angle[2]*180./TMath::Pi()/2.;
  angledeg[3]=angle[3]*180./TMath::Pi()/2.;

  for(int i=0; i<3; i++){
    cout<<angle[i]<<"_[rad]_"<<angledeg[i]<<"_[deg]_"<<endl;
  }

  //  TLegend* leg = new TLegend(0.15, 0.75, 0.4, 0.98,"");
  //TLegend* leg = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm OFF");
  //TLegend* leg2 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm ON");
  //TLegend* leg3 = new TLegend(0.70, 0.20, 0.98, 0.70,"Fe 90 nm (ON-OFF)/(OFF+ON)");
  TLegend* leg = new TLegend(0.6, 0.60, 1.0, 0.80,"");
  TLegend* leg2 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm ON");
  TLegend* leg3 = new TLegend(0.8, 0.20, 1.0, 0.70,"Fe 90 nm");

  TLegend* leg4 = new TLegend(0.6, 0.40, 1.0, 0.60,"");


  leg->SetFillColor(0);

  Int_t nbin = 512;
  Double_t range = 128.;
  Int_t nbin_lambda = 200;
  Double_t lambda_max  = 1.5;
  Int_t nbin_q  = 60;//300 60
  Double_t q_min  = 0.1;//0.6
  Double_t q_max  = 0.50;//0.6
  //Double_t q_max  = 1.0;//0.6
  Int_t nrebinx = 1;
  Int_t nrebiny = 2;

  Int_t LLD  = 500.;
  //  Double_t HLD  = 7400.;

  //  Double_t xbegin=54.;
  /*
  Double_t xbegin=19.;
  Double_t xcenter=29.;
  Double_t xend=40.;
  */
  /*
  Double_t xbegin1=19.;
  Double_t xcenter1=29.;
  Double_t xend1=40.;
  */
  
  //Double_t xbegin1=20.;
  //Double_t xcenter1=35.;
  Double_t xbegin1=40.;
  Double_t xcenter1=55.;

  Double_t xcenter=55.;
  Double_t xend=71.;
  //Double_t xcenter=55.;
  //Double_t xend=71.;
  
  //Double_t ybegin1=43.;
  //Double_t yend1=60.;
  Double_t ybegin1=85.;
  Double_t yend1=102.;
  //Double_t ybegin1=55.;
  //Double_t yend1=71.;

  
  Double_t ybegin=55.;
  Double_t yend=71.;

  //  Double_t ybegin=70.;
  //  Double_t yend=77.;

  TString degstr11[num1];
  TString degstr211[num1];
  TString degstr_p11[num1];

  for(Int_t i=0; i< num1; i++){
    degstr11[i]=Form("SF:ON,  B=%lf mT", vec_H[i]);
    degstr211[i]=Form("SF:OFF, B=%lf mT", vec_H[i]);
    degstr_p11[i]=Form("B=%lf mT", vec_H[i]);

  } 

  TCut cut_rpmt_basic = Form("a>%d && b>%d && c>%d && d>%d && a<%d && b<%d && c< %d && d<%d && f==4",
			     LLD,LLD,LLD,LLD,rpmt_HLD,rpmt_HLD,rpmt_HLD,rpmt_HLD);
  TCut cut_x = Form("x*%f>20 && x*%f<100",range,range);
  TCut cut_y = Form("y*%f>%f && y*%f<%f",range,ybegin,range,yend);
  TCut cut_y1 = Form("y*%f>%f && y*%f<%f",range,ybegin1,range,yend1);
  
  TCut cut_dir = Form("x*%f>%f && x*%f<%f",range,xcenter,range,xend);
  TCut cut_ref = Form("x*%f>%f && x*%f<%f",range,xbegin1,range,xcenter);
  TCut cut_ref1 = Form("x*%f>%f && x*%f<%f",range,xbegin1,range,xcenter1);
  //  TCut cut_tof = Form("tof>1.0e3 && tof<39.9e3");
  TCut cut_tof = "";
  TCut MRcut = "MRflag>0";
  //TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  TCut thecut1 = cut_rpmt_basic && cut_x && cut_y1 && cut_tof;
  TCut thecut = cut_rpmt_basic && cut_x && cut_y && cut_tof;
  if(useMRfirst) thecut = thecut && MRcut;
  TCut thecut0;
  TCut thecut01;

  TCanvas *c1 = new TCanvas("c1","",1200,800);
  c1->Divide(2,2);
  c1->cd(1);

  TString namestr11[num1];
  TString namestr211[num1];
  
  TString namestr_ref= "data/210713_SiFe/20210714193654_list.root"; // file path of the direct data
  //0N
  for(Int_t i=0; i<num1; i++){
    int iscan=i*2;
    namestr11[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  //OFF
  for(Int_t i=0; i<num1; i++){
    int iscan=i*2+1;
    namestr211[i]=Form("data_scans/%s_list_%02d.root",run_id.c_str(), iscan);
  }
  
  tup0 = GetTree1(namestr_ref); // direct data 
    tup0->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); 
    if(useMRfirst) kp0 = tup0->GetMaximum("mp");
      else kp0= (tup0->GetMaximum("kp") - tup0->GetMinimum("kp"));
  cout << "direct data # of kp: "<< kp0 <<endl;


  // for(Int_t i=0; i<2; i++){

  for(Int_t i=0; i<num; i++){
    thecut.Print();
    if(i==0) thecut0=thecut;
    if(i==0) thecut01=thecut1;

    Double_t twopirad = 2*TMath::Pi()*angle[i];
    Double_t twopirad2 = 2*TMath::Pi()*angle2[i];
    Double_t lambda_coeff = 1.e-6*Conversion/Distance;

    Double_t twopirad11 = 2*TMath::Pi()*angle11[i];
    Double_t twopirad211 = 2*TMath::Pi()*angle211[i];

    // Direct data

    hq011[i] = new TH1F(Form("hq011%d",i),Form("%s;q [nm^{-1}];count/bin/25kp","Direct"),nbin_q,q_min,q_max);
    tup0->Draw(Form("%f/(toffo*%f)>>hq011%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    hq011[i]->Scale(25./kp0);

     // Data with SF:ON
    tup11[i] = GetTree1(namestr11[i]);
    tup11[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    if(useMRfirst) kp11[i] = tup11[i]->GetMaximum("mp");
      else kp11[i] = (tup11[i]->GetMaximum("kp") - tup11[i]->GetMinimum("kp"));
    cout << "SF:ON data # of kp: "<< kp11[i] <<endl;

    hq11[i] = new TH1F(Form("hq11%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr11[i].Data()),nbin_q,q_min,q_max);
    // tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup11[i]->Draw(Form("%f/(toffo*%f)>>hq11%d",twopirad11,lambda_coeff,i), thecut && cut_ref,"goff");
    hq11[i]->Scale(25./kp11[i]);
    hq11[i]->Divide(hq011[i]); //Divide by the direct data
    hq11[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq11[i]->SetTitle("Reflectivity (SF ON)");
    //leg->AddEntry(hq1[i],degstr[i],"l");


    // Direct data

    hq0211[i] = new TH1F(Form("hq0211%d",i),Form("%s;q [nm^{-1}];count/bin/25kp","Direct"),nbin_q,q_min,q_max);
    tup0->Draw(Form("%f/(toffo*%f)>>hq0211%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    hq0211[i]->Scale(25./kp0);

     // Data with SF:ON
    tup211[i] = GetTree1(namestr211[i]);
    tup211[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    if(useMRfirst) kp211[i] = tup211[i]->GetMaximum("mp");
      else kp211[i] = (tup211[i]->GetMaximum("kp") - tup211[i]->GetMinimum("kp"));
    cout << "SF:ON data # of kp: "<< kp211[i] <<endl;

    hq211[i] = new TH1F(Form("hq211%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr211[i].Data()),nbin_q,q_min,q_max);
    // tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup211[i]->Draw(Form("%f/(toffo*%f)>>hq211%d",twopirad211,lambda_coeff,i), thecut && cut_ref,"goff");
    hq211[i]->Scale(25./kp211[i]);
    hq211[i]->Divide(hq0211[i]); //Divide by the direct data
    hq211[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq211[i]->SetTitle("Reflectivity (SF ON)");
    //leg->AddEntry(hq1[i],degstr[i],"l");



    tup[i] = GetTree(namestr[i]);
    // tup[i]->SetAlias("toffo","(tof>9.e3)*(tof)+(tof<9.e3)*(tof+40.e3)");
    tup[i]->SetAlias("toffo","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // editted based on suggestion by KM on the August 3rd
    
    tup2[i] = GetTree(namestr2[i]);
    tup2[i]->SetAlias("toffo2","(tof>10.e3)*(tof)+(tof<10.e3)*(tof+40.e3)"); // edited based on suggestion by KM on the August 3rd
    

    //tup[i]->SetAlias("toffo","tof");
    if(useMRfirst) kp[i] = tup[i]->GetMaximum("mp");
    else kp[i] = tup[i]->GetMaximum("kp");
    cout << kp[i]<<endl;

    //hx[i] = new TH1F(Form("hx%d",i),Form("%s;X [mm];count/bin/25kp",degstr[i].Data()),nbin,0.,range);
    hlambda[i] = new TH1F(Form("hlambda%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    
    hlm0[i] = new TH1F(Form("hlm0%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    hlm[i] = new TH1F(Form("hlm%d",i),Form("%s;Wavelength [nm];count/bin/25k",degstr[i].Data()),nbin_lambda,0.,lambda_max);
    
    
    
    hq0[i] = new TH1F(Form("hq0%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq[i] = new TH1F(Form("hq%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq0BG[i] = new TH1F(Form("hq0BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hqBG[i] = new TH1F(Form("hqBG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    
    // hxylambda[i] = new TH3F(Form("hxylambda%d",i),Form("%s;x [mm]; y [mm]; Wavelength [nm];count/bin/25kp",degstr[i].Data()), nbin/nrebinx,0.,range,nbin/nrebiny,0.,range,nbin_lambda,0.,lambda_max);
  
    //tup[i]->Draw(Form("x*%f>>hx%d",range,i), thecut,"goff");
    if(i==0) tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_dir,"goff");
    else tup[i]->Draw(Form("toffo*%f>>hlambda%d",lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0%d",twopirad,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hq%d",twopirad,lambda_coeff,i), thecut && cut_ref,"goff");
    tup[0]->Draw(Form("%f/(toffo*%f)>>hq0BG%d",twopirad,lambda_coeff,i), thecut01 && cut_dir,"goff");
    tup[i]->Draw(Form("%f/(toffo*%f)>>hqBG%d",twopirad,lambda_coeff,i), thecut1 && cut_ref1,"goff");
    
    tup[0]->Draw(Form("toffo*%f>>hlm0%d",lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup[i]->Draw(Form("toffo*%f>>hlm%d",lambda_coeff,i), thecut && cut_ref,"goff");
    


    // tup[i]->Draw(Form("toffo*%f:y*%f:x*%f>>hxylambda%d",lambda_coeff,range,range,i),cut_rpmt_basic && MRcut,"goff");

    
    //hpolratio[i]->Rebin(10);
    //hpolratio2[i]->Rebin(10);

    //hpolratio[i]->Divide(hpolratio2[i]);


    if(i!=0){

      if(i!=3){
        leg->AddEntry(hq[i],degstr[i],"l");
      }
      
    }

   
    if(i==1) leg4->AddEntry(hratio[i],"upspin","l");
    if(i==2) leg4->AddEntry(hratio[i],"downspin","l");


    //leg->AddEntry(hq[i],degstr[i],"l");}
    leg2->AddEntry(hq[i],degstr[i],"l");
    if(i!=0)leg3->AddEntry(hq[i],degstr[i],"l");
    //leg->AddEntry(hq2[i],degstr2[i],"l2");

    //hx[i]->Scale(25./kp[i]);

    //9/6add
    hq[i]->Add(hq[i], hqBG[i],1., -1.);
    hq0[i]->Add(hq0[i], hq0BG[i],1., -1.);

    hlambda[i]->Scale(25./kp[i]);
    hlm0[i]->Scale(25./kp[i]);
    hlm[i]->Scale(25./kp[i]);
    hq[i]->Scale(25./kp[i]);
    hq0[i]->Scale(25./kp[0]);
    hqBG[i]->Scale(25./kp[i]);
    hq0BG[i]->Scale(25./kp[0]);

    for(int i11=0; i11<nbin; i11++){
      double ff[nbin];
      double ff1[nbin];
      if(i==1){
        ff[i11]=hlm[i]->GetBinContent(i11);
        ff1[i11]=hlm0[i]->GetBinContent(i11);
      }
      //cout<<"hlm_"<<ff[i11]<<"_hlm0_"<<ff1[i11]<<endl;
    }

    
    


    
    //hxylambda[i]->Scale(25./kp[i]);
    
    
    
    hratio[i]=(TH1F*)hlambda[i]->Clone(Form("hratio%d",i));
    hratio[i]->Divide(hlambda[0]);
    hratio[i]->GetYaxis()->SetTitle("Reflectivity");
    
    hlm[i]->Divide(hlm0[i]);
    hlm[i]->GetYaxis()->SetTitle("Reflectivity");
    hq[i]->Divide(hq0[i]);
    hq[i]->GetYaxis()->SetTitle("Reflectivity");
    //hq[i]->SetTitle("Reflectivity (SF OFF)");
    //hq2[i]->Divide(hq02[i]);

    if(useMRfirst) kp2[i] = tup2[i]->GetMaximum("mp");
    else kp2[i] = tup2[i]->GetMaximum("kp");
    //cout << kp2[i]<<endl;
    hq02[i] = new TH1F(Form("hq02%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    hq2[i] = new TH1F(Form("hq2%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr2[i].Data()),nbin_q,q_min,q_max);
    hq02BG[i] = new TH1F(Form("hq02BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    hq2BG[i] = new TH1F(Form("hq2BG%d",i),Form("%s;q [nm^{-1}];count/bin/25kp",degstr[i].Data()),nbin_q,q_min,q_max);
    
    tup2[0]->Draw(Form("%f/(toffo2*%f)>>hq02%d",twopirad2,lambda_coeff,i), thecut0 && cut_dir,"goff");
    tup2[i]->Draw(Form("%f/(toffo2*%f)>>hq2%d",twopirad2,lambda_coeff,i), thecut && cut_ref,"goff");
    tup2[0]->Draw(Form("%f/(toffo*%f)>>hq02BG%d",twopirad2,lambda_coeff,i), thecut01 && cut_dir,"goff");
    tup2[i]->Draw(Form("%f/(toffo*%f)>>hq2BG%d",twopirad2,lambda_coeff,i), thecut1 && cut_ref1,"goff");

    hq2[i]->Add(hq2[i], hq2BG[i],1., -1.);
    hq02[i]->Add(hq02[i], hq02BG[i],1., -1.);
    hq2[i]->Scale(25./kp2[i]);
    hq02[i]->Scale(25./kp2[0]);
    hq2BG[i]->Scale(25./kp2[i]);
    hq02BG[i]->Scale(25./kp2[0]);

    hq2[i]->Divide(hq02[i]);

    hq2[i]->GetYaxis()->SetTitle("Reflectivity");
    hq[i]->SetTitle("");
    hq2[i]->SetTitle("");
    //hq02[i]->GetYaxis()->SetTitle("Reflectivity02");

    hpolratio[i]=(TH1F*)hq[i]->Clone(Form("hpolratio%d",i));
    hpolratio[i]->Add(hq[i], hq2[i],1., -1.); // hq: OFF, hq2: ON, calculate Non - Noff
    // hpolratio[i]->Add(hq2[i],-1.);//各binの値=x*hist1の値+y*hist2の値
    hpolratio2[i]=(TH1F*)hq[i]->Clone(Form("hpolratio2%d",i));
    hpolratio2[i]->Add(hq[i], hq2[i],1., 1.); // hq: OFF, hq2: ON, calculate Non + Noff
    hpolratio[i]->Divide(hpolratio2[i]);
    //hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{on}-R_{off})/(R_{on}+R_{off})");
    hpolratio[i]->GetYaxis()->SetTitle("Polarization power (R_{off}-R_{on})/(R_{off}+R_{on})");    
    hpolratio[i]->SetTitle("Polarization power");

    
    

    if(i==9){
      //hx[i]->SetLineColor(i+2);
      hlambda[i]->SetLineColor(i+2);
      hlm0[i]->SetLineColor(i+2);
      hlm[i]->SetLineColor(i+2);
      hratio[i]->SetLineColor(i+2);
      //hq[i]->SetLineColor(i+2);
      hq0[i]->SetLineColor(i+2);
      hq2[i]->SetLineColor(i+2);
      hqBG[i]->SetLineColor(i+2);
      hq0BG[i]->SetLineColor(i+2);
      hq02BG[i]->SetLineColor(i+2);
      hq2BG[i]->SetLineColor(i+2);

      hpolratio[i]->SetLineColor(i+2);
      hpolratio2[i]->SetLineColor(i+2);
    } 
    else {
      //hx[i]->SetLineColor(i+1);
      hlambda[i]->SetLineColor(i+1);
      hlm0[i]->SetLineColor(i+1);
      hlm[i]->SetLineColor(i+1);
      hratio[i]->SetLineColor(i+1);
      hq[i]->SetLineColor(i+1);
      hq0[i]->SetLineColor(i+1);
      hq2[i]->SetLineColor(i+1);
      hpolratio[i]->SetLineColor(i+1);
      hpolratio2[i]->SetLineColor(i+1);
      hqBG[i]->SetLineColor(i+1);
      hq0BG[i]->SetLineColor(i+1);
      hq02BG[i]->SetLineColor(i+1);
      hq2BG[i]->SetLineColor(i+1);
    }

    //hpolratio[i]->Rebin(10);
    //hpolratio[i]->Scale(10);
    //hpolratio3[i]->hpolratio[i]/10.;
/*
    c1->cd(1);
    if(i==1)hpolratio2[i]->Draw("eh");
    else hpolratio2[i]->Draw("ehsames");
    leg->Draw();
    
 */ 
    TH1F* hq11;
    if(i==1)hq11=(TH1F*)hq[i]->Clone(Form("hq%d",i));
    TH1F* hq22;
    if(i==2)hq22=(TH1F*)hq[i]->Clone(Form("hq%d",i));

    ROOT::Math::WrappedMultiTF1 wfB(*f0,1);
    ROOT::Math::WrappedMultiTF1 wfSB(*f1,1);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange rangeB;
    // set the data range
    //rangeB.SetRange(10,90);
    rangeB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataB(opt,rangeB);
    if(i==2)ROOT::Fit::FillData(dataB, hq[i-1]);

    ROOT::Fit::DataRange rangeSB;
    //rangeSB.SetRange(10,50);
    rangeSB.SetRange(0.1,0.47);
    ROOT::Fit::BinData dataSB(opt,rangeSB);
    if(i==2)ROOT::Fit::FillData(dataSB, hq[i]);

    ROOT::Fit::Chi2Function chi2_B(dataB, wfB);
    ROOT::Fit::Chi2Function chi2_SB(dataSB, wfSB);

    GlobalChi2 globalChi2(chi2_B, chi2_SB);

    ROOT::Fit::Fitter fitter;

    const int Npar = 3;
    double par0[Npar] = {89.066,2., 209.0602};
  
    // create before the parameter settings in order to fix or set range on them
    fitter.Config().SetParamsSettings(3,par0);
    // fix 5-th parameter
    //fitter.Config().ParSettings(2).Fix();
    // set limits on the third and 4-th parameter
    //fitter.Config().ParSettings(2).SetLimits(-10,-1.E-4);
    //fitter.Config().ParSettings(3).SetLimits(0,10000);

    fitter.Config().ParSettings(0).Fix();
    //fitter.Config().ParSettings(0).SetLimits(88.,100.);
    fitter.Config().ParSettings(1).SetLimits(1.8,2.2);
    fitter.Config().ParSettings(2).SetLimits(150.,220.);
   
    
    //fitter.Config().ParSettings(3).SetStepSize(5);
  
    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2","Migrad");
  
    // fit FCN function directly
    // (specify optionally data size and flag to indicate that is a chi2 fit)
    fitter.FitFCN(3,globalChi2,0,dataB.Size()+dataSB.Size(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);

    c1->cd(1);
    gStyle->SetOptFit(1111);

    f0->SetNpx(10000);
    f0->SetFitResult( result, iparB);
    f0->SetRange(rangeB().first, rangeB().second);
    f0->SetLineColor(kRed);
    if(i==2)hq[i-1]->GetListOfFunctions()->Add(f0);
    //if(i==1)f0->Draw("same");
    if(i==2)hq[i-1]->Draw("eh");

    f1->SetNpx(10000);
    f1->SetFitResult( result, iparSB);
    f1->SetRange(rangeSB().first, rangeSB().second);
    f1->SetLineColor(kBlue);
    if(i==2)hq[i]->SetLineColor(4);
    if(i==2)hq[i]->GetListOfFunctions()->Add(f1);
    if(i==2)hq[i]->Draw("ehsame");

    
    
    
    leg->Draw();
    //if(i!=0){
    
    
      if(i==1){
        hq[i]->SetLineColor(2);
        hq[i]->Draw("eh");

        //hq11[i]->Draw("ehsame");
      }
      if(i==2){
        hq[i]->SetLineColor(4);
        hq[i]->Draw("ehsame");
        //hq11[i]->Draw("ehsame");
      }
      
      
      //else hq[i]->Draw("ehsames");
      
      //f0->SetParLimits(0,80.e-9,100.e-9);
      
      //f0->SetParameter(0.,30.);//nm 
      //f0->SetParameter(1,2.);
      //f0->SetParLimits(0,.e-9,100.e-9);
      
      //f0->SetParLimits(1,1.8,2.3);
      //f1->SetParLimits(1,1.5,2.2);
      f0->FixParameter(1.,2.);//nm 
      f1->FixParameter(1.,2.);//nm 
      

      //f0->SetParLimits(2,190,211);
      f0->FixParameter(0.,94.37);//nm 


      //f0->FixParameter(1.,2.);//nm 
      f0->FixParameter(2.,209.0602);//nm 

      //f0->SetParLimits(2,190,211);
      //hq[i]->Fit("f0","R","10000",0.173,0.5);
      f0->SetNpx(10000);
      //f0->SetParLimits(1,1.8,1.999);
      //f0->Draw("sames");

      //f1->SetParameter(0.,90.e-9);//m //第１引数が変数の番号、第２引数がその値
      //f1->SetParameter(0.,30.);//nm 
      //f1->SetParameter(1,2.);

      f1->SetParLimits(0,80.,100.);
      //f1->FixParameter(0.,94.37e-9);
      
      //f1->SetParLimits(2,150,220);
      //f1->FixParameter(1.,2.);//nm 
      //
      f1->FixParameter(2.,209.0602);//nm 

      //f0->SetParameter(0.,30.);//nm 
      //f0->SetParameter(1,2.);
      
      //f1->SetParLimits(1,1.8,1.999);
      f1->SetNpx(10000);
      f1->SetLineColor(4);
      //f1->Draw("sames");

      //f2->SetParameter(0.,90.e-9);//m //第１引数が変数の番号、第２引数がその値
      //f1->SetParameter(0.,30.);//nm 
      //f2->SetParameter(1,2.);
      
      //f2->SetNpx(10000);
      //f2->Draw("sames");

      double qq;
      double qq1;
      double EE1;
      qq=sqrt(8*m_nc2nm*(209.062+60.3*1.98))/hbar;
      qq1=sqrt(8*m_nc2nm*(60.3*0.007))/hbar;
      EE1=pow(hbar*0.25,2)/(8*m_nc2nm);
      cout<<"qqqqq"<<qq<<endl;
      cout<<"qqqqq1"<<qq1<<endl;
      cout<<"EE1"<<EE1<<endl;


      if(i==1){
        //hq[1]->Fit("f0","","",0.252,0.5);
        //hq[i]->Fit(f0,"+","",0.1,0.47);
      }
      if(i==2){
        //hq[2]->Fit("f1","","",0.252,0.5);
        //hq[i]->Fit(f1,"+","",0.1,0.47);
      }
      //gStyle->SetOptFit(0101);
      
    

    //}
    
    
    c1->cd(2);
    //gStyle->SetOptFit(1111);
    gPad->SetLogy();
    leg->Draw();
    if(i!=0){
    if(i==1){
      hq[i]->Draw("eh");
      //hq11[i]->Draw("ehsame");
    }
    if(i==2){
      hq[i]->Draw("ehsame");
      //hq11[i]->Draw("ehsame");
    }
    
    if(i==1){
      //hq[1]->Fit("f0","","",0.252,0.5);
      //hq[i]->Fit(f0,"","",0.1,0.47);
    }
    if(i==2){
      //hq[2]->Fit("f1","","",0.252,0.5);
      //hq[i]->Fit(f1,"+","",0.1,0.47);
    }
    //gStyle->SetOptFit(0101);
    
    hq[i]->SetStats(0); //非表示
  

  }
    

    

   
    


    double ymax=1.1;//gPad->GetUymax();
    double ymin=0.;//gPad->GetUymin();
    double xmax1=0.2;
    double xmin1=0.4;


    c1->cd(3);
    hratio[i]->SetStats(0); //非表示
    if(i!=0){
      if(i==1){
        //hlm[i]->Draw("eh");
        
        hratio[i]->Draw("eh");
        hratio[i]->SetLineColor(4);
        //hq11[i]->Draw("ehsame");
      }
      if(i==2){
        hratio[i]->SetLineColor(2);
        hratio[i]->Draw("ehsame");
        //hlm[i]->Draw("ehsame");
        //hq11[i]->Draw("ehsame");
      }
    }
  
    leg->Draw();

    
    c1->cd(4);
    
    leg4->Draw();
    c1->DrawFrame(0.1,0.,0.5,1.1);
    sp2->SetNpx(10000);
    sp2->SetLineColor(4);
    sp2->Draw("");
    sp3->SetNpx(10000);
    //sp3->SetLineColor(2);
    sp3->Draw("same");

    //}
    
    hratio[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hratio[i]->GetYaxis()->SetRangeUser(0.,1.);
    hpolratio2[i]->GetYaxis()->SetRangeUser(-1.,1.);
    hq[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    hq2[i]->GetXaxis()->SetRangeUser(q_min,q_max);
    hq2[i]->GetYaxis()->SetRangeUser(1.e-3,1.1);
    // hq[i]->SaveAs(path_R + Form("hq_off_%d.root", i));
    // hq2[i]->SaveAs(path_R + Form("hq_on_%d.root", i));
    //cout<<"in_"<<angledeg[i]<<"_deg"<<endl;
  }

  
#if 1
  TFile *outfile = TFile::Open(path_R+"pol_check_ichi2_90nm.root","RECREATE");
  //TFile *outfile = TFile::Open(path_R+"pol_check_ichi2_90nm.csv","RECREATE");
  for(Int_t i=0; i<num; i++){
    //hx[i]->Write();
    //hlambda[i]->Write();
    //hratio[i]->Write();
    if(i==1)hq[i]->Write();
    if(i==2)hq[i]->Write();
    //hxylambda[i]->Write();
    //hq2[i]->Write();
  }
  outfile->Close();
#endif
for(int i=0; i<3;i++){
    ofstream ofs(path_R+Form("hq_Rup_Rdown_%s.csv", scan_id.c_str()));  // ファイルパスを指定する
    
    for(int i11=0; i11<nbin_q; i11++){
      double xc[nbin_q];
      double xc1[nbin_q];
      double ff[nbin_q];
      double ffE[nbin_q];
      double ff1[nbin_q];
      double ff1E[nbin_q];
      if(i==1){
        ff[i11]=hq[i]->GetBinContent(i11);
        ffE[i11]=hq[i]->GetBinError(i11);
        xc[i11]=hq[i]->GetXaxis()->GetBinCenter(i11); // 
      }
      if(i==2){
        ff1[i11]=hq[i]->GetBinContent(i11);
        ff1E[i11]=hq[i]->GetBinError(i11);
        xc1[i11]=hq[i]->GetXaxis()->GetBinCenter(i11); // 
      }
      //cout<<"hlm_"<<ff[i11]<<"_hlm0_"<<ff1[i11]<<endl;
      //ofs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
      
      ofs << xc[i11] << ","<<ff[i11] << ","<<ffE[i11] << ","<< ff1[i11] << ","<< ff1E[i11] << endl;
    }
}


  return 0 ;
}

