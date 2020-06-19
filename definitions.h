#define CYL 1 // did we use cylindrification in our simulation?

// General controls: 0 or 1 to tell the code to do or not do something
// LEAVE CHARACTERIZE AT 0 FOR NOW -- NOT WORKING YET
#define CHARACTERIZE 1
#define SHEETS 0

// Threshold values for our criteria to find current sheets
#define SIGMAPHI_THR 10e-5
#define BETAPL_THR 5.0
//#define BU3_THR 10e-3 #0.00775
#define JPEAK_FAC 0.5
#define J_FAC 5.

// Box sizes
#define HALF_BOX_SIZE 5 // to look for local maxima in J_cs
#define SHEET_HALF_BOX_SIZE 5 // to find slope

//#define N1 (512)
//#define N2 (256)
//#define N3 (72)
#define GR (1)
#define NDIM (4)
#define SMALL (1.e-20)
#define DELTA (1.e-5)
#define DELTAR (1.e-4)

// Cylindrification in HARMPI
#define DOCYL 1
#define SINGSMALL 1.e-20
#define COORDSINGFIX 1

#define BETWEEN(value, min, max) (value < max && value >= min)

// all constants in cgs units
#define ME 9.1093826e-28   // electron mass (g)
#define MP 1.67262171e-24  // proton mass (g)
#define CL 2.99792458e10   // speed of light (m/s)
#define GNEWT 6.6742e-8    // gravitational constant
#define KBOL 1.3806505e-16 // Boltzmann constant (erg/K)
#define SIGMA_THOMSON 6.65245873e-25 // Thomson cross-section (cm^2)
#define MSUN 1.989e33      // solar mass (g)
#define LSUN 3.827e33      // solar luminosity (erg/s)
#define YEAR 31536000      // seconds in a year

// temperature and beta-prescription (Mościbrodzka 2016)
#define TPTE_DISK 20. // R_high
#define TPTE_JET 1.   // R_low
#define THETAE_MAX 1000.
#define TP_OVER_TE 100.

#define MBH 4.5e6 * MSUN // Sgr A*
//#define MBH 6.2e9 * MSUN // M87
//#define MBH 5.0e9 * MSUN
//#define MBH 10.0 * MSUN
#define MUNIT (1.e19)

/* some useful macros */
#define DLOOP  for(k=0;k<NDIM;k++)for(l=0;l<NDIM;l++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)

#define FAIL_METRIC (6)

// Constants, definitions and units
extern double M_unit, L_unit, T_unit, RHO_unit, U_unit, B_unit, Ne_unit;

extern char dump_file[100], jcs_output[100], jchar_output[100], jpeak_output[100];
extern char path[100];
extern char gdump[100];

extern double *****grid_gcov, *****grid_gcon, ***grid_gdet, ******grid_conn;

extern double x1in;

extern double t;
extern double gam;
extern double a;
extern double Rin;
extern double Rout;
extern double hslope;
extern double R0;
extern double x1br, cpow2, npow2, rbr;

// cylindrification variables
extern int DOCYLINDRIFYCOORDS, BL;
extern double global_x10, global_x20;
extern double fracdisk, fracjet, r0disk, rdiskend, r0jet, rjetend, jetnu, rsjet, r0grid;

extern double startx[NDIM];
extern double stopx[NDIM];
extern double dx[NDIM];
extern double dt;

extern double ***dJdx, ***dJdy, ***dJdz;

extern double ***J_cs;
extern double ***J_cs_peak;
extern double ***J_cs_char;
extern int ***flag;
extern int ***flag_buffer;

extern int *icoords;
extern int *jcoords;
extern int *kcoords;

extern int N1;
extern int N2;
extern int N3;
extern double ***a_r;
extern double ***a_th;
extern double ***a_phi;
extern double ***rho;
extern double ug;
extern double ****V;
extern double ****B;
extern double divb;
extern double gdet;
extern double ****ucov;
extern double ****ucon;
extern double ****bcov;
extern double ****bcon;
extern double ****jcov;
extern double ****jcon;
extern double jsq;
extern double jdotu;
extern double Jsq;
extern double gJsq;
extern double ***J;        // current
extern double ***sigmaphi; // phi-component of magnetization
extern double ***Sigmaphi;
extern double pg;        // gas pressure
//extern double K;
//extern double EF;       // enthalpy
extern double bsq;
//extern double EE;       // eta
//extern double va2;      // relativistic alfvén speed squared
//extern double va;       // relativistic alfvén speed
//extern double cs2;      // speed of sound squared
//extern double cms2;     // magnetosonic speed squared
//extern double T;
//extern double ptot;
extern double ***betapl;   // plasma beta
//extern double ibetapl;  // inverse plasma beta
//extern double sigma;    // magnetization
//extern double eflem;    // energy flux, EM
//extern double eflma;    // energy flux, matter
//extern double lflem;    // angular momentum flux, EM
//extern double lflma;    // angular momentum flux, matter
