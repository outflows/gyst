// Constants, definitions and units
double M_unit, L_unit, T_unit, RHO_unit, U_unit, B_unit, Ne_unit;

char dump_file[100], jcs_output[100], normal_output[100], jpeak_output[100];
char gdump_file[100];

double gcov[NDIM][NDIM], gcon[NDIM][NDIM], gdet, conn[NDIM][NDIM][NDIM];

double t;
double gam;
double a;
double Rin;
double Rout;
double hslope;
double R0;
double x1br, cpow2, npow2, rbr;

// cylindrification variables
int DOCYLINDRIFYCOORDS, BL;
double global_x10, global_x20;
double fracdisk, fracjet, r0disk, rdiskend, r0jet, rjetend, jetnu, rsjet, r0grid;

double startx[NDIM];
double dx[NDIM];
double dt;

double ***J_cs;
int ***flag;
int ***flag_buffer;

int *icoords;
int *jcoords;
int *kcoords;

int N1;
int N2;
int N3;
double ***r;
double ***th;
double ***phi;
double rho;
double ug;
double v1;
double v2;
double v3;
double B1;
double B2;
double ***B3;
double divb;
double gdet;
double ju0;
double jd0;
double ju1;
double jd1;
double ju2;
double jd2;
double ju3;
double jd3;
//double uu0;
double ud0;
//double uu1;
double ud1;
//double uu2;
double ud2;
//double uu3;
double ud3;
double bu0;
double bd0;
double bu1;
double bd1;
double bu2;
double bd2;
double ***bu3;
double bd3;

double jsq;
double jdotu;
double Jsq;
double gJsq;
double ***J;        // current
double ***sigmaphi; // phi-component of magnetization
double ***Sigmaphi; // phi-component of magnetization
double pg;        // gas pressure
//double K;
//double EF;       // enthalpy
double bsq;
//double EE;       // eta
//double va2;      // relativistic alfvén speed squared
//double va;      // relativistic alfvén speed
//double cs2;      // speed of sound squared
//double cms2;     // magnetosonic speed squared
//double T;
//double ptot;
double ***betapl;   // plasma beta
//double ibetapl;  // inverse plasma beta
//double sigma;    // magnetization
//double eflem;    // energy flux, EM
//double eflma;    // energy flux, matter
//double lflem;    // angular momentum flux, EM
//double lflma;    // angular momentum flux, matter

// characterization.c
void check_box_limits(int i, int j, int k, int box_lower_i, int box_lower_j,
                      int box_lower_k, int box_upper_i, int box_upper_j,
                      int box_upper_k, int halfboxsize);
void find_normal();

// eigen.c
double *get_hessian2D(int i, int j);
double *get_hessian3D(int i, int j, int k);
double *get_evec(double *hessian);
double *get_evec2D(double hessian[4]);
double *get_evec3D(double hessian[9]);

// read_data.c
void assign_units();
void init_storage(void);
void init_storage_metric();
double *malloc_rank1(int n1, int size);
double **malloc_rank2_cont(int n1, int n2);
double ***malloc_rank3_cont(int n1, int n2, int n3);
void read_data(char *fname);
void read_gdump(char *fname);

// sheets.c
void check_adjacency (int i_adj, int j_adj, int k_adj, double J_peak, int mm,
                      int cells, int size_current, int *icoords, int *jcoords,
                      int *kcoords);
void get_current_sheets();
void get_locmax();
void write_current_sheets(char *fname);

// metric.c
void get_connection(double *X, double gcon[][NDIM], double conn[][NDIM][NDIM]);
void gcon_func(double gcov[][NDIM], double gcon[][NDIM]);
void gcov_func(double *X, double gcovp[][NDIM]);
double gdet_func(double gcov[][NDIM]);
void read_metric(char *fname);
void write_metric(char *fname);

void dxdxp_func(double *X, double dxdxp[][NDIM]);
void bl_coord(double *X, double *r, double *th, double *phi);
void bl_coord_vec(double *X, double *V);
void coord(int i, int j, int k, double *X);

void vofx_gammiecoords(double *X, double *V);
void vofx_sjetcoords(double *X, double *V);
double Ftr(double x);
double Ftrgenlin(double x, double xa, double xb, double ya, double yb);
double Fangle(double x);
double limlin(double x, double x0, double dx, double y0);
double mins(double f1, double f2, double df);
double maxs(double f1, double f2, double df);
double minmaxs(double f1, double f2, double df, double dir);
void vofx_cylindrified(double *Xin, void (*vofx)(double*, double*), double *Vout);

int invert_matrix(double Am[][NDIM], double Aminv[][NDIM]);
int LU_decompose(double A[][NDIM], int permute[]);
void LU_substitution(double A[][NDIM], double B[], int permute[]);

double zbrent(double (*func)(double, double), double param1, double lower, double upper, double tol);
double x1_ks_to_x1_mks (double *X_ks);
double x2_ks_to_x2_mks (double *X_ks);
double find_x1_mks(double x1, double radius);
double find_x2_mks(double x2, double theta);
void x_cyl_to_ks(double *X_cyl, double *X_ks);
void x_ks_to_mks(double *X_ks, double *X_mks);
void jac_cyl_to_ks(double *X_cyl, double jac[][NDIM]);
void jac_ks_to_mks(double *X_ks, double jac[][NDIM]);
void fourvec_old_to_new(double *fourvec_old, double jac[][NDIM], double *fourvec_new);
