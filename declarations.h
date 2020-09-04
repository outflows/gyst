// Constants, definitions and units
double M_unit, L_unit, T_unit, RHO_unit, U_unit, B_unit, Ne_unit;

char dump_file[100], jcs_output[100], jchar_output[100], jpeak_output[100];
char gdump_file[100];

double *****grid_gcov, *****grid_gcon, ***grid_gdet, ******grid_conn;

struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
};

double x1in;
int peak_count;

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
double stopx[NDIM];
double dx[NDIM];
double dt;

double ***dJdx, ***dJdy, ***dJdz;

double ***J_cs;
double ***J_cs_peak;
double ***J_cs_char;
int ***flag;

int N1;
int N2;
int N3;
double ***a_r;
double ***a_th;
double ***a_phi;
double ***rho;
double ug;
double ****V;
double ****B;
double divb;
double gdet;
double ****ucov;
double ****ucon;
double ****bcov;
double ****bcon;
double ****jcov;
double ****jcon;

double jsq;
double jdotu;
double Jsq;
double gJsq;
double ***J;        // current
double ***sigma; // phi-component of magnetization
double ***sigmaphi; // phi-component of magnetization
double ***Sigmaphi; // phi-component of magnetization
double pg;        // gas pressure
//double K;
double EF;       // enthalpy
double bsq;
double EE;       // eta
double va2;      // relativistic alfvén speed squared
double ***va;      // relativistic alfvén speed
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
double ***bernoulli; // bernoulli


// characterization.c


// eigen.c
void get_hessian(int i, int j, int k, double Hess[3][3]);
void get_evec(double Hess[3][3], double eigvec[9]);
void get_1stderivative(int i, int j, int k, double ***J_cs);


void characterize2D();
void get_1stderivative2D(int i, int j, int k, double ***J_cs);
void get_evec2D(double Hess[2][2], double eigvec[4]);
void get_hessian2D(int i, int j, int k, double Hess2D[2][2]);
void xy_to_rth(double x, double y, double *r, double *th);
void rth_to_xy(double r, double th, double *x, double *y);
void vec_pol_to_cart(double eig_pol[2], double eig_cart[2], double theta);

// read_data.c
void assign_units();
void read_data(char *fname);
void read_gdump(char *fname);
void read_current_sheets(char *fname, double ***sheets);
void write_current_sheets(char *fname, double ***sheets);

// malloc.c
void init_storage_data();
void init_storage_metric();
double *malloc_rank1(int n1, int size);
double **malloc_rank2_cont(int n1, int n2);
double ***malloc_rank3_cont(int n1, int n2, int n3);
double ****malloc_rank4_cont(int n1, int n2, int n3, int n4);
double *****malloc_rank5_cont(int n1, int n2, int n3, int n4, int n5);
double ******malloc_rank6_cont(int n1, int n2, int n3, int n4, int n5, int n6);

// sheets.c
void get_current_sheets();
void fill(int i, int j, int k, double threshold);
void get_locmax();

// characterization.c
void characterize();
void xyz_to_rthphi(double x, double y, double z, double *r, double *th, double *phi);
void rthphi_to_xyz(double r, double th, double phi, double *x, double *y, double *z);
void vec_sph_to_cart(double eig_sph[3], double eig_cart[3], double theta, double phi);
void reverse_array(double *array, int n);
void merge_arrays(double *a1, double *a2, double *newa, int size_a1, int size_a2);

// metric.c
void get_connection(double *X, struct of_geom *geom, double ******conn, int i, int j, int k);
void get_geometry(int ii, int jj, int kk, struct of_geom *geom);
void gcon_func(double gcov[][NDIM], double gcon[][NDIM]);
void gcov_func(double *X, double gcovp[][NDIM]);
double gdet_func(double gcov[][NDIM]);
void read_metric(char *fname);
void write_metric(char *fname);
void dxdxp_func(double *X, double dxdxp[][NDIM]);
void bl_coord(double *X, double *r, double *th, double *phi);
void bl_coord_vec(double *X, double *V);
void coord(int i, int j, int k, double *X);
void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
double dot3(double *v1, double *v2);

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
double find_x1_cyl(double x1, double radius);
double find_x2_cyl(double x2, double theta);
double calcrmks(double x1);
void to1stquadrant_single(double x2in, double x2out, int *ismirrored);
double func2_single(double r0, double rr, double x20, double x2);
double sinth1in_single(double r0, double rr, double x20, double x2);
double th2in_single(double r0, double rr, double x20, double x2);
double sinth0_single(double x20, double r0, double rr);
double calcth_cylindrified(double x2in);

//double x1_ks_to_x1_mks (double *X_ks);
//double x2_ks_to_x2_mks (double *X_ks);
//double find_x1_mks(double x1, double radius);
//double find_x2_mks(double x2, double theta);
//void x_cyl_to_ks(double *X_cyl, double *X_ks);
//void x_ks_to_mks(double *X_ks, double *X_mks);
//void jac_cyl_to_ks(double *X_cyl, double jac[][NDIM]);
//void jac_ks_to_mks(double *X_ks, double jac[][NDIM]);
//void fourvec_old_to_new(double *fourvec_old, double jac[][NDIM], double *fourvec_new);



void characterize2D_single();
//void characterize2D_single2();
void reverse_array_int(int *array, int n);
void merge_arrays_int(int *a1, int *a2, int *newa, int size_a1, int size_a2);
void Xtoijk_new(double X[NDIM], int *i, int *j, int *k, double del[NDIM]);
//void characterize2D2();





double find_normal(int i, int j, int k);
void check_box_limits(int i, int j, int k,
                      int box_lower_i, int box_lower_j, int box_lower_k,
                      int box_upper_i, int box_upper_j, int box_upper_k,
                      int halfboxsize);
void characterize2D_normal();
