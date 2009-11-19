/* Advection_diffusion.c */
void advection_diffusion_parameters(struct All_variables *E);
void advection_diffusion_allocate_memory(struct All_variables *E);
void PG_timestep_particle(struct All_variables *E);
void PG_timestep(struct All_variables *E);
void predictor(struct All_variables *E, float *field, float *fielddot);
void corrector(struct All_variables *E, float *field, float *fielddot, float *Dfielddot);
void pg_solver(struct All_variables *E, float *T, float *Tdot, float *DTdot, float **V, struct SOURCES Q0, float diff, int bc, float **TBC, unsigned int *FLAGS);
void pg_shape_fn(struct All_variables *E, int el, struct Shape_function *PG, float **V, double rtf[4][9], float diffusion);
void element_residual(struct All_variables *E, int el, struct Shape_function PG, float **vel, float *field, float *fielddot, struct SOURCES Q0, double Eres[9], double rtf[4][9], float diff, float **BC, unsigned int *FLAGS);
void std_timestep(struct All_variables *E);
void process_heating(struct All_variables *E);
/* Boundary_conditions.c */
void velocity_boundary_conditions(struct All_variables *E);
void temperature_boundary_conditions(struct All_variables *E);
void velocity_refl_vert_bc(struct All_variables *E);
void temperature_refl_vert_bc(struct All_variables *E);
void temperature_imposed_botm_bcs(struct All_variables *E, float *BC[], int dirn);
void horizontal_bc(struct All_variables *E, float *BC[], int ROW, int dirn, float value, unsigned int mask, char onoff, int level);
void velocity_apply_periodic_bcs(struct All_variables *E);
void temperature_apply_periodic_bcs(struct All_variables *E);
void strip_bcs_from_residual(struct All_variables *E, double *Res, int level);
void temperatures_conform_bcs(struct All_variables *E);
void velocities_conform_bcs(struct All_variables *E, double *U);
void equalize_id_ien_lm(struct All_variables *E);
/* Citcom.c */
int main(int argc, char **argv);
/* Composition_adv.c */
void Runge_Kutta(struct All_variables *E, float *C, float *V[4], int on_off);
void Euler(struct All_variables *E, float *C, float *V[4], int on_off);
void transfer_markers_processors(struct All_variables *E, int on_off);
void unify_markers_array(struct All_variables *E, int no_tran, int no_recv);
void prepare_transfer_arrays(struct All_variables *E);
int locate_processor(struct All_variables *E, double XMC1, double XMC2, double XMC3);
void get_C_from_markers(struct All_variables *E, float *C);
void element_markers(struct All_variables *E, int con);
void velocity_markers(struct All_variables *E, float *V[4], int con);
int get_element(struct All_variables *E, double XMC1, double XMC2, double XMC3, double dX[4]);
int in_the_domain(struct All_variables *E, double r, double t, double f);
float area_of_4node1(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4);
float area_of_3node(float x1, float y1, float x2, float y2, float x3, float y3);
float mean_of_3node(int a, float x1, float y1, float x2, float y2, float x3, float y3);
float mean_of_4node(int a, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4);
float mean_of_5node(int a, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, float x5, float y5);
float dist1(float XO[4], float XN[4]);
/* Construct_arrays.c */
void construct_ien(struct All_variables *E);
void construct_id(struct All_variables *E);
void construct_lm(struct All_variables *E);
void construct_node_maps(struct All_variables *E);
void construct_node_ks(struct All_variables *E);
void construct_masks(struct All_variables *E);
void construct_sub_element(struct All_variables *E);
void construct_elt_ks(struct All_variables *E);
void construct_elt_gs(struct All_variables *E);
void construct_mat_group(struct All_variables *E);
void construct_stiffness_B_matrix(struct All_variables *E);
void rebuild_BI_on_boundary(struct All_variables *E);
/* Convection.c */
void set_convection_defaults(struct All_variables *E);
void read_convection_settings(struct All_variables *E);
void convection_derived_values(struct All_variables *E);
void convection_allocate_memory(struct All_variables *E);
void convection_initial_fields(struct All_variables *E);
void convection_boundary_conditions(struct All_variables *E);
void convection_initial_temperature(struct All_variables *E);
void process_restart_tc(struct All_variables *E);
void convection_initial_markers1(struct All_variables *E);
void convection_initial_markers(struct All_variables *E, int use_element_nodes_for_init_c);
void setup_plume_problem(struct All_variables *E);
void PG_process(struct All_variables *E, int ii);
/* Drive_solvers.c */
void general_stokes_solver(struct All_variables *E);
/* Element_calculations.c */
void assemble_forces(struct All_variables *E, int penalty);
void get_elt_k(struct All_variables *E, int el, double elt_k[24 * 24], int lev, int iconv);
void assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs);
void e_assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs);
void n_assemble_del2_u(struct All_variables *E, double *u, double *Au, int level, int strip_bcs);
void build_diagonal_of_K(struct All_variables *E, int el, double elt_k[24 * 24], int level);
void build_diagonal_of_Ahat(struct All_variables *E);
void assemble_div_u(struct All_variables *E, double *U, double *divU, int level);
void assemble_grad_p(struct All_variables *E, double *P, double *gradP, int lev);
double assemble_dAhatp_entry(struct All_variables *E, int e, int level);
void get_elt_g(struct All_variables *E, int el, higher_precision elt_del[24][1], int lev);
void get_elt_h(struct All_variables *E, int el, double elt_h[1], int penalty);
void get_elt_f(struct All_variables *E, int el, double elt_f[24], int penalty, int bcs);
void get_aug_k(struct All_variables *E, int el, double elt_k[24 * 24], int level);
/* General_matrix_functions.c */
double **dmatrix(int nrl, int nrh, int ncl, int nch);
float **fmatrix(int nrl, int nrh, int ncl, int nch);
void dfree_matrix(double **m, int nrl, int nrh, int ncl, int nch);
void ffree_matrix(float **m, int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
float *fvector(int nl, int nh);
void dfree_vector(double *v, int nl, int nh);
void ffree_vector(float *v, int nl, int nh);
int *sivector(int nl, int nh);
void sifree_vector(int *v, int nl, int nh);
double pdot(struct All_variables *E, double *A, double *B, int lev);
double pselfdot(struct All_variables *E, double *A);
double vdot(struct All_variables *E, double *A, double *B, int level);
double vselfdot(struct All_variables *E, double *A, int level);
double vfselfdot(struct All_variables *E, float *A, int level);
float fdot(float *A, float *B, int n1, int n2);
float fselfdot(float *A, int n1, int n2);
float dot(struct All_variables *E, float *A, float *B);
float selfdot(struct All_variables *E, float *A);
void dvcopy(double *A, double *B, int a, int b);
void vcopy(float *A, float *B, int a, int b);
void vprod(double *R, double *A, double *B, int a, int b);
float fnmax(struct All_variables *E, float *A, int a, int b);
int solve_del2_u(struct All_variables *E, double *d0, double *F, double acc, int high_lev, int ic);
double multi_grid(struct All_variables *E, double *d1, double *F, double *Au, double acc, int hl);
double conj_grad(struct All_variables *E, double *d0, double *F, double *Au, double acc, int *cycles, int level);
void jacobi(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess);
void element_gauss_seidel(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess);
void gauss_seidel1(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess);
void gauss_seidel(struct All_variables *E, double *d0, double *F, double *Ad, double acc, int *cycles, int level, int guess);
void print_elt_k(struct All_variables *E, double a[24 * 24]);
double cofactor(double A[4][4], int i, int j, int n);
double determinant(double A[4][4], int n);
float area_of_4node(float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4);
double modified_plgndr_a(int l, int m, double t);
double sqrt_multis(int jj, int ii);
double multis(int ii);
int int_multis(int ii);
double plgndr_a(int l, int m, double t);
double sphere_h(int l, int m, double t, double f, int ic);
/* Geometry_cartesian.c */
void set_2dc_defaults(struct All_variables *E);
void set_2pt5dc_defaults(struct All_variables *E);
void set_3ds_defaults(struct All_variables *E);
void set_3dc_defaults(struct All_variables *E);
/* Ggrd_handling.c */
void convection_initial_temperature_ggrd(struct All_variables *E);
/* Global_operations.c */
void remove_horiz_ave(struct All_variables *E, float *X, float *H, int store_or_not);
void return_horiz_sum(struct All_variables *E, float *X, float *H, int nn);
void return_horiz_ave(struct All_variables *E, float *X, float *H);
float return_bulk_value(struct All_variables *E, float *Z, float z_thld, int average);
double global_vdot(struct All_variables *E, double *A, double *B, int lev);
double global_pdot(struct All_variables *E, double *A, double *B, int lev);
float global_tdot(struct All_variables *E, float *A, float *B, int lev);
float Tmax(struct All_variables *E, float *T);
float global_fmin(struct All_variables *E, float a);
float global_fmax(struct All_variables *E, float a);
void sum_across_depth_sph1(struct All_variables *E, float *sphc, float *sphs);
void sum_across_surface(struct All_variables *E, float *data, int total);
void sum_across_surf_sph1(struct All_variables *E, float *sphc, float *sphs);
void gather_TG_to_me0(struct All_variables *E, float *TG);
void propogator_down_process(struct All_variables *E, float *Tadi);
double sum_across_depth(struct All_variables *E, double temp1);
/* Instructions.c */
void read_instructions(struct All_variables *E, int argc, char **argv);
void allocate_common_vars(struct All_variables *E);
void interruption(int signal_number);
void global_default_values(struct All_variables *E);
void global_derived_values(struct All_variables *E);
void read_initial_settings(struct All_variables *E);
void check_bc_consistency(struct All_variables *E);
void set_up_nonmg_aliases(struct All_variables *E);
void report(struct All_variables *E, char *string);
void record(struct All_variables *E, char *string);
void common_initial_fields(struct All_variables *E);
void initial_pressure(struct All_variables *E);
void initial_velocity(struct All_variables *E);
/* Nodal_mesh.c */
void node_locations(struct All_variables *E);
void pre_interpolation(struct All_variables *E);
void dlogical_mesh_to_real(struct All_variables *E, double *data, int level);
void flogical_mesh_to_real(struct All_variables *E, float *data, int level);
void p_to_nodes(struct All_variables *E, double *P, float *PN, int lev);
void p_to_centres(struct All_variables *E, float *PN, double *P, int lev);
void v_to_intpts(struct All_variables *E, float *VN, float *VE, int lev);
void v_to_nodes(struct All_variables *E, float *VE, float *VN, int lev);
void visc_to_intpts(struct All_variables *E, float *VN, float *VE, int lev);
void visc_to_nodes(struct All_variables *E, float *VE, float *VN, int lev);
void visc_from_ele_to_gint(struct All_variables *E, float *VN, float *VE, int lev);
void visc_from_gint_to_ele(struct All_variables *E, float *VE, float *VN, int lev);
void visc_from_gint_to_nodes(struct All_variables *E, float *VE, float *VN, int lev);
void visc_from_nodes_to_gint(struct All_variables *E, float *VN, float *VE, int lev);
/* Output.c */
void output_velo_related(struct All_variables *E, int file_number);
void output_velo_related_binary(struct All_variables *E, int file_number);
void output_temp(struct All_variables *E, int file_number);
void process_restart(struct All_variables *E);
void print_field_spectral_regular(struct All_variables *E, float *TG, float *sphc, float *sphs, int proc_loc, char *filen);
/* Output_gzdir.c */
void output_velo_related_gzdir(struct All_variables *E, int file_number);
void process_restart_tc_gzdir(struct All_variables *E);
/* Pan_problem_misc_functions.c */
int get_process_identifier(void);
void unique_copy_file(struct All_variables *E, char *name, char *comment);
void thermal_buoyancy(struct All_variables *E);
double SIN_D(double x);
double COT_D(double x);
void *Malloc1(int bytes, char *file, int line);
int read_previous_field(struct All_variables *E, float *field, char *name, char *abbr);
void fcopy_interpolating(struct All_variables *E, float *X, float *Z, float *Y, int nx, int nz, int ny, float *T, float *TT);
float cross2d(float x11, float x12, float x21, float x22, int D);
void field_arbitrary_rectangle_file(struct All_variables *E, int parse_and_apply, struct Rect *RECT, char *name, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
void field_arbitrary_rectangle(struct All_variables *E, struct Rect *RECT, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
void field_arbitrary_circle_file(struct All_variables *E, int parse_and_apply, struct Circ *CIRC, char *name, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
void field_arbitrary_circle(struct All_variables *E, struct Circ *CIRC, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
void field_arbitrary_harmonic_file(struct All_variables *E, int parse_and_apply, struct Harm *HARM, char *name, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
void field_arbitrary_harmonic(struct All_variables *E, struct Harm *HARM, float *field, int BC, unsigned int *bcbitf, unsigned int bcmask_on, unsigned int bcmask_off);
double myatan(double y, double x);
FILE *safe_fopen(char *name, char *mode);
void *safe_malloc(size_t size);
void calc_cbase_at_tp(float theta, float phi, float *base);
void convert_pvec_to_cvec(float vr, float vt, float vp, float *base, float *cvec);
void rtp2xyz(float r, float theta, float phi, float *xout);
void myerror(char *message, struct All_variables *E);
/* Parallel_related.c */
void parallel_process_initilization(struct All_variables *E, int argc, char **argv);
void parallel_process_termination(void);
void parallel_domain_decomp1(struct All_variables *E);
void parallel_shuffle_ele_and_id(struct All_variables *E);
void parallel_shuffle_ele_and_id_bc1(struct All_variables *E);
void parallel_shuffle_ele_and_id_bc2(struct All_variables *E);
void parallel_communication_routs(struct All_variables *E);
void parallel_communication_routs1(struct All_variables *E);
void parallel_communication_routs2(struct All_variables *E);
void parallel_communication_routs3(struct All_variables *E);
void parallel_communication_routs4(struct All_variables *E);
void exchange_number_rec_markers(struct All_variables *E);
void exchange_markers(struct All_variables *E);
void exchange_id_d20(struct All_variables *E, double *U, int lev);
void exchange_node_f20(struct All_variables *E, float *U, int lev);
double CPU_time0(void);
void parallel_process_sync(void);
/* Parsing.c */
void setup_parser(struct All_variables *E, char *filename);
void shutdown_parser(struct All_variables *E);
void add_to_parameter_list(char *name, char *value);
int compute_parameter_hash_table(char *s);
int input_int(char *name, int *value, char *interpret, int m);
int input_string(char *name, char *value, char *Default, int m);
int input_boolean(char *name, int *value, char *interpret, int m);
int input_float(char *name, float *value, char *interpret, int m);
int input_double(char *name, double *value, char *interpret, int m);
int input_int_vector(char *name, int number, int *value, int m);
int input_char_vector(char *name, int number, char *value, int m);
int input_float_vector(char *name, int number, float *value, int m);
int input_double_vector(char *name, int number, double *value, int m);
int interpret_control_string(char *interpret, int *essential, double *Default, double *minvalue, double *maxvalue);
/* Phase_change.c */
void phase_change(struct All_variables *E, float *B6, float *B_b6, float *B4, float *B_b4);
/* Process_buoyancy.c */
void process_temp_field(struct All_variables *E, int ii);
void heat_flux(struct All_variables *E);
void heat_flux1(struct All_variables *E);
void plume_buoyancy_flux(struct All_variables *E);
/* Process_velocity.c */
void process_new_velocity(struct All_variables *E, int ii);
void get_surface_velo(struct All_variables *E, float *SV);
void get_ele_visc(struct All_variables *E, float *EV);
void get_surf_stress(struct All_variables *E, float *SXX, float *SYY, float *SZZ, float *SXY, float *SXZ, float *SZY);
void averages(struct All_variables *E);
/* Profiling.c */
float CPU_time(void);
/* Shape_functions.c */
void construct_shape_functions(struct All_variables *E);
double lpoly(int p, double y);
double lpolydash(int p, double y);
/* Size_does_matter.c */
void twiddle_thumbs(struct All_variables *yawn, int scratch_groin);
void get_global_shape_fn(struct All_variables *E, int el, int pressure, double rtf[4][9], int sphere, int level);
void form_rtf_bc(int k, double x[4], double rtf[4][9], double bc[4][4]);
void get_rtf(struct All_variables *E, int el, int pressure, double rtf[4][9], int lev);
void construct_c3x3matrix_el(struct All_variables *E, int el, struct CC *cc, struct CCX *ccx, int lev, int pressure);
void get_global_1d_shape_fn(struct All_variables *E, int el, struct Shape_function1 *GM, struct Shape_function1_dA *dGammax, int top);
void get_global_1d_shape_fn1(struct All_variables *E, int el, struct Shape_function1 *GM, struct Shape_function1_dA *dGammax, int top);
void mass_matrix(struct All_variables *E);
/* Solver_conj_grad.c */
void set_cg_defaults(struct All_variables *E);
void cg_allocate_vars(struct All_variables *E);
void assemble_forces_iterative(struct All_variables *E);
/* Solver_multigrid.c */
void set_mg_defaults(struct All_variables *E);
void mg_allocate_vars(struct All_variables *E);
void project_vector(struct All_variables *E, int start_lev, double *AU, double *AD, int ic);
void interp_vector(struct All_variables *E, int start_lev, double *AD, double *AU);
void project_scalar_e(struct All_variables *E, int start_lev, float *AU, float *AD);
void project_scalar(struct All_variables *E, int start_lev, float *AU, float *AD);
void project_viscosity(struct All_variables *E);
void inject_node_fvector(struct All_variables *E, int start_lev, float **AU, float **AD);
void inject(struct All_variables *E, int start_lev, double *AU, double *AD);
void un_inject_vector(struct All_variables *E, int start_lev, double *AD, double *AU);
void inject_scalar(struct All_variables *E, int start_lev, float *AU, float *AD);
void inject_scalar_e(struct All_variables *E, int start_lev, float *AU, float *AD);
/* Sphere_harmonics.c */
void set_sphere_harmonics(struct All_variables *E);
void sphere_harmonics_layer(struct All_variables *E, float **T, float *sphc, float *sphs, int iprint, char *filen);
void sphere_interpolate(struct All_variables *E, float **T, float *TG);
void sphere_expansion(struct All_variables *E, float *TG, float *sphc, float *sphs);
void inv_sphere_harmonics(struct All_variables *E, float *sphc, float *sphs, float *TG, int proc_loc);
void compute_sphereh_table(struct All_variables *E);
/* Stokes_flow_Incomp.c */
void solve_constrained_flow_iterative(struct All_variables *E);
float solve_Ahat_p_fhat(struct All_variables *E, double *V, double *P, double *F, double imp, int *steps_max);
void v_from_vector(struct All_variables *E, float **V, double *F);
/* Topo_gravity.c */
void get_CBF_topo(struct All_variables *E, float *H, float *HB);
void get_STD_topo(struct All_variables *E, float *tpg, float *tpgb, int ii);
/* Viscosity_structures.c */
void viscosity_parameters(struct All_variables *E);
void get_viscosity_option(struct All_variables *E);
void viscosity_for_system(struct All_variables *E);
void get_system_viscosity(struct All_variables *E, int propogate, float *evisc, float *visc);
void apply_viscosity_smoother(struct All_variables *E, float *visc, float *evisc);
void visc_from_mat(struct All_variables *E, float *Eta, float *EEta);
void visc_from_T(struct All_variables *E, float *Eta, float *EEta, int propogate);
void visc_from_S(struct All_variables *E, float *Eta, float *EEta, int propogate);
void strain_rate_2_inv(struct All_variables *E, float *EEDOT, int SQRT);
int layers(struct All_variables *E, float x3);
int weak_zones(struct All_variables *E, int node, float t_b);
float boundary_thickness(struct All_variables *E, float *H);
