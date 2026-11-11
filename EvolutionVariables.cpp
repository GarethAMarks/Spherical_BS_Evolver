#ifndef EVOLUTIONVARIABLES_CPP_
#define EVOLUTIONVARIABLES_CPP_

#ifndef SPACEDIM
#define SPACEDIM 3
#endif

#include "EvolutionVariables.h"
#include "mathutils.h"
#include "mpi_helpers.h"
#include "spline.h"
#include <iomanip>
#include<algorithm>
#include<filesystem>
#include <complex.h>
#include <array>
#include <chrono>

// reflection offset across z = 0: 0 for vertex-centred grids, 1 for cell-centred (set in Spacetime::initialize)
int BSSNSlice::origin_mirror_shift = 0;
bool BSSNSlice::one_sided_outer = false;
int BSSNSlice::outer_last_j = -1;
int BSSNSlice::forced_res_fac = 0;
int BSSNSlice::forced_below_j = -1;
const std::vector<int>* BSSNSlice::level_table = nullptr;

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::rotate;
// Spacetime constructor: initialize reusable CSF model with this context
Spacetime::Spacetime() : csf_model(this), rsf_model(this) {

    sub_min_time = 10.0;
    adaptive_timestep = false;
    A_thresh = 10000000.0; // Default to a very large value if not provided in params.
    // Initialize scalar-field energies
    E_phi = 0.0;
    E_psi = 0.0;

    // Cache rank/size up front so early error paths on any rank can print
    // consistently. local_N / global_lo / global_hi are populated later in
    // init_domain_decomposition() once n_gridpoints is known.
    mpi_rank = bss_mpi::rank();
    mpi_size = bss_mpi::size();
}

// Contiguously split n_gridpoints across mpi_size ranks. Remainder points
// are distributed one-per-rank starting at rank 0, so |local_N_r - local_N_s| <= 1.
// Precondition: n_gridpoints must already be set (initialize() takes care of this).
void Spacetime::init_domain_decomposition()
{
    if (n_gridpoints <= 0)
    {
        if (mpi_rank == 0)
            cerr << "ERROR: init_domain_decomposition called with n_gridpoints=" << n_gridpoints << endl;
        bss_mpi::abort_all(1);
    }

    // Require at least NGHOST+1 owned points per rank so a rank's owned zone is
    // never fully covered by its own ghost width (which would break the 7-point
    // KO stencil and any exchange that assumes a positive owned interior).
    const int min_per_rank = NGHOST + 1;
    if (n_gridpoints < mpi_size * min_per_rank)
    {
        if (mpi_rank == 0)
            cerr << "ERROR: n_gridpoints=" << n_gridpoints
                 << " too small for " << mpi_size << " MPI ranks "
                 << "(need at least " << mpi_size * min_per_rank << ")." << endl;
        bss_mpi::abort_all(1);
    }

    const int base = n_gridpoints / mpi_size;
    const int rem  = n_gridpoints % mpi_size;

    local_N   = base + (mpi_rank < rem ? 1 : 0);
    global_lo = mpi_rank * base + std::min(mpi_rank, rem);
    global_hi = global_lo + local_N;

    if (mpi_rank == 0 && mpi_size > 1)
        cout << "MPI domain decomposition: " << mpi_size << " ranks, "
             << n_gridpoints << " total gridpoints, ~" << base
             << " points/rank, NGHOST=" << NGHOST << endl;
}

#ifdef USE_MPI
// Compute (counts, displs) for MPI_Gatherv/Allgatherv given a per-point double
// stride. Uses the same base/rem split as init_domain_decomposition so the
// layout matches what each rank actually owns. (Only used by the MPI build.)
static void compute_gather_layout(int n_gridpoints, int mpi_size, int doubles_per_point,
                                  std::vector<int>& counts, std::vector<int>& displs)
{
    counts.assign(mpi_size, 0);
    displs.assign(mpi_size, 0);
    const int base = n_gridpoints / mpi_size;
    const int rem  = n_gridpoints % mpi_size;
    for (int r = 0; r < mpi_size; ++r)
    {
        const int lN  = base + (r < rem ? 1 : 0);
        const int gLo = r * base + std::min(r, rem);
        counts[r] = lN  * doubles_per_point;
        displs[r] = gLo * doubles_per_point;
    }
}
#endif

void Spacetime::gather_owned_to_root(const std::vector<double>& local_owned,
                                     std::vector<double>& full_buf,
                                     int doubles_per_point) const
{
#ifdef USE_MPI
    if (mpi_size <= 1) { full_buf = local_owned; return; }
    std::vector<int> counts, displs;
    compute_gather_layout(n_gridpoints, mpi_size, doubles_per_point, counts, displs);
    double* recv_ptr = nullptr;
    if (mpi_rank == 0)
    {
        full_buf.assign(n_gridpoints * doubles_per_point, 0.0);
        recv_ptr = full_buf.data();
    }
    MPI_Gatherv(local_owned.data(), local_N * doubles_per_point, MPI_DOUBLE,
                recv_ptr, counts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
#else
    (void)doubles_per_point;
    full_buf = local_owned;
#endif
}

void Spacetime::allgather_owned(const std::vector<double>& local_owned,
                                std::vector<double>& full_buf,
                                int doubles_per_point) const
{
#ifdef USE_MPI
    if (mpi_size <= 1) { full_buf = local_owned; return; }
    std::vector<int> counts, displs;
    compute_gather_layout(n_gridpoints, mpi_size, doubles_per_point, counts, displs);
    full_buf.assign(n_gridpoints * doubles_per_point, 0.0);
    MPI_Allgatherv(local_owned.data(), local_N * doubles_per_point, MPI_DOUBLE,
                   full_buf.data(), counts.data(), displs.data(), MPI_DOUBLE,
                   MPI_COMM_WORLD);
#else
    (void)doubles_per_point;
    full_buf = local_owned;
#endif
}

// Serialize the caller-provided slice's state (and theta) at each rank's owned
// range, gather to rank 0, and on rank 0 write via a scratch BSSNSlice so the
// input slice is not mutated anywhere. Non-root ranks participate in the
// gather and then return.
void Spacetime::write_checkpoint(BSSNSlice& slice, const std::string& file_name)
{
    // 16 doubles per point (theta as 0 when not CCZ4 — but we handle theta
    // as a separate gather below when CCZ4 is on).
    constexpr int NCOL = 16;
    std::vector<double> local(local_N * NCOL);
    for (int k = 0; k < local_N; ++k)
    {
        const int j = owned_begin() + k;
        double* p = &local[k * NCOL];
        p[0]  = slice.states2[j].bssn.chi;
        p[1]  = slice.states2[j].bssn.h_zz;
        p[2]  = slice.states2[j].bssn.h_ww;
        p[3]  = slice.states2[j].bssn.A_zz;
        p[4]  = slice.states2[j].bssn.A_ww;
        p[5]  = slice.states2[j].bssn.K;
        p[6]  = slice.states2[j].bssn.c_chris_Z;
        p[7]  = slice.states2[j].csf.phi_re;
        p[8]  = slice.states2[j].csf.phi_im;
        p[9]  = slice.states2[j].csf.K_phi_re;
        p[10] = slice.states2[j].csf.K_phi_im;
        p[11] = slice.states2[j].bssn.alpha;
        p[12] = slice.states2[j].bssn.beta;
        p[13] = use_CCZ4 ? slice.theta[j] : 0.0;
        p[14] = slice.states2[j].rsf.psi;
        p[15] = slice.states2[j].rsf.K_psi;
    }

    std::vector<double> full;
    gather_owned_to_root(local, full, NCOL);

    if (mpi_rank != 0) return;

    // Build a scratch slice with the gathered data and write via the
    // existing BSSNSlice::write_slice formatter to keep the on-disk format
    // identical to the serial case (so read_checkpoint stays interoperable).
    BSSNSlice scratch;
    scratch.R = slice.R;
    scratch.use_CCZ4 = use_CCZ4;
    scratch.states2.resize(n_gridpoints);
    if (use_CCZ4) scratch.theta.resize(n_gridpoints);
    for (int j = 0; j < n_gridpoints; ++j)
    {
        const double* p = &full[j * NCOL];
        scratch.states2[j].bssn.chi       = p[0];
        scratch.states2[j].bssn.h_zz      = p[1];
        scratch.states2[j].bssn.h_ww      = p[2];
        scratch.states2[j].bssn.A_zz      = p[3];
        scratch.states2[j].bssn.A_ww      = p[4];
        scratch.states2[j].bssn.K         = p[5];
        scratch.states2[j].bssn.c_chris_Z = p[6];
        scratch.states2[j].csf.phi_re     = p[7];
        scratch.states2[j].csf.phi_im     = p[8];
        scratch.states2[j].csf.K_phi_re   = p[9];
        scratch.states2[j].csf.K_phi_im   = p[10];
        scratch.states2[j].bssn.alpha     = p[11];
        scratch.states2[j].bssn.beta      = p[12];
        if (use_CCZ4) scratch.theta[j]    = p[13];
        scratch.states2[j].rsf.psi        = p[14];
        scratch.states2[j].rsf.K_psi      = p[15];
    }
    scratch.write_slice(file_name);
}

// Fill NGHOST ghost slots on each side of this rank's owned range with data
// from neighboring ranks. Under the replicated-storage model slice.states2
// (and slice.theta if CCZ4) is sized n_gridpoints on every rank; each rank
// only *writes* [global_lo, global_hi) but the array is indexed globally, so
// the neighbor's data lands at global positions:
//   left ghosts  = [global_lo - NGHOST, global_lo)   from left  neighbor
//   right ghosts = [global_hi, global_hi + NGHOST)   from right neighbor
//
// Wire format per point: NDOF_STATE consecutive doubles laid out as
//   [BSSNState (9)] [CSF (4)] [RSF (2)]
// matching the struct member ordering below. Theta is a separate double per
// point sent in a second message when use_CCZ4 is on.
//
// In serial (no USE_MPI, or mpi_size==1) this is a no-op: owns_origin() and
// owns_outer() are both true so nothing is sent.
void Spacetime::exchange_ghosts(BSSNSlice& slice)
{
#ifdef USE_MPI
    if (mpi_size <= 1)
        return;

    // Sanity: the caller must have arrays sized to the full grid.
    if ((int)slice.states2.size() < n_gridpoints)
    {
        cerr << "ERROR (rank " << mpi_rank << "): exchange_ghosts on slice with "
             << slice.states2.size() << " states2 entries; expected >= "
             << n_gridpoints << endl;
        bss_mpi::abort_all(1);
    }
    if (use_CCZ4 && (int)slice.theta.size() < n_gridpoints)
    {
        cerr << "ERROR (rank " << mpi_rank << "): exchange_ghosts on slice with "
             << slice.theta.size() << " theta entries; expected >= "
             << n_gridpoints << endl;
        bss_mpi::abort_all(1);
    }

    const int nbuf = NGHOST * NDOF_STATE;
    std::vector<double> send_left (nbuf);
    std::vector<double> send_right(nbuf);
    std::vector<double> recv_left (nbuf);
    std::vector<double> recv_right(nbuf);

    // Pack the NGHOST owned points nearest each interior boundary.
    //   Leftmost  owned globals: [global_lo, global_lo + NGHOST)         -> send to left  neighbor
    //   Rightmost owned globals: [global_hi - NGHOST, global_hi)         -> send to right neighbor
    auto pack = [&](std::vector<double>& buf, int global_start)
    {
        for (int k = 0; k < NGHOST; ++k)
        {
            const State& s = slice.states2[global_start + k];
            double* p = &buf[k * NDOF_STATE];
            p[0]  = s.bssn.chi;
            p[1]  = s.bssn.h_zz;
            p[2]  = s.bssn.h_ww;
            p[3]  = s.bssn.A_zz;
            p[4]  = s.bssn.A_ww;
            p[5]  = s.bssn.K;
            p[6]  = s.bssn.c_chris_Z;
            p[7]  = s.bssn.alpha;
            p[8]  = s.bssn.beta;
            p[9]  = s.csf.phi_re;
            p[10] = s.csf.phi_im;
            p[11] = s.csf.K_phi_re;
            p[12] = s.csf.K_phi_im;
            p[13] = s.rsf.psi;
            p[14] = s.rsf.K_psi;
        }
    };
    auto unpack = [&](const std::vector<double>& buf, int global_start)
    {
        for (int k = 0; k < NGHOST; ++k)
        {
            State& s = slice.states2[global_start + k];
            const double* p = &buf[k * NDOF_STATE];
            s.bssn.chi       = p[0];
            s.bssn.h_zz      = p[1];
            s.bssn.h_ww      = p[2];
            s.bssn.A_zz      = p[3];
            s.bssn.A_ww      = p[4];
            s.bssn.K         = p[5];
            s.bssn.c_chris_Z = p[6];
            s.bssn.alpha     = p[7];
            s.bssn.beta      = p[8];
            s.csf.phi_re     = p[9];
            s.csf.phi_im     = p[10];
            s.csf.K_phi_re   = p[11];
            s.csf.K_phi_im   = p[12];
            s.rsf.psi        = p[13];
            s.rsf.K_psi      = p[14];
        }
    };

    if (!owns_origin()) pack(send_left,  global_lo);
    if (!owns_outer())  pack(send_right, global_hi - NGHOST);

    const int left_nbr  = mpi_rank - 1;
    const int right_nbr = mpi_rank + 1;
    const int tag_l2r = 100;
    const int tag_r2l = 101;

    MPI_Request reqs[4];
    int nreq = 0;
    if (!owns_origin())
    {
        MPI_Irecv(recv_left.data(),  nbuf, MPI_DOUBLE, left_nbr,  tag_l2r, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(send_left.data(),  nbuf, MPI_DOUBLE, left_nbr,  tag_r2l, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    if (!owns_outer())
    {
        MPI_Irecv(recv_right.data(), nbuf, MPI_DOUBLE, right_nbr, tag_r2l, MPI_COMM_WORLD, &reqs[nreq++]);
        MPI_Isend(send_right.data(), nbuf, MPI_DOUBLE, right_nbr, tag_l2r, MPI_COMM_WORLD, &reqs[nreq++]);
    }
    MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

    if (!owns_origin()) unpack(recv_left,  global_lo - NGHOST);
    if (!owns_outer())  unpack(recv_right, global_hi);

    // Theta (CCZ4 only) — separate small message.
    if (use_CCZ4)
    {
        std::vector<double> tsl(NGHOST), tsr(NGHOST), trl(NGHOST), trr(NGHOST);
        if (!owns_origin())
            for (int k = 0; k < NGHOST; ++k) tsl[k] = slice.theta[global_lo + k];
        if (!owns_outer())
            for (int k = 0; k < NGHOST; ++k) tsr[k] = slice.theta[global_hi - NGHOST + k];

        nreq = 0;
        if (!owns_origin())
        {
            MPI_Irecv(trl.data(), NGHOST, MPI_DOUBLE, left_nbr,  tag_l2r + 2, MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Isend(tsl.data(), NGHOST, MPI_DOUBLE, left_nbr,  tag_r2l + 2, MPI_COMM_WORLD, &reqs[nreq++]);
        }
        if (!owns_outer())
        {
            MPI_Irecv(trr.data(), NGHOST, MPI_DOUBLE, right_nbr, tag_r2l + 2, MPI_COMM_WORLD, &reqs[nreq++]);
            MPI_Isend(tsr.data(), NGHOST, MPI_DOUBLE, right_nbr, tag_l2r + 2, MPI_COMM_WORLD, &reqs[nreq++]);
        }
        MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);

        if (!owns_origin())
            for (int k = 0; k < NGHOST; ++k) slice.theta[global_lo - NGHOST + k] = trl[k];
        if (!owns_outer())
            for (int k = 0; k < NGHOST; ++k) slice.theta[global_hi + k] = trr[k];
    }
#else
    (void)slice;
#endif
}

//overloads for addition/ scalar multiplication of BSSNState sets and the slice arrays containing them
// Forward declarations for BSSNState operators used by State operators below
BSSNState operator+(const BSSNState& s1, const BSSNState& s2);
BSSNState operator-(const BSSNState& s1, const BSSNState& s2);
BSSNState operator*(const BSSNState& s1, const BSSNState& s2);
BSSNState operator*(double c, const BSSNState& s);

// Composite State elementwise operations (applies to both csf and bssn)
State operator+(const State& a, const State& b) {
    return State{a.csf + b.csf, a.bssn + b.bssn, a.rsf + b.rsf};
}

State operator-(const State& a, const State& b) {
    return State{a.csf - b.csf, a.bssn - b.bssn, a.rsf - b.rsf};
}

State operator*(const State& a, const State& b) {
    return State{a.csf * b.csf, a.bssn * b.bssn, a.rsf * b.rsf};
}

State operator*(double c, const State& a) {
    return State{c * a.csf, c * a.bssn, c * a.rsf};
}

inline State operator/(const State& a, double c) {
    return (1.0 / c) * a;
}

BSSNState operator+(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi + s2.chi, s1.h_zz + s2.h_zz, s1.h_ww + s2.h_ww, s1.A_zz + s2.A_zz, s1.A_ww + s2.A_ww, s1.K + s2.K, s1.c_chris_Z + s2.c_chris_Z,
                      s1.alpha + s2.alpha, s1.beta + s2.beta};
}

BSSNState operator-(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi - s2.chi, s1.h_zz - s2.h_zz, s1.h_ww - s2.h_ww, s1.A_zz - s2.A_zz, s1.A_ww - s2.A_ww, s1.K - s2.K, s1.c_chris_Z - s2.c_chris_Z,
                      s1.alpha - s2.alpha, s1.beta - s2.beta};
}

//termwise multiplication for convenient use with characteristic speeds etc.
BSSNState operator*(const BSSNState& s1, const BSSNState& s2)
{
    return (BSSNState){s1.chi * s2.chi, s1.h_zz * s2.h_zz, s1.h_ww * s2.h_ww, s1.A_zz * s2.A_zz, s1.A_ww * s2.A_ww, s1.K * s2.K, s1.c_chris_Z * s2.c_chris_Z,
                      s1.alpha * s2.alpha, s1.beta * s2.beta};
}

BSSNState operator*(double c, const BSSNState& s)
{
    return (BSSNState){c * s.chi, c * s.h_zz, c * s.h_ww , c * s.A_zz, c * s.A_ww, c * s.K, c * s.c_chris_Z,
                      c * s.alpha, c * s.beta};
}

inline BSSNState operator/(const BSSNState& s, double c)
{
    return (1. / c) * s;
}

BSSNSlice operator+(const BSSNSlice& slice1, const BSSNSlice& slice2)
{
    BSSNSlice return_slice;
    unsigned int length = slice1.states2.size();

    if (length != slice2.states2.size() || slice1.R != slice2.R)
    {
        cerr << "ERROR: attempted to add data from slices of different sizes!" << endl;
        exit(1);
    }

    return_slice.states2.resize(length);

    if (slice1.use_CCZ4)
        return_slice.theta.resize(length);


    return_slice.R = slice1.R;
    return_slice.has_BH = slice1.has_BH;
    return_slice.use_CCZ4 = slice1.use_CCZ4;
    return_slice.refinement_points = slice1.refinement_points;

    for (unsigned int j = 0; j < length; j++)
    {
        return_slice.states2[j] = slice1.states2[j] + slice2.states2[j];

        if (slice1.use_CCZ4)
            return_slice.theta[j] = slice1.theta[j] + slice2.theta[j];
    }

    return return_slice;
}

BSSNSlice operator*(double c, const BSSNSlice& slice)
{
    BSSNSlice return_slice;

    int length = slice.states2.size();
    return_slice.states2.resize(length);

    if (slice.use_CCZ4)
        return_slice.theta.resize(length);

    return_slice.R = slice.R;
    return_slice.has_BH = slice.has_BH;
    return_slice.refinement_points = slice.refinement_points;
    return_slice.use_CCZ4 = slice.use_CCZ4;

    for (int j = 0; j < length; j++)
    {
        return_slice.states2[j] = c * slice.states2[j];

        if (slice.use_CCZ4)
            return_slice.theta[j] = c * slice.theta[j];
    }

    return return_slice;
}

//radial partial derivative of a given bssn var at index. Order is an optional argument that allows higher z-derivatives to be taken, default is 1.
//should add checks on refinement levels
//Builds the 5-point stencil (indices, spacing, parity flips, outer one-sided position) at the given index
void BSSNSlice::build_stencil(int index, DerivStencil& st) const
{
    // Use composite state vector size for grid length to ensure compatibility with RK temp slices
    int n_gridpoints = states2.size();
    // Use cell-centred grid spacing to match read_BS_data (innermost point at z = 0.5*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    //refinement level; starts at 1 for no refinement and halves every time (table lookup when available)
    const int ref_level = (level_table && (int)level_table->size() == n_gridpoints && index >= 0 && index < n_gridpoints)
                          ? (*level_table)[index] : const_cast<BSSNSlice*>(this)->get_refinement_level(index, const_cast<std::vector<int>&>(refinement_points));

    // Use bit-shift instead of pow for 2^(ref_level-1)
    int res_fac = 1 << (ref_level - 1);
    if (forced_res_fac > 0 && index <= forced_below_j) //Berger-Oliger buffer points: coarse spacing inside the finer region
        res_fac = forced_res_fac;

    //check index is valid and error out if not
   if (index < 0 || index >= n_gridpoints )
    {
        cerr << "ERROR: invalid index requested in radial derivative" << endl;
        exit(1);
    }

    //set of indices using which derivative will be evaluated
    int* J = st.J;
    J[0] = index - 2 * res_fac; J[1] = index - res_fac; J[2] = index; J[3] = index + res_fac; J[4] = index + 2 * res_fac;

    //enforce BC at z = 0 by using symmetry: the ghost at index -k mirrors index k (vertex-centred, z_j = j dr)
    //or index k - 1 (cell-centred, z_j = (j + 1/2) dr)
    for (int k = 0; k < 2; k++)
        if (J[k] < 0) J[k] = -J[k] - origin_mirror_shift;

    //At outer edge just use central value for rightmost two; this is cut off in initialize()/ overwritten with BCs anyway.
    //With one_sided_outer (constraint-preserving BCs), use one-sided 5-point stencils on the outermost 5 points instead.
    st.outer_pos = 0; // 0: centred stencil; 4 / 5: index is the 2nd-outermost / outermost point of the one-sided stencil
    const int jN = (outer_last_j >= 0) ? outer_last_j : n_gridpoints - 1; // outermost active point (differs from n_gridpoints - 1 with refinement)
    if (one_sided_outer && index > jN - 2 * res_fac)
    {
        st.outer_pos = (index == jN) ? 5 : 4;
        for (int k = 0; k < 5; k++)
            J[k] = jN - (4 - k) * res_fac;
    }
    else
    {
    if (index >= n_gridpoints - res_fac)
        J[3] = J[2];
    if (index >= n_gridpoints - 2 * res_fac)
        J[4] = J[3];
    }

    //odd-parity variables (contracted Christoffel, beta) change sign across z = 0 at the mirrored indices
    st.flip0 = (index <= 1);
    st.flip1 = (index == 0);
    st.h = res_fac * dr;
    st.h2 = pow(st.h, 2);
}

//Derivative of a variable from a prebuilt stencil (same arithmetic as the original d_z)
double BSSNSlice::d_z_with(bssn_var var, const DerivStencil& st, int order) const
{
    const int* J = st.J;
    double F[5];

    //Fill out five values; note that at outer boundary this gives unreliable results
    switch (var)
    {
        case v_chi:      for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.chi;      break;
        case v_h_zz:     for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.h_zz;     break;
        case v_h_ww:     for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.h_ww;     break;
        case v_A_zz:     for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.A_zz;     break;
        case v_A_ww:     for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.A_ww;     break;
        case v_K:        for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.K;        break;
        case v_c_chris_Z:for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.c_chris_Z;break;
        case v_phi_re:   for (int j = 0; j < 5; j++) F[j] = states2[J[j]].csf.phi_re;    break;
        case v_phi_im:   for (int j = 0; j < 5; j++) F[j] = states2[J[j]].csf.phi_im;    break;
        case v_K_phi_re: for (int j = 0; j < 5; j++) F[j] = states2[J[j]].csf.K_phi_re;  break;
        case v_K_phi_im: for (int j = 0; j < 5; j++) F[j] = states2[J[j]].csf.K_phi_im;  break;
        case v_alpha:    for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.alpha;    break;
        case v_beta:     for (int j = 0; j < 5; j++) F[j] = states2[J[j]].bssn.beta;     break;
        case v_theta:    for (int j = 0; j < 5; j++) F[j] = theta[J[j]];                 break;
        case v_psi:      for (int j = 0; j < 5; j++) F[j] = states2[J[j]].rsf.psi;       break;
        case v_K_psi:    for (int j = 0; j < 5; j++) F[j] = states2[J[j]].rsf.K_psi;     break;
        default:
            cerr << "ERROR: invalid variable requested for differentiation" << endl;
            exit(0);
    }

    const bool parity_is_odd = ((var == v_c_chris_Z || var == v_beta)/*|| (has_BH && (var == v_chi || var == v_alpha))*/);

    //account for odd parity of contracted christoffel symbols and beta across z = 0
    if (st.flip0 && parity_is_odd)
        F[0] = -F[0];
    if (st.flip1 && parity_is_odd)
        F[1] = -F[1];

    if (st.outer_pos != 0)
    {
        const double h = st.h;
        if (order == 1)
            return (st.outer_pos == 5) ? p5_stencil(h, F[0], F[1], F[2], F[3], F[4]) : p4_stencil(h, F[0], F[1], F[2], F[3], F[4]);
        if (order == 2) // one-sided 5-point second derivatives
            return (st.outer_pos == 5) ? (11. * F[0] - 56. * F[1] + 114. * F[2] - 104. * F[3] + 35. * F[4]) / (12. * h * h)
                                       : (-F[0] + 4. * F[1] + 6. * F[2] - 20. * F[3] + 11. * F[4]) / (12. * h * h);
        cerr << "ERROR: one-sided outer stencil only implemented for first and second derivatives" << endl;
        exit(1);
    }
    return fivePointDeriv(st.h, order, F[0],F[1],F[2],F[3],F[4]);
}

double BSSNSlice::d_z(bssn_var var, int index, int order = 1 )
{
    DerivStencil st;
    build_stencil(index, st);
    return d_z_with(var, st, order);
}

//second z-derivative, basically syntactic sugar for d_z(var, index, 2)
double BSSNSlice::d_zz(bssn_var var, int index )
{
    return d_z(var, index, 2);
}

//converts the current slice to a tangherlini BH of mass m; must have set states vector size already. Also returns mass so it can be set in spacetime object
double BSSNSlice::make_tangherlini (double m, double min_chi, double D)
{
    int n_gridpoints = states2.size();

    // Use cell-centred grid spacing to match slice convention (z = (j+0.5)*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    double psi_power = -4. / (D - 3.); //12. / ((D - 3.) * (1. - D));

    for (int j = 0; j < n_gridpoints; j++)
    {
        double r = ( (j + 0.5) * dr);

        //conformal factor in isotropic convention
        double psi = 1. + pow(m, D - 3.) / (4 * pow(r, D - 3.));

        //conformal factor and conformally rescaled metric components are all powers of X in this gauge
    states2[j].bssn.chi = max(min_chi, pow(psi, psi_power));
    states2[j].bssn.h_zz = 1.;
    states2[j].bssn.h_ww = 1.;

    states2[j].bssn.A_zz = 0.;
    states2[j].bssn.A_ww = 0.;
    states2[j].bssn.K = 0.;
    states2[j].bssn.c_chris_Z = 0.;
    states2[j].csf.phi_re = 0.;
    states2[j].csf.phi_im = 0.;
    states2[j].csf.K_phi_re = 0.;
    states2[j].csf.K_phi_im = 0.;

    states2[j].bssn.alpha = pow(states2[j].bssn.chi, 0.5);
    states2[j].bssn.beta = 0.;

    }
    return m;
}

// Helper to apply high-order damping at a single gridpoint j
void Spacetime::apply_damping_at_point(int j,
                                       BSSNSlice* slice_ptr,
                                       std::vector<State>& rhs_states2,
                                       std::vector<double>& rhs_theta)
{
    // refinement factor local to j
    int res_fac = 1 << (refinement_levels[j] - 1);
    if (BSSNSlice::forced_res_fac > 0 && j <= BSSNSlice::forced_below_j) //Berger-Oliger buffer points
        res_fac = BSSNSlice::forced_res_fac;

    // indices for 7-point centered stencil with symmetry fix near the origin
    std::array<int, 7> J = {j - 3 * res_fac, j - 2 * res_fac, j - 1 * res_fac, j, j + 1 * res_fac, j + 2 * res_fac, j + 3 * res_fac};
    if (j < 3) {
        for (int k = 0; k < 7; ++k)
            if (J[k] < 0) J[k] = -J[k] - BSSNSlice::origin_mirror_shift;
    }

    const double h = dr * res_fac;
    const double d_mult = damping_factor * ((h * h * h * h * h) / 64.0);

    // references for State stencil
    State& sJ0 = slice_ptr->states2[J[0]];
    State& sJ1 = slice_ptr->states2[J[1]];
    State& sJ2 = slice_ptr->states2[J[2]];
    State& sJ3 = slice_ptr->states2[J[3]];
    State& sJ4 = slice_ptr->states2[J[4]];
    State& sJ5 = slice_ptr->states2[J[5]];
    State& sJ6 = slice_ptr->states2[J[6]];

    State damping_corr = d_mult * sevenPointDeriv(h, 6, sJ0, sJ1, sJ2, sJ3, sJ4, sJ5, sJ6);

    // the stencil above reflects every variable as even across z = 0; redo the odd variables (the contracted
    // Christoffel symbol and the shift) with the correct sign flip for reflected points
    if (j < 3)
    {
        std::array<double, 7> sgn;
        for (int k = 0; k < 7; ++k) sgn[k] = (j + (k - 3) * res_fac < 0) ? -1. : 1.;
        auto odd_corr = [&](auto get)
        {
            return d_mult * sevenPointDeriv<double>(h, 6, sgn[0] * get(sJ0), sgn[1] * get(sJ1), sgn[2] * get(sJ2),
                                                    get(sJ3), sgn[4] * get(sJ4), sgn[5] * get(sJ5), sgn[6] * get(sJ6));
        };
        damping_corr.bssn.c_chris_Z = odd_corr([](const State& st) { return st.bssn.c_chris_Z; });
        damping_corr.bssn.beta      = odd_corr([](const State& st) { return st.bssn.beta; });
    }

    // CCZ4: taper the dissipation of chi (a zero-speed field) to zero in the outer region
    if (use_CCZ4 && KO_chi_outer_frac > 0.)
        damping_corr.bssn.chi *= chi_ko_w[j];

    rhs_states2[j] = rhs_states2[j] + damping_corr;

    if (use_CCZ4 && rhs_theta.size() == static_cast<size_t>(n_gridpoints))
    {
        double& tJ0 = slice_ptr->theta[J[0]];
        double& tJ1 = slice_ptr->theta[J[1]];
        double& tJ2 = slice_ptr->theta[J[2]];
        double& tJ3 = slice_ptr->theta[J[3]];
        double& tJ4 = slice_ptr->theta[J[4]];
        double& tJ5 = slice_ptr->theta[J[5]];
        double& tJ6 = slice_ptr->theta[J[6]];
        double theta_corr = d_mult * sevenPointDeriv<double>(h, 6, tJ0, tJ1, tJ2, tJ3, tJ4, tJ5, tJ6);
        rhs_theta[j] += theta_corr;
    }
}

// Apply Sommerfeld-like outer boundary conditions to rhs
// Updates rhs.states2 at [last_active_j] and [last_active_j - res_fac]
void Spacetime::sommerfeld_BC(BSSNSlice& rhs, int res_fac_outer, BSSNSlice* slice_ptr)
{
    // Sommerfeld BC writes rhs at global indices last_active_j - res_fac and
    // last_active_j (the outermost 2 active points). Only the rank that owns
    // the outer edge should apply it.
    if (!owns_outer())
        return;

    int res_fac = res_fac_outer; //number of points skipped per slot in stencil

    // constraint-preserving BCs: keep the interior (one-sided) RHS of the BSSN variables computed by slice_rhs
    BSSNState cpbc_interior_rhs[2];
    double cpbc_interior_theta_rhs[2] = {0., 0.};
    if (constraint_preserving_BCs)
    {
        cpbc_interior_rhs[0] = rhs.states2[last_active_j - res_fac].bssn;
        cpbc_interior_rhs[1] = rhs.states2[last_active_j].bssn;
        if (use_CCZ4)
        {
            cpbc_interior_theta_rhs[0] = rhs.theta[last_active_j - res_fac];
            cpbc_interior_theta_rhs[1] = rhs.theta[last_active_j];
        }
    }

    // asymptotic states and their derivatives where relevant
    State asymp_state{}; // matter part zero by default
    // radii of the two outermost active points (including the cell-centring offset)
    const double r4 = (last_active_j - res_fac + grid_offset) * dr;
    const double r5 = (last_active_j + grid_offset) * dr;
    asymp_state.bssn = (BSSNState){chi_asymp(r4), h_zz_asymp(r4), h_ww_asymp(r4), 0., 0., 0., chris_Z_asymp(r4), alpha_asymp(r4), 0.};
    State asymp_deriv{}; // r-derivative of asymptotic expansion (matter zero)
    asymp_deriv.bssn = (BSSNState){d_chi_asymp(r4), d_h_zz_asymp(r4), d_h_ww_asymp(r4), 0., 0., 0., d_chris_Z_asymp(r4), d_alpha_asymp(r4), 0.};

    //version without spatially varying asymptotics
    if (!spatially_varying_BC)
    {
        asymp_state.bssn = (BSSNState){1., 1., 1., 0., 0., 0., 0., 1., 0.};
        asymp_deriv.bssn = (BSSNState){0., 0., 0., 0., 0., 0., 0., 0., 0.};
    }

    // Falloff powers p in d_t f = -v [d_r (f - f_inf) + p (f - f_inf) / r]: the BC is then exact for f = f_inf + c / r^p.
    // Metric variables: static (Coulomb-like) falloff p = D - 3, so a residual c / r^{D-3} mass mismatch is tolerated.
    // Matter fields: radiative falloff p = (D - 2) / 2. Both reduce to p = 1 in D = 4.
    const double p_metric = D - 3.;
    const double p_matter = 0.5 * (D - 2.);
    State N{};
    N.bssn = (BSSNState){p_metric, p_metric, p_metric, p_metric, p_metric, p_metric, p_metric, p_metric, p_metric};
    N.csf = (CSF){p_matter, p_matter, p_matter, p_matter};
    
    if (add_real_field)
        N.rsf = (RSF){p_matter, p_matter};

    //outermost 5 active gridpoints
    State& s1 = slice_ptr->states2[last_active_j - 4 * res_fac];
    State& s2 = slice_ptr->states2[last_active_j - 3 * res_fac];
    State& s3 = slice_ptr->states2[last_active_j - 2 * res_fac];
    State& s4 = slice_ptr->states2[last_active_j - 1 * res_fac];
    State& s5 = slice_ptr->states2[last_active_j];

    //limiting characteristic speeds at infinity
    // characteristic speeds
    State char_speeds{};
    double cs = sqrt(2. / s4.bssn.alpha);
    char_speeds.bssn = (BSSNState){cs, 1., 1., 1., 1., cs, 1., cs, 1.};
    char_speeds.csf = (CSF){1., 1., 1., 1.};

    if (add_real_field)
        char_speeds.rsf = (RSF){1., 1.};

    rhs.states2[last_active_j - res_fac] = (-1.) * char_speeds * ( p4_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                         + N * (s4 - asymp_state) / r4 - asymp_deriv);

    //update variable asymptotic states to outermost edge
    if (spatially_varying_BC)
    {
        asymp_state.bssn = (BSSNState){chi_asymp(r5), h_zz_asymp(r5), h_ww_asymp(r5), 0., 0., 0., chris_Z_asymp(r5), alpha_asymp(r5), 0.};
        asymp_deriv.bssn = (BSSNState){d_chi_asymp(r5), d_h_zz_asymp(r5), d_h_ww_asymp(r5), 0., 0., 0., d_chris_Z_asymp(r5), d_alpha_asymp(r5), 0.};
    }

    rhs.states2[last_active_j] = (-1) * char_speeds * ( p5_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                 + N * (s5 - asymp_state) / r5 - asymp_deriv);

    // Gamma-constraint-preserving condition on Gamma^z at the two outermost points (gamma_constraint_BC): replace the
    // Sommerfeld RHS by d/dt of the algebraic expression
    //   Gamma^z(h) = h_zz' / (2 h_zz^2) - (n/2) h_ww' / (h_zz h_ww) + n (1/h_ww - 1/h_zz) / z,
    // built from the (interior + Sommerfeld) RHS of h_zz, h_ww with one-sided z-derivatives. This keeps
    // G^z = Gamma^z - Gamma^z(h) fixed at the boundary instead of driving Gamma^z independently of h, which otherwise
    // reflects outgoing gauge pulses as inward-propagating G^z violations.
    if (!constraint_preserving_BCs && gamma_constraint_BC)
    {
        const double n = D - 2.;
        const double h = dr * res_fac;
        const int jj[5] = {last_active_j - 4 * res_fac, last_active_j - 3 * res_fac, last_active_j - 2 * res_fac,
                           last_active_j - res_fac, last_active_j};
        double a[5], b[5], a_t[5], b_t[5]; // h_zz, h_ww and their time derivatives on the 5 outermost points
        for (int k = 0; k < 5; k++)
        {
            a[k]   = slice_ptr->states2[jj[k]].bssn.h_zz;  b[k]   = slice_ptr->states2[jj[k]].bssn.h_ww;
            a_t[k] = rhs.states2[jj[k]].bssn.h_zz;         b_t[k] = rhs.states2[jj[k]].bssn.h_ww;
        }
        for (int q = 3; q < 5; q++)
        {
            const double z = (q == 3) ? r4 : r5;
            auto d1 = [&](const double* f) { return (q == 3) ? p4_stencil(h, f[0], f[1], f[2], f[3], f[4])
                                                             : p5_stencil(h, f[0], f[1], f[2], f[3], f[4]); };
            const double A = a[q], B = b[q], A_z = d1(a), B_z = d1(b), A_t = a_t[q], B_t = b_t[q], A_zt = d1(a_t), B_zt = d1(b_t);
            rhs.states2[jj[q]].bssn.c_chris_Z = 0.5 * A_zt / (A * A) - A_z * A_t / (A * A * A)
                                              - 0.5 * n * (B_zt / (A * B) - B_z * (A_t / (A * A * B) + B_t / (A * B * B)))
                                              + n * (A_t / (A * A) - B_t / (B * B)) / z;
        }
    }

    if (use_CCZ4)
    {
        double& t1 = slice_ptr->theta[last_active_j - 4 * res_fac];
        double& t2 = slice_ptr->theta[last_active_j - 3 * res_fac];
        double& t3 = slice_ptr->theta[last_active_j - 2 * res_fac];
        double& t4 = slice_ptr->theta[last_active_j - 1 * res_fac];
        double& t5 = slice_ptr->theta[last_active_j];

        rhs.theta[last_active_j - res_fac] = -(p4_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                             + p_metric * t4 / r4);

        rhs.theta[last_active_j] = -(p5_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                   + p_metric * t5 / r5);
    }

    // constraint-preserving BCs overwrite all BSSN / CCZ4 variables at the outer points (matter keeps Sommerfeld)
    if (constraint_preserving_BCs)
        constraint_preserving_BC(rhs, res_fac, slice_ptr, cpbc_interior_rhs, cpbc_interior_theta_rhs);
}

//enforce tracelessness of A
void Spacetime::make_A_traceless(BSSNSlice* slice_ptr)
{
    double n = D - 2.;

    // MPI: only clamp owned points; ghost values are refreshed by exchange_ghosts.
    for (int j = owned_begin(); j < owned_end(); j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

        const double h_zz = slice_ptr->states2[j].bssn.h_zz;
        const double h_ww = slice_ptr->states2[j].bssn.h_ww;
        const double A_zz = slice_ptr->states2[j].bssn.A_zz;
        const double A_ww = slice_ptr->states2[j].bssn.A_ww;

            double A = A_zz / h_zz + n * A_ww / h_ww;
        slice_ptr->states2[j].bssn.A_zz = A_zz - h_zz * A / (D - 1.);
        slice_ptr->states2[j].bssn.A_ww = A_ww - h_ww * A / (D - 1.);
    }

}


//reads data from a boson_star object we have already solved for into initial BSSN slice.
void BSSNSlice::read_BS_data (BosonStar& boson_star, int BS_resolution_factor, bool isotropic, bool cell_centered)
{
    //prepare to fill states array

    if ((boson_star.n_gridpoints + BS_resolution_factor - 1) % BS_resolution_factor != 0 )
        cout << "WARNING: Incompatible resolution factor used in read_BS_data" << endl;

    int n_gridpoints = (boson_star.n_gridpoints + BS_resolution_factor - 1) / BS_resolution_factor;
    states2.resize(n_gridpoints); // initialize composite State vector

    if (use_CCZ4)
       theta.resize(n_gridpoints);

    R = boson_star.R;
    // Grid spacing; when cell_centered==true, cell centers are at z=(j+0.5)*dr
    double dr = R / ((double)n_gridpoints - 1.0);

    double D = boson_star.D;

    //
    // Precompute BosonStar spacing and sizes for interpolation
    int N_bs = boson_star.n_gridpoints;
    double dr_bs = boson_star.R / (double)(max(1, N_bs - 1));

    // Build x-grid for spline (uniform nodes j*dr_bs)
    std::vector<double> X_nodes(N_bs);
    for (int k = 0; k < N_bs; ++k) X_nodes[k] = k * dr_bs;

    // Build contiguous arrays and their splines
    tk::spline s_X, s_A, s_phi, s_psi_iso, s_A_iso, s_phi_iso, s_pert, s_pert_iso;
    tk::spline s_rsf_iso, s_rsf;
    if (isotropic) {
        if (!boson_star.psi_iso_array.empty()) {
            s_psi_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_psi_iso.set_points(X_nodes, boson_star.psi_iso_array, tk::spline::cspline);
        }
        if (!boson_star.A_iso_array.empty())   {
            s_A_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_A_iso.set_points(X_nodes, boson_star.A_iso_array, tk::spline::cspline);
        }
        if (!boson_star.phi_iso_array.empty()) {
            s_phi_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_phi_iso.set_points(X_nodes, boson_star.phi_iso_array, tk::spline::cspline);
        }
        if (boson_star.add_real_field && !boson_star.rsf_iso_array.empty()) {
            s_rsf_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_rsf_iso.set_points(X_nodes, boson_star.rsf_iso_array, tk::spline::cspline);
        }
        if (boson_star.perturb && !boson_star.pert_iso_array.empty()) {
            s_pert_iso.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_pert_iso.set_points(X_nodes, boson_star.pert_iso_array, tk::spline::cspline);
        }
    } else {
        std::vector<double> bs_X(N_bs), bs_A(N_bs), bs_phi(N_bs);
        for (int k = 0; k < N_bs; ++k) {
            bs_X[k] = boson_star.state[k].X;
            bs_A[k] = boson_star.state[k].A;
            bs_phi[k] = boson_star.state[k].phi;
        }
        s_X.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_A.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_phi.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        s_X.set_points(X_nodes, bs_X, tk::spline::cspline);
        s_A.set_points(X_nodes, bs_A, tk::spline::cspline);
        s_phi.set_points(X_nodes, bs_phi, tk::spline::cspline);
        if (boson_star.perturb && !boson_star.pert_array.empty()) {
            s_pert.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_pert.set_points(X_nodes, boson_star.pert_array, tk::spline::cspline);
        }
        if (boson_star.add_real_field && !boson_star.rsf_array.empty()) {
            std::vector<double> bs_rsf(N_bs);
            for (int k = 0; k < N_bs; ++k) bs_rsf[k] = boson_star.rsf_array[k];
            s_rsf.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
            s_rsf.set_points(X_nodes, bs_rsf, tk::spline::cspline);
        }
    }

    for (int j = 0; j < n_gridpoints; j++)
    {
        // cell- or vertex-centered coordinate
        double z = (j + 0.5 * (double)cell_centered) * dr;

        // Interpolate X (areal conformal factor) or psi (isotropic conformal factor)
        double X_interp = 1.;
        double psi_interp = 1.;
        if (isotropic)
        {
            // interpolate isotropic psi (conformal factor)
            psi_interp = s_psi_iso(z);
            states2[j].bssn.chi = pow(psi_interp, -4.);
            states2[j].bssn.h_zz = 1.; //isotropic gauge is conformally flat
            states2[j].bssn.h_ww = 1.;
        }
        else
        {
            X_interp = s_X(z);
            states2[j].bssn.chi = pow(X_interp, -2. / (D - 1.));
            states2[j].bssn.h_zz = pow(X_interp, (2. * D - 4.) / (D - 1.));
            states2[j].bssn.h_ww = pow(X_interp, -2. / (D - 1.));
        }

        //for time-independent static, spherically symmetric boson stars K_ij = 0-- may need to do more work when we consider perturbations.
        states2[j].bssn.A_zz = 0.;
        states2[j].bssn.A_ww = 0.;
        states2[j].bssn.K = 0.;

        // start at t = 0 so the scalar field is real regardless of phase
        double A_interp = 0.;
        double phi_interp = 0.;
        if (isotropic)
        {
            A_interp = s_A_iso(z);
            phi_interp = s_phi_iso(z);
            states2[j].csf.phi_re = A_interp;
            states2[j].csf.phi_im = 0.;
            states2[j].bssn.alpha = exp(phi_interp);
            states2[j].bssn.beta = 0.;
        }
        else
        {
            A_interp = s_A(z);
            phi_interp = s_phi(z);
            states2[j].csf.phi_re = A_interp;
            states2[j].csf.phi_im = 0.;
            states2[j].bssn.alpha = exp(phi_interp);
            states2[j].bssn.beta = 0.;
        }

        // Real scalar field initial data (psi). If boson_star requests a real field, interpolate
        // from its isotropic or areal arrays into states2.rsf.psi. Leave K_psi = 0.
        if (boson_star.add_real_field) {
            if (isotropic) {
                if (!boson_star.rsf_iso_array.empty())
                    states2[j].rsf.psi = s_rsf_iso(z);
                else
                    states2[j].rsf.psi = 0.0;
            } else {
                if (!boson_star.rsf_array.empty())
                    states2[j].rsf.psi = s_rsf(z);
                else
                    states2[j].rsf.psi = 0.0;
            }
        } else {
            states2[j].rsf.psi = 0.0;
        }
        states2[j].rsf.K_psi = 0.0;

        //starting BS real means its momentum is imaginary (with 0 starting shift)
    states2[j].csf.K_phi_re = 0.;

        double pert_correction = 0.;
        if (boson_star.perturb)
        {
            double pert_interp = 0.;
            if (isotropic)
                pert_interp = s_pert_iso(z);
            else
                pert_interp = s_pert(z);

            pert_correction = 1. * (boson_star.omega / states2[j].bssn.alpha) * pert_interp; //conserves Noether charge at 1st order in perturbation
        }

    states2[j].csf.K_phi_im = - boson_star.omega * states2[j].csf.phi_re / (2. * states2[j].bssn.alpha) + pert_correction;

        if (use_CCZ4)
            theta[j] = 0.;
    }

    //separate loop to fill in values for the contracted Christoffel symbols as they require derivatives of h (hence info off j value)
    //(the origin mirror used by d_z must match the grid: set it here since initialize() only sets it after reading the data)
    origin_mirror_shift = cell_centered ? 1 : 0;
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j + 0.5 * (double)cell_centered) * dr; // cell- or vertex-centred coordinate, as in the loop above

        // inverse metric components, for convenience
        double h_ZZ = 1 / states2[j].bssn.h_zz;
        double h_WW = 1 / states2[j].bssn.h_ww;

        states2[j].bssn.c_chris_Z = 0.5 * h_ZZ * h_ZZ * d_z(v_h_zz,j) + (D - 2.) * (-0.5 * h_ZZ * h_WW * d_z(v_h_ww,j));

        // extra geometric term, which goes to zero at z = 0
        if (z > 1e-8)
            states2[j].bssn.c_chris_Z += (D - 2.) * (h_WW * (1 - h_ZZ * states2[j].bssn.h_ww) / z);

        // hard enforce 0 for isotropic data
        if (isotropic)
            states2[j].bssn.c_chris_Z = 0.;
    }
}

//writes BSSN evolution variables from a particular slice to text file
void BSSNSlice::write_slice(std::string file_name)
{
    int length = states2.size();
    double dr = R / ((double)length - 1.0);
    std::ofstream data_file{file_name};

     if (!data_file)
    {
        cerr << "SliceData.dat could not be opened for writing!\n";
        exit(1);
    }

    for (int j = 0; j < length; j++)
    {
        //if (!active_points[j]) //perform no computations on inactive points
            //continue;

        double theta0 = 0.; //if using CCZ4, write theta at the end
        if (use_CCZ4)
            theta0 = theta[j];

    data_file <<  std::setprecision (16) << (j) * dr << "   " << states2[j].bssn.chi << "    " << states2[j].bssn.h_zz << "    " << states2[j].bssn.h_ww  << "    " << states2[j].bssn.A_zz
        << "   " << states2[j].bssn.A_ww << "    " << states2[j].bssn.K << "    " << states2[j].bssn.c_chris_Z  << "    " << states2[j].csf.phi_re << "    " << states2[j].csf.phi_im
        << "   " << states2[j].csf.K_phi_re << "    " << states2[j].csf.K_phi_im << "    " << states2[j].bssn.alpha << "    " << states2[j].bssn.beta << "    " <<  theta0
        << "   " << states2[j].rsf.psi << "    " << states2[j].rsf.K_psi << endl;
    }
}

//reads slice data from checkpoint file of form checkpoint(time).dat w/o parentheses
void BSSNSlice::read_checkpoint(int time, int n_gridpoints)
{
    states2.resize(n_gridpoints);

    std::ifstream checkpoint_file("checkpoint" + std::to_string(time) + ".dat");
    if (!checkpoint_file.is_open())
    {
        cerr << "Could not open checkpoint file!" << endl;
        exit(1);
    }

    string line;
    int j = 0;
    while (std::getline(checkpoint_file, line))
    {
        std::istringstream iss(line);
    double z;
    // read into temporaries to populate both legacy states and new states2
    double chi,h_zz,h_ww,A_zz,A_ww,K,c_chris_Z,phi_re,phi_im,K_phi_re,K_phi_im,alpha,beta;

        if (use_CCZ4)
        {
            double& t = theta[j];

            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta >> t)
            {
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                j++;
            }
        }
        else
        {
            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta)
            {
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                j++;
            }
        }
    }
}

//determines which refinement level j belongs to. 1 is highest resolution; for each level up the resolution halves.
int BSSNSlice::get_refinement_level(int j, std::vector<int>& refinement_points)
{
    //cout << "Size is " <<  refinement_points.size() << endl;

    if (refinement_points.empty()) //always return 1 if no refinement
        return 1;

    int n_refinements = refinement_points.size();
    int level = 1;
    int k = 0;

    while ( k < n_refinements && j >= refinement_points[k] - (1 << k))//check this difference -- meant to ensure we can use a stencil at points spaced by 2^(k + 1) safely at j
    {
        level++;
        k++;
    }

    return level;
}

//convergence test for three slices. Resolutions must differ by factor of 2. Also, currently does not support refinement.
void slice_convergence_test (BSSNSlice& sl, BSSNSlice& sm, BSSNSlice& sh)
{
    if (sl.R != sm.R || sm.R != sh.R)
    {
        cout << "ERROR: Convergence test requested on slices of different radii";
        exit(1);
    }

    if (2*sl.states2.size() - 1 != sm.states2.size()  || 2*sm.states2.size() - 1 != sh.states2.size()   )
    {
        cout << "ERROR: Convergence test requested on slices with incompatible sizes";
        exit(1);
    }

    //write to file
    std::ofstream conv_file{ "conv.dat" };

    // If we couldn't open the output file stream for writing
    if (!conv_file)
    {
        // Print an error and exit
        cerr << "Error: conv.dat could not be opened for writing\n";
        exit(1);
    }

    //write change in A between med/low and high/med resolution to file.
    for (unsigned int j = 0; j < sl.states2.size(); j++)
    {
    // For 4th-order convergence the refinement factor scaling is 2^4 = 16
    conv_file << j * sl.R /(sl.states2.size() - 1) << "   " << sm.states2[2*j].csf.phi_re -  sl.states2[j].csf.phi_re  << "    " << 16.*(sh.states2[4*j].csf.phi_re  -  sm.states2[2*j].csf.phi_re ) << endl;
    }
}


//these let us call d_z and d_zz on the current slice without explicitly referencing it
//must remember to update current_slice_ptr appropriately!!!
double Spacetime::d_z(bssn_var var, int index, int order = 1)
{
    const BSSNSlice* sl = current_slice_ptr;
    const int n = (int)sl->states2.size();
    if (sl != deriv_cache_slice || index != deriv_cache_index || n != deriv_cache_n || sl->R != deriv_cache_R
        || BSSNSlice::forced_res_fac != deriv_cache_frf || BSSNSlice::forced_below_j != deriv_cache_fbj)
    {
        //new point: build the stencil once and gather the five values of every variable (the RHS takes ~24 derivatives
        //per point, so each call below is then a short dot product with no indirection)
        sl->build_stencil(index, deriv_cache);
        deriv_cache_slice = sl; deriv_cache_index = index; deriv_cache_n = n; deriv_cache_R = sl->R;
        deriv_cache_frf = BSSNSlice::forced_res_fac; deriv_cache_fbj = BSSNSlice::forced_below_j;

        const bool has_theta = (sl->theta.size() == sl->states2.size());
        for (int k = 0; k < 5; k++)
        {
            const State& st = sl->states2[deriv_cache.J[k]];
            deriv_vals[v_chi][k] = st.bssn.chi;           deriv_vals[v_h_zz][k] = st.bssn.h_zz;
            deriv_vals[v_h_ww][k] = st.bssn.h_ww;         deriv_vals[v_A_zz][k] = st.bssn.A_zz;
            deriv_vals[v_A_ww][k] = st.bssn.A_ww;         deriv_vals[v_K][k] = st.bssn.K;
            deriv_vals[v_c_chris_Z][k] = st.bssn.c_chris_Z;
            deriv_vals[v_phi_re][k] = st.csf.phi_re;      deriv_vals[v_phi_im][k] = st.csf.phi_im;
            deriv_vals[v_K_phi_re][k] = st.csf.K_phi_re;  deriv_vals[v_K_phi_im][k] = st.csf.K_phi_im;
            deriv_vals[v_alpha][k] = st.bssn.alpha;       deriv_vals[v_beta][k] = st.bssn.beta;
            deriv_vals[v_theta][k] = has_theta ? sl->theta[deriv_cache.J[k]] : 0.;
            deriv_vals[v_psi][k] = st.rsf.psi;            deriv_vals[v_K_psi][k] = st.rsf.K_psi;
        }
        //odd parity of the contracted Christoffel symbol and the shift across z = 0
        if (deriv_cache.flip0) { deriv_vals[v_c_chris_Z][0] = -deriv_vals[v_c_chris_Z][0]; deriv_vals[v_beta][0] = -deriv_vals[v_beta][0]; }
        if (deriv_cache.flip1) { deriv_vals[v_c_chris_Z][1] = -deriv_vals[v_c_chris_Z][1]; deriv_vals[v_beta][1] = -deriv_vals[v_beta][1]; }
    }

    const double* F = deriv_vals[var];
    const BSSNSlice::DerivStencil& st = deriv_cache;
    if (st.outer_pos != 0)
    {
        const double h = st.h;
        if (order == 1)
            return (st.outer_pos == 5) ? p5_stencil(h, F[0], F[1], F[2], F[3], F[4]) : p4_stencil(h, F[0], F[1], F[2], F[3], F[4]);
        if (order == 2) // one-sided 5-point second derivatives
            return (st.outer_pos == 5) ? (11. * F[0] - 56. * F[1] + 114. * F[2] - 104. * F[3] + 35. * F[4]) / (12. * h * h)
                                       : (-F[0] + 4. * F[1] + 6. * F[2] - 20. * F[3] + 11. * F[4]) / (12. * h * h);
        cerr << "ERROR: one-sided outer stencil only implemented for first and second derivatives" << endl;
        exit(1);
    }
    //same expressions as fivePointDeriv (h2 = pow(h, 2) precomputed per point)
    if (order == 1)
        return ((1.0 / 12.0) * F[0] - (2.0 / 3.0) * F[1] + (2.0 / 3.0) * F[3] - (1.0 / 12.0) * F[4]) / st.h;
    if (order == 2)
        return ( (-1.0 / 12.0) * F[0] + (4.0 / 3.0) * F[1] - (5.0 / 2.0) * F[2] + (4.0 / 3.0) * F[3] - (1.0 / 12.0) * F[4]) / st.h2;
    return fivePointDeriv(st.h, order, F[0], F[1], F[2], F[3], F[4]);
}

double Spacetime::d_zz(bssn_var var, int index)
{
    return d_z(var, index, 2);
}

double Spacetime::V(const double A)
{
    if (!solitonic)
        return mu * mu * A * A + lambda * pow(A,4);

    else
        return mu * mu * A * A * pow((1. - 2. * pow(A / sigma, 2)), 2) + lambda * pow(A,4) ;
}

double Spacetime::dV(const double A)
{
    if (!solitonic)
        return mu * mu + 2 * A * A * lambda;

    else
        return mu * mu - 8. * mu * mu * pow(A / sigma, 2) + 12. * mu * mu * pow(A / sigma, 4) + 2. * A * A * lambda;
}

// Constraint-preserving outer boundary conditions for BSSN in spherical symmetry, following Alcubierre & Torres,
// arXiv:1407.8529, translated to SBSE variables and generalised to D dimensions (originally D = 4): a = h_zz, b = h_ww, e^{4 phi} = 1/chi, A_a = A_zz / h_zz,
// Delta^r = Gamma^z. All BSSN variables are evolved to the outer boundary with their own equations using one-sided
// stencils (slice_rhs + BSSNSlice::one_sided_outer); at the outermost point the two incoming modes are then replaced:
//  * gauge mode  w^alpha_- = W K - d_z ln(alpha), W = sqrt(h_zz f / chi): Sommerfeld-type condition
//    w^alpha_- = (alpha - 1) / (alpha z), imposed through its time derivative, which fixes d_t K;
//  * light-speed mode w^l_- = sqrt(h_zz / chi) (A_a - 2K/3) - (2/3)(...): the incoming constraint mode
//    Omega_- = h_zz H / (6 chi) + d_z C_Delta / 2 + sqrt(h_zz / chi) M  (Omega_- ~ d_z w^l_-) is set to zero via the
//    paper's time-derivative method, d_t w^l_- = (d_t w^l_-)_interior - v_l Omega_-, v_l = alpha sqrt(chi / h_zz),
//    which fixes d_t A_a (A_ww follows from tracelessness, A_b = -A_a / 2).
// Here H is the BSSN Hamiltonian constraint, M = d_z A_a - (2/3) d_z K - (3/2) A_a d_z chi / chi + (3/2) A_a (2/z + d_z h_ww / h_ww)
// - 8 pi j_z the radial momentum constraint, and C_Delta = h_zz (Gamma^z - Gamma^z(h)). Matter fields keep Sommerfeld BCs.
// The derivation assumes negligible shift at the boundary.
// CCZ4 extension (D = 4, zero shift, principal part; c^2 = chi / h_zz, Sigma = h_zz Gamma^z - d_a - 2 d_b - 8 d_z phi):
//   d_t Theta ~ (alpha c^2 / 2) d_z Sigma,  d_t Sigma ~ 2 alpha d_z Theta  =>  Theta -/+ (c/2) Sigma are light-speed
//   outgoing / incoming fields; the gauge and A_a light-speed fields are the BSSN ones with K -> Khat = K - 2 Theta and
//   Delta -> Gamma-hat. Three incoming fields are fixed: Theta by the Sommerfeld condition d_t Theta = -v_l (d_z Theta + Theta/z);
//   the gauge field as above with Khat; and w^l_- through Omega_- + (4/3) sqrt(h_zz / chi) d_z Theta = 0, with H built from
//   the Z-corrected Ricci tensor (which d_z w^l_- matches at principal order).
void Spacetime::constraint_preserving_BC(BSSNSlice& rhs, int res_fac, BSSNSlice* slice_ptr, const BSSNState interior_rhs[2],
                                         const double interior_theta_rhs[2])
{
    const double n = D - 2.;
    const int jN = last_active_j;
    const double h = dr * res_fac;
    const int jj[5] = {jN - 4 * res_fac, jN - 3 * res_fac, jN - 2 * res_fac, jN - res_fac, jN};

    // all BSSN (and CCZ4 Theta) variables at both outer points: interior equations with one-sided stencils
    rhs.states2[jN - res_fac].bssn = interior_rhs[0];
    rhs.states2[jN].bssn = interior_rhs[1];
    if (use_CCZ4)
    {
        rhs.theta[jN - res_fac] = interior_theta_rhs[0];
        rhs.theta[jN] = interior_theta_rhs[1];
    }

    const BSSNState& s = slice_ptr->states2[jN].bssn;
    BSSNState& r = rhs.states2[jN].bssn;
    const double z = (jN + grid_offset) * dr;
    const double chi = std::max(min_chi, s.chi);
    const double a = s.h_zz, b = s.h_ww, alpha = s.alpha, K = s.K;
    const double A_a = s.A_zz / a;

    // ---- CCZ4 Theta mode (light speed): Sommerfeld condition d_t Theta = -v_l (d_z Theta + Theta / z)
    const double v_l = alpha * sqrt(chi / a);
    const double Theta = use_CCZ4 ? slice_ptr->theta[jN] : 0.;
    double dz_Theta = 0., dTheta_t = 0.;
    if (use_CCZ4)
    {
        double th[5];
        for (int k = 0; k < 5; k++)
            th[k] = slice_ptr->theta[jj[k]];
        dz_Theta = p5_stencil(h, th[0], th[1], th[2], th[3], th[4]);
        // keep Theta's (lower-order) kappa1 damping in the boundary condition: without it, the condition is inconsistent with
        // the interior equation at the boundary and slowly drives a global Theta / K mode, faster for larger kappa1
        const double k1_b = covariant_Z4 ? c1 : c1 * alpha;
        const double p_Theta = 0.5 * (D - 2.); // radiative falloff for the Theta (constraint) mode; = 1 in D = 4
        const double Theta_t_new = theta_off ? 0. : -v_l * (dz_Theta + p_Theta * Theta / z) - 0.5 * k1_b * (D + c2 * (D - 2.)) * Theta;
        dTheta_t = Theta_t_new - rhs.theta[jN];
        rhs.theta[jN] = Theta_t_new;
    }

    // ---- gauge mode: d_t Khat from d_t [W Khat - d_z ln alpha] = d_t [(alpha - 1) / (alpha z)], Khat = K - 2 Theta
    // (in BSSN Theta = 0 and Khat = K)
    const double Khat = K - 2. * Theta;
    const double f  = shock_gauge ? (1. + shock_fac / (alpha * alpha)) : (one_log_fac / alpha);
    const double fp = shock_gauge ? (-2. * shock_fac / (alpha * alpha * alpha)) : (-one_log_fac / (alpha * alpha));
    const double W = sqrt(a * f / chi);
    double lt[5]; // d_t ln(alpha) on the outer stencil
    for (int k = 0; k < 5; k++)
        lt[k] = rhs.states2[jj[k]].bssn.alpha / slice_ptr->states2[jj[k]].bssn.alpha;
    const double d_t_q = p5_stencil(h, lt[0], lt[1], lt[2], lt[3], lt[4]); // d_t d_z ln(alpha)
    const double Wt_over_W = 0.5 * (r.h_zz / a + fp * r.alpha / f - r.chi / chi);
    const double p_alpha = D - 3.; // static falloff alpha - 1 ~ 1/z^{D-3} in the Sommerfeld-type gauge condition
    const double Khat_t_new = (p_alpha * r.alpha / (alpha * alpha * z) + d_t_q - Khat * W * Wt_over_W) / W;
    const double K_t_new = Khat_t_new + 2. * (use_CCZ4 ? rhs.theta[jN] : 0.);
    const double dK_t = K_t_new - r.K;
    const double dKhat_t = dK_t - 2. * dTheta_t;
    r.K = K_t_new;

    // ---- light-speed constraint mode Omega_- = 0
    double Aa[5], Kv[5], chiv[5], hwwv[5], CD[5];
    for (int k = 0; k < 5; k++)
    {
        const BSSNState& sk = slice_ptr->states2[jj[k]].bssn;
        const double zk = (jj[k] + grid_offset) * dr;
        Aa[k] = sk.A_zz / sk.h_zz; Kv[k] = sk.K; chiv[k] = sk.chi; hwwv[k] = sk.h_ww;
        const double dhzz = slice_ptr->d_z(v_h_zz, jj[k], 1), dhww = slice_ptr->d_z(v_h_ww, jj[k], 1);
        const double Gam_h = 0.5 * dhzz / (sk.h_zz * sk.h_zz) - 0.5 * n * dhww / (sk.h_zz * sk.h_ww)
                           + n * (1. / sk.h_ww - 1. / sk.h_zz) / zk;
        CD[k] = sk.h_zz * (sk.c_chris_Z - Gam_h);
    }
    auto d5 = [&](const double* v) { return p5_stencil(h, v[0], v[1], v[2], v[3], v[4]); };

    const double Ham_N = chi * h_ZZ[jN] * R_zz[jN] + n * chi * h_WW[jN] * R_ww[jN]
                       - (A_a * A_a + n * pow(s.A_ww / b, 2)) + n * K * K / (D - 1.) - 16. * M_PI * rho[jN];
    // radial momentum constraint in D dimensions: d_z A_a - (n/d) d_z K + (n+1) A_a (1/z + d_z ln h_ww / 2 - d_z ln chi / 2)
    const double dsp = D - 1.; // number of spatial dimensions
    const double Mom_N = d5(Aa) - (n / dsp) * d5(Kv) + (n + 1.) * A_a * (1. / z + 0.5 * d5(hwwv) / b - 0.5 * d5(chiv) / chi)
                       - 8. * M_PI * j_z[jN];
    // in CCZ4: R_zz, R_ww include the Z terms (Ham_N is the hatted Hamiltonian), Gamma^z is Gamma-hat, and the
    // light-speed field is w^l_- = sqrt(h_zz / chi) (A_a - 2 Khat / 3) - (2/3)(...), whose radial derivative matches
    // Omega_- + (4/3) sqrt(h_zz / chi) d_z Theta at principal order
    // Omega_- = (n-1) h_zz H / (2 d chi) + d_z C_Delta / 2 + sqrt(h_zz/chi) (M + (2n/d) d_z Theta), whose principal part is
    // d_z of the incoming light-speed field w^l_- = sqrt(h_zz/chi) (A_a - (n/d) Khat) + ... (characteristic analysis in D dims;
    // reduces to h_zz H / (6 chi) + ... + (4/3) d_z Theta in D = 4)
    const double Omega_m = (n - 1.) * a * Ham_N / (2. * dsp * chi) + 0.5 * d5(CD) + sqrt(a / chi) * (Mom_N + (2. * n / dsp) * dz_Theta);

    const double dA_a_t = (n / dsp) * dKhat_t - alpha * chi / a * Omega_m;
    r.A_zz += a * dA_a_t;
    r.A_ww += -(b / n) * dA_a_t; // tracelessness: A_b = -A_a / n
}

//Asymptotic (Schwarzschild-Tangherlini) values of chi and alpha and their radial derivatives, used by the Sommerfeld BCs.
//Areal coordinates:     f = 1 - mu / r^{D-3}, gamma_rr = 1/f  =>  chi = f^{1/(D-1)},  alpha = sqrt(f)
//Isotropic coordinates: psi = 1 + q, q = mu / (4 r^{D-3})  =>  chi = psi^{-4/(D-3)}, alpha = (1 - q) / (1 + q)
//with mu = 16 pi M_bc / ((D-2) Omega_{D-2}), i.e. mu = 2 M_bc in D = 4, where M_bc is the Misner-Sharp mass at the
//extraction radius (slice_mass is wrong for non-conformally-flat slices). make_tangherlini data is isotropic.
double Spacetime::tangherlini_mu() const
{
    const double Omega = 2. * pow(M_PI, 0.5 * (D - 1.)) / tgamma(0.5 * (D - 1.)); // area of unit (D-2)-sphere
    return 16. * M_PI * M_bc / ((D - 2.) * Omega);
}

double Spacetime::chi_asymp(double r)
{
    const double mu_T = tangherlini_mu();
    if (isotropic || make_tangherlini)
        return pow(1. + 0.25 * mu_T / pow(r, D - 3.), -4. / (D - 3.));

    return pow(1. - mu_T / pow(r, D - 3.), 1. / (D - 1.));
}

double Spacetime::alpha_asymp(double r)
{
    const double mu_T = tangherlini_mu();
    if (isotropic || make_tangherlini)
    {
        const double q = 0.25 * mu_T / pow(r, D - 3.);
        return (1. - q) / (1. + q);
    }

    return sqrt(1. - mu_T / pow(r, D - 3.));
}

double Spacetime::d_chi_asymp(double r)
{
    const double mu_T = tangherlini_mu();
    if (isotropic || make_tangherlini)
    {
        const double q = 0.25 * mu_T / pow(r, D - 3.);
        return 4. * q * pow(1. + q, -4. / (D - 3.) - 1.) / r;
    }

    const double f = 1. - mu_T / pow(r, D - 3.);
    return (D - 3.) * mu_T * pow(r, 2. - D) * pow(f, 1. / (D - 1.) - 1.) / (D - 1.);
}

double Spacetime::d_alpha_asymp(double r)
{
    const double mu_T = tangherlini_mu();
    if (isotropic || make_tangherlini)
    {
        const double q = 0.25 * mu_T / pow(r, D - 3.);
        return 2. * (D - 3.) * q / (r * (1. + q) * (1. + q));
    }

    const double f = 1. - mu_T / pow(r, D - 3.);
    return 0.5 * (D - 3.) * mu_T * pow(r, 2. - D) / sqrt(f);
}

//Asymptotic conformal metric: identity in isotropic coordinates; in areal coordinates gamma_zz = 1/f, gamma_ww = 1 give
//h_zz = f^{-(D-2)/(D-1)}, h_ww = f^{1/(D-1)}. chris_Z_asymp uses the same expression as the initial-data construction.
double Spacetime::h_zz_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 1.;
    return pow(1. - tangherlini_mu() / pow(r, D - 3.), -(D - 2.) / (D - 1.));
}

double Spacetime::h_ww_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 1.;
    return pow(1. - tangherlini_mu() / pow(r, D - 3.), 1. / (D - 1.));
}

double Spacetime::d_h_zz_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 0.;
    const double mu_T = tangherlini_mu();
    const double f = 1. - mu_T / pow(r, D - 3.);
    const double df = (D - 3.) * mu_T * pow(r, 2. - D);
    return -(D - 2.) / (D - 1.) * pow(f, -(D - 2.) / (D - 1.) - 1.) * df;
}

double Spacetime::d_h_ww_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 0.;
    const double mu_T = tangherlini_mu();
    const double f = 1. - mu_T / pow(r, D - 3.);
    const double df = (D - 3.) * mu_T * pow(r, 2. - D);
    return pow(f, 1. / (D - 1.) - 1.) * df / (D - 1.);
}

double Spacetime::chris_Z_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 0.;
    const double h_ZZ = 1. / h_zz_asymp(r);
    const double h_WW = 1. / h_ww_asymp(r);
    return 0.5 * h_ZZ * h_ZZ * d_h_zz_asymp(r) - 0.5 * (D - 2.) * h_ZZ * h_WW * d_h_ww_asymp(r)
           + (D - 2.) * h_WW * (1. - h_ZZ * h_ww_asymp(r)) / r;
}

double Spacetime::d_chris_Z_asymp(double r)
{
    if (isotropic || make_tangherlini)
        return 0.;
    const double eps = 1e-4 * r; // smooth analytic function: centred difference is ample
    return (chris_Z_asymp(r + eps) - chris_Z_asymp(r - eps)) / (2. * eps);
}

 // returns RHS as {chi, eta} where eta is the radial derivative of chi
 vector<double> Spacetime::ham_init_rhs(double r, double chi, double eta)
 {
    int k = floor(r / dr) - 1;
    while (k * dr < r)
        k++;

    const BSSNSlice& s = slices[0]; //initial slice
    current_slice_ptr = &slices[0];

    int j0 = bound(k - 1, 0, n_gridpoints - 4);

    //cubic interpolation to get state values off gridpoints (for midstep solving)
    double phi_re = cubic_interp(r, s.states2[j0].csf.phi_re, s.states2[j0 + 1].csf.phi_re, s.states2[j0 + 2].csf.phi_re, s.states2[j0 + 3].csf.phi_re, j0, dr);
    double phi_im = cubic_interp(r, s.states2[j0].csf.phi_im, s.states2[j0 + 1].csf.phi_im, s.states2[j0 + 2].csf.phi_im,s.states2[j0 + 3].csf.phi_im, j0, dr);
    double K_phi_re = cubic_interp(r, s.states2[j0].csf.K_phi_re, s.states2[j0 + 1].csf.K_phi_re, s.states2[j0 + 2].csf.K_phi_re, s.states2[j0 + 3].csf.K_phi_re,j0, dr);
    double K_phi_im = cubic_interp(r, s.states2[j0].csf.K_phi_im, s.states2[j0 + 1].csf.K_phi_im, s.states2[j0 + 2].csf.K_phi_im, s.states2[j0 + 3].csf.K_phi_im,j0, dr);

    //WARNING: these are very similarly named to Spacetime member variables introduced to avoid re-computing derivatives during evolution; do not confuse!
    double dz_phi_re = cubic_interp(r, d_z(v_phi_re, j0), d_z(v_phi_re, j0 + 1), d_z(v_phi_re, j0 + 2), d_z(v_phi_re, j0 + 3),j0,  dr);
    double dz_phi_im = cubic_interp(r, d_z(v_phi_im, j0), d_z(v_phi_im, j0 + 1), d_z(v_phi_im, j0 + 2), d_z(v_phi_im, j0 + 3),j0,  dr);

    double mod_phi = sqrt(phi_re * phi_re + phi_im + phi_im);
    double rho0 = 2. * (K_phi_im * K_phi_im + K_phi_re * K_phi_re) +  0.5 * chi * (pow(dz_phi_re,2) + pow(dz_phi_im,2)) + 0.5 * V(mod_phi);
    
    //account for real scalar field if needed
    if (add_real_field)
    {   
        //double psi = cubic_interp(r, s.states2[j0].rsf.psi, s.states2[j0 + 1].rsf.psi, s.states2[j0 + 2].rsf.psi, s.states2[j0 + 3].rsf.psi,j0, dr);
        double dz_psi= cubic_interp(r, d_z(v_psi, j0), d_z(v_psi, j0 + 1), d_z(v_psi, j0 + 2), d_z(v_psi, j0 + 3),j0,  dr);
        double K_psi = cubic_interp(r, s.states2[j0].rsf.K_psi, s.states2[j0 + 1].rsf.K_psi, s.states2[j0 + 2].rsf.K_psi, s.states2[j0 + 3].rsf.K_psi,j0, dr);

        rho0 += 0.5 * chi * dz_psi * dz_psi + 2. * K_psi * K_psi;
    }

    double chi_rhs = eta;
    double eta_rhs = 5. * eta * eta / (4. * chi) + 8. * M_PI * rho0;

    if (r > 0.) eta_rhs -= 2 * eta / r;
    else eta_rhs /= 3.; // in the limit r -> 0, eta / r is replaced by d_r(eta), so net effect is to divide the usual RHS by 3!

    vector<double> rhs = {chi_rhs, eta_rhs};
    return rhs;

 }

 void Spacetime::solve_initial_ham(bool quiet)
 {
    if (!isotropic)
    {
        cout << "WARNING: Spacetime Hamiltonian solver called in non-isotropic coordinates; aborting request" << endl;
        return;
    }

    BSSNSlice& s = slices[0]; //initial slice
    double c1, c2, c3, c4, e1, e2, e3, e4;
    double eta = 0;

    if (!quiet)
        cout << "\n Starting initial Spacetime Hamiltonian constraint solver..." << endl;

   for (int j = 0; j < n_gridpoints - 1; j++) //for (int j = n_gridpoints - 1; j > 0; j--)
    {
        double r = j * dr;

        //RK4 solver for Ham. constraint
        //pretty inefficient but OK for now
        c1 = ham_init_rhs(r, s.states2[j].bssn.chi, eta)[0];
        e1 = ham_init_rhs(r, s.states2[j].bssn.chi, eta)[1];

        c2 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[0];
        e2 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c1, eta + 0.5 * dr * e1)[1];

        c3 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[0];
        e3 = ham_init_rhs(r + dr / 2., s.states2[j].bssn.chi + 0.5 * dr * c2, eta + 0.5 * dr * e2)[1];

        c4 = ham_init_rhs(r + dr, s.states2[j].bssn.chi + dr * c3, eta + dr * e3)[0];
        e4 = ham_init_rhs(r + dr, s.states2[j].bssn.chi + dr * c3, eta + dr * e3)[1];

        s.states2[j + 1].bssn.chi = s.states2[j].bssn.chi + (dr / 6.) * (c1 + 2 * c2 + 2 * c3 + c4);
        eta = eta + (dr / 6.) * (e1 + 2 * e2 + 2 * e3 + e4);

        if (isnan(s.states2[j + 1].bssn.chi))
        {
            cout << "ERROR: constraint solver produced nan's on step " << j  << endl;
            return; //exit(1); //should probably make a param about whether this fails
        }
    }

    //initializes alpha = sqrt(chi)-- mostly for formation runs where static lapse solver may have failed.
    for (int j = 0; j < n_gridpoints; j++)
        s.states2[j].bssn.alpha = sqrt(s.states2[j].bssn.chi);


    if (!quiet)
        cout << "Successfully ran constraint solver" << endl;

 }

//compute auxiliary quantities at jth point
void Spacetime::auxiliary_quantities_at_point(BSSNSlice* slice_ptr, int j)
{
    const double n = D - 2.;
    const double z = (j + grid_offset) * dr;

    const double invDm1 = 1.0 / (D - 1.);

    // shorthand versions of the state on desired slice for convenience
    const double chi = max(min_chi, slice_ptr->states2[j].bssn.chi);
    const double h_zz = slice_ptr->states2[j].bssn.h_zz;
    const double h_ww = slice_ptr->states2[j].bssn.h_ww;
    const double c_chris_Z = slice_ptr->states2[j].bssn.c_chris_Z;
    const double phi_re = slice_ptr->states2[j].csf.phi_re;
    const double phi_im = slice_ptr->states2[j].csf.phi_im;
    const double beta = slice_ptr->states2[j].bssn.beta;

    // inverse metric components
    h_ZZ[j] = 1.0 / h_zz;
    h_WW[j] = 1.0 / h_ww;
    const double inv_hzz = h_ZZ[j];
    const double inv_hww = h_WW[j];
    const double inv_chi = 1.0 / chi;

    //Christoffel symbols that are only needed for other aux variables
    double chris_zzz = 0.5 * d_z_h_zz;
    double chris_wwz = 0.5 * d_z_h_ww;
    double chris_zww = -0.5 * d_z_h_ww;
    if (z > min_z)
        chris_zww += (h_zz - h_ww) / z; // extra term that vanishes as z -> 0
    double chris_Wwz = inv_hww * chris_wwz;

    // christoffel symbols stored for evolution equations
    chris_Zzz[j] = inv_hzz * chris_zzz;
    chris_Zww[j] = inv_hzz * chris_zww;
    const double cZzz = chris_Zzz[j];
    const double cZww = chris_Zww[j];

    // auxiliary constraint
    if (z > min_z)
        aux_constraint[j] = (d_z_beta + n * beta / z) * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);
    else
        aux_constraint[j] = (n + 1.) * d_z_beta * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);

    // 3-covariant derivatives of alpha and tracefree parts
    D_zz_alpha[j] = d_zz(v_alpha, j) - cZzz * d_z_alpha + 0.5 * d_z_chi * d_z_alpha * inv_chi;
    
    const double d_alpha_z = (z <= min_z) ? d_zz(v_alpha, j) : (d_z_alpha / z); // regularized alpha_z/z at the origin
    D_ww_alpha[j] = d_alpha_z - (cZww + 0.5 * h_ww * inv_hzz * d_z_chi * inv_chi) * d_z_alpha;

    D_zz_alpha_TF[j] = n * (D_zz_alpha[j] - h_zz * inv_hww * D_ww_alpha[j]) * invDm1;
    D_ww_alpha_TF[j] = (D_ww_alpha[j] - h_ww * inv_hzz * D_zz_alpha[j]) * invDm1;

    // cached derivatives
    const double dzz_chi = d_zz(v_chi, j);
    const double cD_zz_chi = dzz_chi - cZzz * d_z_chi;
    double d_chi_z = (z <= min_z) ? dzz_chi : (d_z_chi / z); // d_z(chi) / z replaced with 2nd deriv at 0
    double c_chris_Z_overz = (z <= min_z) ? d_z(v_c_chris_Z, j) : (c_chris_Z / z); // c_chris_Z / z replaced with 2nd deriv at 0
    const double dzz_h_zz = d_zz(v_h_zz, j);
    const double dzz_h_ww = d_zz(v_h_ww, j);
    double d_h_ww_z = (z <= min_z) ? dzz_h_ww : (d_z_h_ww / z); // d_z(h_ww) / z replaced with 2nd deriv at 0
    const double dzchi_over_chi = d_z_chi * inv_chi;

    // Ricci chi-part
    double R_zz_chi = n * (cD_zz_chi + (0.5 * inv_hww * d_z_h_ww * d_z_chi + d_chi_z) - d_z_chi * dzchi_over_chi) * (0.5 * inv_chi);
    double R_ww_chi = h_ww * inv_hzz * (cD_zz_chi + (2. * D - 5.) * (0.5 * inv_hww * d_z_h_ww * d_z_chi + d_chi_z) - (D - 1.) * d_z_chi * (0.5 * dzchi_over_chi)) * (0.5 * inv_chi);

    // Ricci conformal-part
    double R_zz_c_t1 = (z <= min_z) ? (-0.5 * dzz_h_ww) : ((h_zz - h_ww) / (z * z) - 0.5 * d_z_h_zz / z);
    double h_zw_diff = (z <= min_z) ? (0.5 * (dzz_h_zz - dzz_h_ww)) : ((h_zz - h_ww) / (z * z));

    const double dz_Gamma = d_z(v_c_chris_Z, j); // derivative of the evolved (BSSN Gamma-tilde or CCZ4 Gamma-hat) variable
    double R_zz_c = n * inv_hww * (R_zz_c_t1 + chris_Wwz * chris_wwz + 2. * chris_Wwz * chris_zww)
                    - 0.5 * inv_hzz * dzz_h_zz + h_zz * dz_Gamma
                    + c_chris_Z * chris_zzz + 3. * inv_hzz * chris_zzz * cZzz;

    double R_ww_c = -0.5 * inv_hzz * dzz_h_ww + 3. * inv_hzz * chris_Wwz * chris_wwz - 0.5 * n * inv_hww * d_h_ww_z
                    + inv_hww * (cZww * chris_zww + 2. * cZww * chris_wwz)
                    + h_ww * c_chris_Z_overz + c_chris_Z * chris_wwz - inv_hww * h_zw_diff;

    R_zz[j] = R_zz_c + R_zz_chi;
    R_ww[j] = R_ww_c + R_ww_chi;

    // Difference between the Ricci terms built from the evolved Gamma variable (used in the evolution equations and in
    // Ham[j]) and those built from the metric alone, Gamma^z(h) = h^zz Gamma^z_zz + n h^ww Gamma^z_ww. Used to report the
    // pure ADM Hamiltonian constraint (constraints_from_metric): in BSSN the evolved Gamma-tilde differs from Gamma(h) by
    // the Gamma constraint, in CCZ4 by 2 Z^z / chi, and both leak into Ham[j] otherwise (in CCZ4, Ham[j] is then the
    // "hatted" Hamiltonian H + 2 div Z that the Theta equation damps, which overstates the constraint improvement).
    {
        const double Gamma_h = inv_hzz * cZzz + n * inv_hww * cZww;
        const double G_diff = c_chris_Z - Gamma_h; // = 2 Z^z / chi in CCZ4
        // d_z Gamma^z(h) analytically: Gamma(h) = (u + n v / h_ww) / h_zz, u = d h_zz / (2 h_zz), v = chris_zww
        const double du = 0.5 * dzz_h_zz * inv_hzz - 0.5 * d_z_h_zz * d_z_h_zz * inv_hzz * inv_hzz;
        const double dv = (z <= min_z) ? (-0.5 * dzz_h_ww + 0.5 * (dzz_h_zz - dzz_h_ww))
                                       : (-0.5 * dzz_h_ww + (d_z_h_zz - d_z_h_ww) / z - (h_zz - h_ww) / (z * z));
        const double dGamma_h = -d_z_h_zz * inv_hzz * Gamma_h + (du + n * (dv * inv_hww - chris_zww * d_z_h_ww * inv_hww * inv_hww)) * inv_hzz;
        const double dG = dz_Gamma - dGamma_h;
        const double G_overz = (z <= min_z) ? dG : (G_diff / z);
        const double dR_zz = h_zz * dG + G_diff * chris_zzz;
        const double dR_ww = h_ww * G_overz + G_diff * chris_wwz;
        ham_corr[j] = chi * inv_hzz * dR_zz + n * chi * inv_hww * dR_ww;
    }

    if (use_CCZ4) //if CCZ4 is used, add supplementary terms to the (extended) Ricci tensor and fill out theta_Z
    {
        double theta = 0.5 * chi * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);
        theta_Z[j] = theta;
        const double dR_zz_Z = theta * (chi * d_z_h_zz + h_zz * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_zzz * inv_chi;
        const double dR_ww_Z = theta * (chi * d_z_h_ww - h_ww * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_wwz * inv_chi;
        R_zz[j] += dR_zz_Z;
        R_ww[j] += dR_ww_Z;
        ham_corr[j] += chi * inv_hzz * dR_zz_Z + n * chi * inv_hww * dR_ww_Z; // explicit 2 D_(i Z_j) terms
    }

    // traceless part of Ricci components
    R_zz_TF[j] = n * (R_zz[j] - h_zz * inv_hww * R_ww[j]) * invDm1;
    R_ww_TF[j] = (R_ww[j] - h_ww * inv_hzz * R_zz[j]) * invDm1;

    // |phi|
    double mod_phi = sqrt(phi_re * phi_re + phi_im * phi_im);
    if (j == 0) test_ctr = mod_phi;

    // matter quantities via ComplexScalarField helpers (reuse member to avoid per-point construction)
    CSF csf = slice_ptr->states2[j].csf;
    const BSSNState& bssn = slice_ptr->states2[j].bssn;
    // Only phi derivatives are used; K_phi derivatives and second derivatives not used in stress-energy
    CSF d_z_csf{d_z_phi_re, d_z_phi_im, 0.0, 0.0};
    CSF d_zz_csf{0.0, 0.0, 0.0, 0.0};

    rho[j] = csf_model.rho(csf, bssn, d_z_csf, d_zz_csf);
    j_z[j] = csf_model.j_z(csf, bssn, d_z_csf, d_zz_csf);
    S_zz[j] = csf_model.S_zz(csf, bssn, d_z_csf, d_zz_csf);
    S_ww[j] = csf_model.S_ww(csf, bssn, d_z_csf, d_zz_csf);

    // If configured, add massless real scalar field contributions
    if (add_real_field)
    {
        RSF rsf = slice_ptr->states2[j].rsf;
        RSF d_z_rsf{d_z(v_psi, j), d_z(v_K_psi, j)};
        RSF d_zz_rsf{0.0, 0.0};

        rho[j] += rsf_model.rho(rsf, bssn, d_z_rsf, d_zz_rsf);
        j_z[j] += rsf_model.j_z(rsf, bssn, d_z_rsf, d_zz_rsf);
        S_zz[j] += rsf_model.S_zz(rsf, bssn, d_z_rsf, d_zz_rsf);
        S_ww[j] += rsf_model.S_ww(rsf, bssn, d_z_rsf, d_zz_rsf);
    }

    // Trace S = gamma^{ij} S_ij = chi * (h^{zz} S_zz + n h^{ww} S_ww)
    S[j] = chi * (inv_hzz * S_zz[j] + n * inv_hww * S_ww[j]);

    S_zz_TF[j] = S_zz[j] - S[j] * h_zz / (chi * (D - 1.));
    S_ww_TF[j] = S_ww[j] - S[j] * h_ww / (chi * (D - 1.));
}

//compute auxiliary quantities on a slice. Should no longer be called in slice_rhs (pointwise evaluator called directly for efficiency to avoid re-computing derivatives)
void Spacetime::compute_auxiliary_quantities(BSSNSlice* slice_ptr, bool derivatives_computed)
{
    //double n = D - 2.;
    dr = R / ((double)n_gridpoints - 1.0);

    // MPI: only compute aux at owned points; each rank's aux[] arrays are
    // therefore only valid in [owned_begin, owned_end). Downstream users
    // (compute_diagnostics, slice_mass, ...) already restrict their reads.
    for (int j = owned_begin(); j < owned_end(); j++)
    {
        d_z_chi = d_z(v_chi,j);
        d_z_h_zz = d_z(v_h_zz,j);
        d_z_h_ww = d_z(v_h_ww,j);
        d_z_phi_re = d_z(v_phi_re,j);
        d_z_phi_im = d_z(v_phi_im,j);
        d_z_alpha = d_z(v_alpha,j);
        d_z_beta = d_z(v_beta,j);

        auxiliary_quantities_at_point(slice_ptr, j);
    }
}

//returns a slice corresponding to the RHS of the BSSN evolution equations. Also computes needed auxiliary quantities first
BSSNSlice Spacetime::slice_rhs(BSSNSlice* slice_ptr, int skipped_levels, const std::vector<int>* point_list, bool do_outer_BC)
{
    BSSNSlice rhs;
    slice_rhs_into(slice_ptr, rhs, skipped_levels, point_list, do_outer_BC);
    return rhs;
}

void Spacetime::slice_rhs_into(BSSNSlice* slice_ptr, BSSNSlice& rhs, int skipped_levels, const std::vector<int>* point_list, bool do_outer_BC)
{
    double n = D - 2.;

    //define return slice and size its states array appropriately
    rhs.R = R;
    rhs.has_BH = slice_ptr->has_BH;
    rhs.refinement_points = slice_ptr->refinement_points;

    rhs.use_CCZ4 = use_CCZ4;
    if ((int)rhs.states2.size() != n_gridpoints) rhs.states2.resize(n_gridpoints); // composite State vector (zero-initialised when first allocated)

    int max_ref = (refinement_levels.size() == 0 ) ? 1 : refinement_levels[refinement_levels.size() - 1]; //refinement level at outermost bdry
    const int res_fac_outer = (max_ref > 0) ? (1 << (max_ref - 1)) : 1;

    // MPI: only compute RHS at points this rank owns. Under replicated storage,
    // rhs.states2 is still sized n_gridpoints; non-owned slots are left zero
    // (initialized above) and are never read by anyone downstream because the
    // corresponding rank owns them. The upper bound n_gridpoints - 2 stays
    // because the last two active gridpoints are set by sommerfeld_BC on the
    // rank that owns the outer edge.
    const int j_lo = std::max(owned_begin(), 0);
    // With constraint-preserving BCs the interior equations are also evaluated (with one-sided stencils) at the two
    // outermost points; constraint_preserving_BC then modifies only the incoming modes there.
    const int j_hi = std::min(owned_end(), constraint_preserving_BCs ? n_gridpoints : n_gridpoints - 2);
    const int n_eval = point_list ? (int)point_list->size() : (j_hi - j_lo);
    for (int jj = 0; jj < n_eval; jj++)
    {
        const int j = point_list ? (*point_list)[jj] : (j_lo + jj);
        if (!point_list && (!active_points[j] || refinement_levels[j] > max_ref - skipped_levels))//perform no calculations and skip on inactive points
            continue;

        double z = (j + grid_offset) * dr;

        //shorthand versions of the state vars on desired slice for convenience
        const double chi = max(min_chi, slice_ptr->states2[j].bssn.chi);
        const double h_zz = slice_ptr->states2[j].bssn.h_zz;
        const double h_ww = slice_ptr->states2[j].bssn.h_ww;
        const double A_zz = slice_ptr->states2[j].bssn.A_zz;
        const double A_ww = slice_ptr->states2[j].bssn.A_ww;
        const double K = slice_ptr->states2[j].bssn.K;
        const double c_chris_Z = slice_ptr->states2[j].bssn.c_chris_Z;

        const double alpha = slice_ptr->states2[j].bssn.alpha;
        const double beta = slice_ptr->states2[j].bssn.beta;


        if (wave_mode) //converts to wave eq'n solver with beta the time derivative of phi_re for testing purposes
        {
            rhs.states2[j].bssn = (BSSNState){0.,0.,0.,0.,0.,0.,0.,0.,0.};
            rhs.states2[j].csf.phi_re = beta;
            rhs.states2[j].bssn.beta = d_zz(v_phi_re, j);
            continue;
        }

        //store local variables for commonly-used derivatives to avoid unnecessary re-computation
        d_z_chi = d_z(v_chi,j);
        d_z_h_zz = d_z(v_h_zz,j);
        d_z_h_ww = d_z(v_h_ww,j);
        d_z_phi_re =  d_z(v_phi_re,j);
        d_z_phi_im =  d_z(v_phi_im,j);
        d_z_alpha =  d_z(v_alpha,j);
        d_z_beta = d_z(v_beta,j);
        d_z_K = d_z(v_K,j);

        auxiliary_quantities_at_point(slice_ptr, j);

        double beta_z = ((z <= min_z) ? d_z_beta : (beta / z)); //beta / z replaced by z-derivative at z = 0

        // precompute common inverses/constants for this point
        const double invDm1 = 1.0 / (D - 1.);
        const double inv_hzz = h_ZZ[j];
        const double inv_hww = h_WW[j];
        const double cZzz = chris_Zzz[j];
        const double cZww = chris_Zww[j];

        //fill out RHS of the BSSN evolution system
        rhs.states2[j].bssn.chi  = beta * d_z_chi  + 2.* chi * (alpha * K - d_z_beta - n * beta_z ) * invDm1;
        rhs.states2[j].bssn.h_zz = beta * d_z_h_zz + 2. * h_zz * d_z_beta - 2. * h_zz * (d_z_beta + n * beta_z) * invDm1 - 2. * alpha * A_zz;
        rhs.states2[j].bssn.h_ww = beta * d_z_h_ww - 2. * h_ww * (d_z_beta - beta_z) * invDm1 - 2. * alpha * A_ww;

        if (!use_CCZ4) //different evolution eq'n used in CCZ4
        {
            rhs.states2[j].bssn.K = beta * d_z_K - chi * inv_hzz * D_zz_alpha[j] + alpha * inv_hzz * inv_hzz * A_zz * A_zz + alpha * K * K * invDm1
                            + n * inv_hww * (alpha * A_ww * A_ww / h_ww - chi * D_ww_alpha[j]) + 8. * M_PI * alpha * (S[j] + (D - 3.) * rho[j]) / n;
        }

        rhs.states2[j].bssn.A_zz = beta * d_z(v_A_zz, j) + 2. * A_zz * d_z_beta - 2. * A_zz * (d_z_beta + n * beta_z) * invDm1 + alpha * K * A_zz
                                - 2. * alpha * A_zz * inv_hzz * A_zz + chi * (alpha * R_zz_TF[j] - D_zz_alpha_TF[j] - 8. * M_PI * alpha * S_zz_TF[j]);
        rhs.states2[j].bssn.A_ww = beta * d_z(v_A_ww, j) - 2. * A_ww * (d_z_beta -  beta_z) * invDm1 + alpha * A_ww * (K - 2 * inv_hww * A_ww)
                                + chi * (alpha * R_ww_TF[j] - D_ww_alpha_TF[j] - 8. * M_PI * alpha * S_ww_TF[j]);

        rhs.states2[j].bssn.c_chris_Z = beta * d_z(v_c_chris_Z,j) + 2. * c_chris_Z * (d_z_beta + n * beta_z) * invDm1 + inv_hzz * d_zz(v_beta,j) - c_chris_Z * d_z_beta
                                    + (D - 3.) * inv_hzz * d_zz(v_beta,j) * invDm1 - 2. * (D - 2.) * alpha * inv_hzz * d_z_K * invDm1
                                    - A_zz * inv_hzz * inv_hzz * ( (D - 1.) * alpha * d_z_chi / chi + 2 * d_z_alpha)
                                    + 2 * alpha * (cZzz * inv_hzz * inv_hzz * A_zz + n * cZww * inv_hww * inv_hww * A_ww)
                                    - sigma_BSSN * aux_constraint[j] - 16. * M_PI * alpha * j_z[j]* inv_hzz;

        if (z >= min_z) //add terms in d_z(beta) / z - beta /z^2 when z =/= 0
            rhs.states2[j].bssn.c_chris_Z += n * (inv_hww + (D - 3.) * inv_hzz * invDm1) * (d_z_beta / z - beta / (z * z));
    
        //Gauge variable update using moving puncture evolution, unless evolve_shift is off in which case do not update beta
        double f = (shock_gauge) ? (1. + shock_fac / (alpha * alpha) ): (one_log_fac/ alpha); //bona-masso f
        rhs.states2[j].bssn.alpha = beta * d_z_alpha - f * (alpha * alpha) * K;
        double gamma_eff = gamma_fac;
        if (constraint_preserving_BCs && shift_taper_frac > 0.)
            gamma_eff *= shift_taper_w[j];
        rhs.states2[j].bssn.beta =  (evolve_shift)? (beta * d_z_beta + gamma_eff * c_chris_Z - eta * beta): 0.;
       
        // ComplexScalarField RHS for matter variables (uses internal origin regularization and geometric couplings)
        CSF d_z_csf{d_z_phi_re, d_z_phi_im, d_z(v_K_phi_re, j), d_z(v_K_phi_im, j)};
        // Optimization: K_phi second derivatives not used by RHS; pass zeros for those entries
        CSF d_zz_csf{d_zz(v_phi_re, j), d_zz(v_phi_im, j), 0.0, 0.0};

        const CSF csf_rhs = csf_model.rhs_CSF(
            slice_ptr->states2[j].csf,
            slice_ptr->states2[j].bssn,
            d_z_csf,
            d_zz_csf,
            d_z_alpha,
            d_z_chi,
            chris_Zzz[j],
            chris_Zww[j],
            z);

        rhs.states2[j].csf = csf_rhs;

        // RealScalarField RHS for massless real scalar (optional)
        if (add_real_field)
        {
            RSF d_z_rsf{ d_z(v_psi, j), d_z(v_K_psi, j) };
            // Second derivative needed only for psi
            RSF d_zz_rsf{ d_zz(v_psi, j), 0.0 };

            const RSF rsf_rhs = rsf_model.rhs_RSF(
                slice_ptr->states2[j].rsf,
                slice_ptr->states2[j].bssn,
                d_z_rsf,
                d_zz_rsf,
                d_z_alpha,
                d_z_chi,
                chris_Zzz[j],
                chris_Zww[j],
                z);

            rhs.states2[j].rsf = rsf_rhs;
        }
        else
        {
            rhs.states2[j].rsf = RSF{0.0, 0.0};
        }

        //add supplementary terms if CCZ4 evolution is on
        if (use_CCZ4)
        {
            const double theta = slice_ptr->theta[j];
            if ((int)rhs.theta.size() != n_gridpoints) rhs.theta.resize(n_gridpoints);

            // kappa1 damping: covariant_Z4 drops the lapse factor (kappa1 alpha -> kappa1) so the damping does not
            // switch off near punctures where alpha -> 0
            const double k1 = covariant_Z4 ? c1 : c1 * alpha;
            // optionally scale the kappa1 damping of the Theta equation inside an apparent horizon: the energy non-conservation
            // of that term acting on the puncture's (truncation-error) Theta is the main driver of the horizon-mass drift in
            // under-resolved puncture evolutions; switching it off inside ~0.6 R_AH reduces the drift by ~30% (lower-order
            // change only, so hyperbolicity is untouched). Scaling the K- or Gamma-equation damping instead does not help.
            double k1_T = k1;
            if (k1_theta_inner_factor != 1. && ah_radius > 0.)
                k1_T *= 1. - (1. - k1_theta_inner_factor) * 0.5 * (1. - tanh((z - inner_damp_frac * ah_radius) / (0.1 * ah_radius)));

            rhs.states2[j].bssn.K = beta * d_z_K - chi * inv_hzz * D_zz_alpha[j]  - n * inv_hww * chi * D_ww_alpha[j]
                            + alpha * (chi * inv_hzz * R_zz[j] + n * chi * inv_hww * R_ww[j] + K * (K - 2. * theta)
                            + 8. * M_PI * (S[j] - (D - 1.) * rho[j]) / (D - 2.)) - (D - 1.) * k1 * (1. + c2) * theta;

            rhs.states2[j].bssn.A_zz += -2. * alpha * theta * A_zz;
            rhs.states2[j].bssn.A_ww += -2. * alpha * theta * A_ww;

            //note this is now \hat{\Gamma}^z, not \tilde{\Gamma}^z. So we first remove BSSN specific damping then add new terms.
            //The BSSN shift terms above use \hat{\Gamma}^z = \tilde{\Gamma}^z(h) + 2 Z^z / chi, which already contains the kappa3 = 1
            //term 2 [(2/(D-1)) Z^z d_k beta^k - Z^k d_k beta^z] / chi; the last line here corrects it to kappa3 = c3.
            rhs.states2[j].bssn.c_chris_Z +=   sigma_BSSN * aux_constraint[j] - 2. * d_z_alpha * theta * inv_hzz + 2. * alpha * inv_hzz * d_z(v_theta, j)
                                        - 2. * k1 * theta_Z[j] / chi - 4. * alpha * K * theta_Z[j] * invDm1 / chi
                                        + 2. * (c3 - 1.) * invDm1 * theta_Z[j] * ((3. - D) * d_z_beta + 2. * n * beta_z) / chi;

            rhs.states2[j].bssn.alpha +=  2. * alpha * alpha * f * theta; //4. * alpha * theta; //corrects for K -> K - 2 * theta in 1+log slicing
            //rhs.states2[j].bssn.beta +=  (evolve_shift)? ( -2. * gamma_fac * theta_Z[j] / chi): 0.; //replace theta hat -> theta tilde in CCZ4 moving puncture gauge; seems to produce instability...

            if (!theta_off)
                rhs.theta[j] = beta * d_z(v_theta, j) + 0.5 * alpha * (chi * inv_hzz * R_zz[j] - inv_hzz * inv_hzz * A_zz * A_zz + n * K * K * invDm1
                            + n * (chi * inv_hww * R_ww[j] - A_ww * A_ww / (h_ww * h_ww))  - 2. * K * theta - 2. * theta_Z[j] * d_z_alpha / alpha
                            - 16. * M_PI * rho[j]) - 0.5 * k1_T * (D + c2 * n) * theta;
            else
                rhs.theta[j] = 0.;
        }

        // add damping away from the outer edge. The centred 7-point KO stencil must not reach the outermost points, which are
        // set by the boundary conditions and not themselves dissipated: with the cut-off at 3 points (stencil touching the
        // last point) the truncated operator is not dissipative and drives a slow global Theta / K mode in CCZ4.
        if (damping_factor != 0. && j <= last_active_j - std::max(3, KO_outer_cutoff) * res_fac_outer)
            apply_damping_at_point(j, slice_ptr, rhs.states2, rhs.theta);
    }

    //apply Sommerfeld boundary conditions at outer boundary
    if (do_outer_BC)
        sommerfeld_BC(rhs, res_fac_outer, slice_ptr);
}

//TODO: make these work with MR
double Spacetime::slice_mass(BSSNSlice* slice_ptr)
{
    const double n = D - 2;
    double mass = 0;

    // Owned-range partial sum; global result via bss_mpi::sum.
    const int j_hi = std::min(owned_end(), (int)(0.95 * n_gridpoints) - 1);
    for (int j = owned_begin(); j < j_hi; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;
    const double& chi = slice_ptr->states2[j].bssn.chi;

        // M = (Omega_n / 16 pi) int psi^{(D+1)/(D-3)} R^(3) z^{D-2} dz, with psi^{(D+1)/(D-3)} = chi^{-(D+1)/4} (chi^{-5/4} in D = 4)
        // and R^(3) replaced via the Hamiltonian constraint. Exact only for conformally flat slices (h_ij = delta_ij) that satisfy
        // the constraints; otherwise it misses -(Omega_n / 16 pi) int psi R[h] z^{D-2} dz plus a Gamma^z surface term, and
        // picks up constraint violation inside BHs. Prefer slice_ms_mass for BH spacetimes.
        // Prefactor Omega_n / 16 pi = pi^{(D-3)/2} / (8 Gamma((D-1)/2)), = 1/4 in D = 4.
        mass += -0.125 * pow(M_PI, 0.5 * (D - 3.))  * quad_weights[j] * (Ham[j]- h_ZZ[j] * chi * R_zz[j] - n * h_WW[j] * chi * R_ww[j]) * pow(z, D - 2.) * pow(chi, -0.25 * (D + 1.)) / tgamma(0.5 * (D - 1.));
    }
    return bss_mpi::sum(mass);
}

// Misner-Sharp mass at gridpoint j (local; caller must own j and have valid ghosts for the d_z stencil):
//   m = (n Omega_n / 16 pi) R_ar^{n-1} (1 - nabla_a R_ar nabla^a R_ar),   n = D - 2,
//   R_ar = z sqrt(h_ww / chi),
//   nabla_a R_ar nabla^a R_ar = (chi / h_zz) (d_z R_ar)^2 - (R_ar K^theta_theta)^2,  K^theta_theta = A_ww / h_ww + K / (D - 1).
// Equals R_ar/2 (1 - ...) in D = 4, and the ADM mass of Schwarzschild-Tangherlini in any slicing. Returns 0 at z = 0.
double Spacetime::ms_mass_at_point(BSSNSlice* slice_ptr, int j)
{
    const double n = D - 2.;
    const double dr_local = R / ((double)n_gridpoints - 1.0);
    const double z = (j + grid_offset) * dr_local;
    if (z <= 0.)
        return 0.;

    const BSSNState& b = slice_ptr->states2[j].bssn;
    const double R_ar = z * sqrt(b.h_ww / b.chi);
    const double dz_R_ar = R_ar * (1. / z + 0.5 * slice_ptr->d_z(v_h_ww, j, 1) / b.h_ww - 0.5 * slice_ptr->d_z(v_chi, j, 1) / b.chi);
    const double K_th = b.A_ww / b.h_ww + b.K / (D - 1.);
    const double grad_R_sq = (b.chi / b.h_zz) * dz_R_ar * dz_R_ar - R_ar * R_ar * K_th * K_th;

    const double Omega_n = 2. * pow(M_PI, 0.5 * (D - 1.)) / tgamma(0.5 * (D - 1.)); // area of unit n-sphere
    return n * Omega_n / (16. * M_PI) * pow(R_ar, n - 1.) * (1. - grad_R_sq);
}

// Misner-Sharp mass at the active gridpoint nearest z = r_frac * R. MPI: collective; only the rank owning the
// extraction point contributes before the sum-reduce.
double Spacetime::slice_ms_mass(BSSNSlice* slice_ptr, double r_frac)
{
    const double dr_local = R / ((double)n_gridpoints - 1.0);

    // nearest gridpoint to the target radius, kept clear of the outer-boundary stencils, then stepped down to an active point
    int j = (int)std::lround(r_frac * R / dr_local - grid_offset);
    j = std::max(1, std::min(j, last_active_j - 6));
    while (j > 1 && !active_points[j])
        j--;

    const double contribution = (j >= global_lo && j < global_hi) ? ms_mass_at_point(slice_ptr, j) : 0.;
    return bss_mpi::sum(contribution);
}

double Spacetime::slice_charge(BSSNSlice* slice_ptr)
{
    const double n = D - 2;
    double charge = 0;

    const int j_hi = std::min(owned_end(), (int)(0.95 * n_gridpoints) - 1);
    for (int j = owned_begin(); j < j_hi; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;

    const double& chi = slice_ptr->states2[j].bssn.chi;
    const double& phi_re = slice_ptr->states2[j].csf.phi_re;
    const double& phi_im = slice_ptr->states2[j].csf.phi_im;
    const double& K_phi_re = slice_ptr->states2[j].csf.K_phi_re;
    const double& K_phi_im = slice_ptr->states2[j].csf.K_phi_im;

        if (D == 4.)
            charge += -8. * M_PI  * quad_weights[j] * (phi_re * K_phi_im - phi_im * K_phi_re) * z * z * pow(chi, -1.5);
        else
             charge += -4. * pow(M_PI, 0.5 * (D - 1.)) * quad_weights[j] * (phi_re * K_phi_im - phi_im * K_phi_re) * pow(z, n) * pow(chi, 0.5 * (1. - D)) / tgamma(0.5 * (D - 1.));
    }
    return bss_mpi::sum(charge);
}

// Compute energies E_phi (complex scalar) and E_psi (real scalar) by integrating
// the corresponding rho contributions over the slice using the spherical volume measure.
void Spacetime::compute_scalar_energies(BSSNSlice* slice_ptr)
{
    if (slice_ptr == nullptr || slice_ptr->states2.empty())
    {
        E_phi = 0.0;
        E_psi = 0.0;
        return;
    }

    // Local grid characteristics from the provided slice to avoid reliance on globals
    const int N = static_cast<int>(slice_ptr->states2.size());
    const double dr_local = slice_ptr->R / (static_cast<double>(N) - 1.0);
    const double offset = 0.5 * static_cast<double>(cell_centered); // match grid_offset convention

    // Surface area of the unit (D-2)-sphere
    const double sphere_area = (D == 4.)
        ? (4.0 * M_PI)
        : (2.0 * std::pow(M_PI, 0.5 * (D - 1.)) / std::tgamma(0.5 * (D - 1.)));

    double e_phi = 0.0;
    double e_psi = 0.0;

    // Integrate over inner 95% of the grid (to mirror other diagnostics that avoid boundary artifacts).
    // MPI: partial sum over owned range, then bss_mpi::sum below.
    const int Jmax = static_cast<int>(0.95 * N) - 1;
    const int j_lo = std::max(owned_begin(), 0);
    const int j_hi = std::min(owned_end(), Jmax);
    for (int j = j_lo; j < j_hi; ++j)
    {
        if (!active_points[j])
            continue;

        const double z = (j + offset) * dr_local;
        const BSSNState& bssn = slice_ptr->states2[j].bssn;
        const double chi = std::max(min_chi, bssn.chi);

        // Complex scalar field contribution
        {
            const CSF csf = slice_ptr->states2[j].csf;
            const CSF d_z_csf{slice_ptr->d_z(v_phi_re, j), slice_ptr->d_z(v_phi_im, j), 0.0, 0.0};
            const CSF d_zz_csf{0.0, 0.0, 0.0, 0.0};
            const double rho_phi = csf_model.rho(csf, bssn, d_z_csf, d_zz_csf);
            e_phi += quad_weights[j] * sphere_area * rho_phi * std::pow(z, D - 2.) * std::pow(chi, 0.5 * (1.0 - D));
        }

        // Real scalar field (optional)
        if (add_real_field)
        {
            const RSF rsf = slice_ptr->states2[j].rsf;
            const RSF d_z_rsf{slice_ptr->d_z(v_psi, j), slice_ptr->d_z(v_K_psi, j)};
            const RSF d_zz_rsf{0.0, 0.0};
            const double rho_psi = rsf_model.rho(rsf, bssn, d_z_rsf, d_zz_rsf);
            e_psi += quad_weights[j] * sphere_area * rho_psi * std::pow(z, D - 2.) * std::pow(chi, 0.5 * (1.0 - D));
        }
    }

    E_phi = bss_mpi::sum(e_phi);
    E_psi = add_real_field ? bss_mpi::sum(e_psi) : 0.0;
}

//compute the 4D Ricci scalar at a grid point on the current slice
double Spacetime::ricci_4_val(int j) const
{
    if (current_slice_ptr == nullptr || n_gridpoints <= 0)
        return 0.0;

    const BSSNSlice& s = *current_slice_ptr;

    const double& chi0 = s.states2[j].bssn.chi;
    const double& h_zz0 = s.states2[j].bssn.h_zz;
    const double& h_ww0 = s.states2[j].bssn.h_ww;
    const double& h_ZZ0 = h_ZZ[j];
    const double& h_WW0 = h_WW[j];
    const double& K0 = s.states2[j].bssn.K;
    const double& A_zz0 = s.states2[j].bssn.A_zz;
    const double& A_ww0 = s.states2[j].bssn.A_ww;
    const double& alpha0 = s.states2[j].bssn.alpha;
    double d_z_K0 = current_slice_ptr->d_z(v_K, j);
    const double& beta0 = s.states2[j].bssn.beta;

    const double K_zz0 = (A_zz0 + h_zz0 * K0 / (D - 1.)) / chi0;
    const double K_ww0 = (A_ww0 + h_ww0 * K0 / (D - 1.)) / chi0;

    double d_tK0 = beta0 * d_z_K0 - chi0 * h_ZZ0 * D_zz_alpha[j] + alpha0 * h_ZZ0 * h_ZZ0 * A_zz0 * A_zz0 
                    + alpha0 * K0 * K0  / (D - 1.)
                    + (D - 2.) * h_WW0 * (alpha0 * A_ww0 * A_ww0 / h_ww0 - chi0 * D_ww_alpha[j]) 
                    + 8. * M_PI * alpha0 * (S[j] + (D - 3.) * rho[j]) / (D - 2.);

    //beta * d_z_K - chi * inv_hzz * D_zz_alpha[j] + alpha * inv_hzz * inv_hzz * A_zz * A_zz + alpha * K * K * invDm1
                            //+ n * inv_hww * (alpha * A_ww * A_ww / h_ww - chi * D_ww_alpha[j]) + 8. * M_PI * alpha * (S[j] + (D - 3.) * rho[j]) / n;
        
    const double ricci_3 = chi0 * (R_zz[j] * h_ZZ0 + (D - 2.) * R_ww[j] * h_WW0);
    const double K_ijK_IJ = chi0 * chi0 * (K_zz0 * K_zz0 * h_ZZ0 * h_ZZ0 + (D - 2.) * K_ww0 * K_ww0 * h_WW0 * h_WW0);
    const double lap_alpha = chi0 * (h_ZZ0 * D_zz_alpha[j] + (D - 2.) * h_WW0 * D_ww_alpha[j]);

    return ricci_3 + K0 * K0 + K_ijK_IJ
        + 2. * (-lap_alpha - d_tK0 + beta0 * d_z_K0) / alpha0;
}

double Spacetime::ricci_4_matter(int j) const
{
    if (current_slice_ptr == nullptr || n_gridpoints <= 0)
        return 0.0;

    double T = S[j] - rho[j];

    return cabs(-8. * M_PI * T);
}

double Spacetime::get_ah_mass() const
{
    if (ah_radius < 0.0 || current_slice_ptr == nullptr || n_gridpoints <= 0)
        return 0.0;

    const double dr_local = R / ((double)n_gridpoints - 1.0);
    const double frac_idx = ah_radius / dr_local - grid_offset;
    const int j0 = std::max(0, std::min((int)std::floor(frac_idx), n_gridpoints - 2));
    const int j1 = j0 + 1;
    const double lam = frac_idx - j0;

    // MPI: only the rank owning j0 has valid state at j0. j1 is either
    // also owned (j0 not at owned_end - 1) or in that rank's right ghost slot
    // (populated by exchange_ghosts). Others contribute 0 and we sum-reduce.
    double contribution = 0.0;
    const bool own_j0 = (j0 >= global_lo && j0 < global_hi);
    if (own_j0)
    {
        const double chi0  = (1.0 - lam) * current_slice_ptr->states2[j0].bssn.chi  + lam * current_slice_ptr->states2[j1].bssn.chi;
        const double h_ww0 = (1.0 - lam) * current_slice_ptr->states2[j0].bssn.h_ww + lam * current_slice_ptr->states2[j1].bssn.h_ww;
        // irreducible mass from the horizon area: (n Omega_n / 16 pi) R_AH^{n-1}, = R_AH / 2 in D = 4
        const double n = D - 2.;
        const double Omega_n = 2. * pow(M_PI, 0.5 * (D - 1.)) / tgamma(0.5 * (D - 1.));
        contribution = n * Omega_n / (16. * M_PI) * pow(std::sqrt(h_ww0 / chi0) * ah_radius, n - 1.);
    }
    return bss_mpi::sum(contribution);
}

//computes hamiltonian and momentum constraints and conformal metric determinant
void Spacetime::compute_diagnostics (BSSNSlice* slice_ptr)
{
    double n = D - 2.;
    dr = R / ((double)n_gridpoints - 1.0);

    Ham_L2 = 0.;
    Mom_L2 = 0.;

    // Owned-range partial accumulation; global L2 norms combined via bss_mpi::sum below.
    const int j_lo_d = std::max(owned_begin(), 0);
    const int j_hi_d = std::min(owned_end(), n_gridpoints - 2);
    for (int j = j_lo_d; j < j_hi_d; j++)
    {
        if (!active_points[j]) //perform no computations on inactive points
            continue;

    double z = (j + grid_offset) * dr;

    //shorthand versions of the BSSN vars on desired slice for convenience
    double chi = slice_ptr->states2[j].bssn.chi;
    const double h_zz = slice_ptr->states2[j].bssn.h_zz;
    const double h_ww = slice_ptr->states2[j].bssn.h_ww;
    const double A_zz = slice_ptr->states2[j].bssn.A_zz;
    const double A_ww = slice_ptr->states2[j].bssn.A_ww;
    const double K = slice_ptr->states2[j].bssn.K;
    //const double beta = slice_ptr->states2[j].bssn.beta;

        if (chi < min_chi)
            chi = min_chi;

        Ham[j] = chi * h_ZZ[j] * R_zz[j] - h_ZZ[j] * h_ZZ[j] * A_zz * A_zz + n * K * K / (D - 1.)
                 + n * (chi * h_WW[j] * R_ww[j] - A_ww * A_ww / (h_ww * h_ww)) - 16. * M_PI * rho[j];
        // pure ADM Hamiltonian constraint (Ricci scalar from the metric only); Ham[j] itself keeps the evolved Gamma
        // variable since slice_mass and the constraint-preserving BCs rely on that form
        Ham_adm[j] = constraints_from_metric ? Ham[j] - ham_corr[j] : Ham[j];


        double A_zzww_diff = (z <= min_z) ? 0. : ((A_zz - A_ww) / z) ; //difference between A_zz and A_ww, vanishing when z = 0 as required

        Mom_Z[j] = - n * d_z(v_K,j) /(D - 1.) + n * h_WW[j] * (A_zzww_diff - 0.5 * h_WW[j] * A_ww * d_z(v_h_ww,j))
               - n * h_WW[j] * chris_Zww[j] * A_zz + h_ZZ[j] * (d_z(v_A_zz,j) - 2. * chris_Zzz[j] * A_zz - (D - 1.) * A_zz * d_z(v_chi,j) / (2. * chi)) - 8. * M_PI * j_z[j];

        det_h[j] =  h_zz * pow(h_ww,n);

        double start_val = (make_tangherlini ? pow(0.5, 1. / (D - 3.)) : 0.);
        double end_val = R * (only_BS_violation ? (r_99 / R) : 0.9);

        //add contribution to L2 norms of Ham/ Mom constraints; z^2 factor for violation on sphere. Ad hoc cutoff radiui to avoid
        //violation being dominated by non-propagating boundary noise/ BH center spike.
       if (z > start_val && z < end_val)
        {
            Ham_L2 += quad_weights[j] * Ham_adm[j] * Ham_adm[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
            Mom_L2 += quad_weights[j] * Mom_Z[j] * Mom_Z[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
        }
    }

    //take square root to get L^2 norms; normalize by 16 * pi * central energy density.
    // MPI: combine per-rank sums-of-squares before sqrt.
    Ham_L2 = sqrt(bss_mpi::sum(Ham_L2)) / (16. * M_PI * rho0_init);
    Mom_L2 = sqrt(bss_mpi::sum(Mom_L2)) / (16. * M_PI * rho0_init);
}

//writes diagnostic quantities to output file. maybe just make appropriate call to write_slice() instead?
// MPI: collective. All ranks pack their owned rows; rank 0 gathers and writes.
void Spacetime:: write_current_slice(std::string file_name)
{
    const BSSNSlice& s = *current_slice_ptr;
    dr = R / ((double)n_gridpoints - 1.0);

    // 16 doubles per point (all state fields; theta packed as 0 if not CCZ4).
    constexpr int NCOL = 16;
    std::vector<double> local(local_N * NCOL);
    for (int k = 0; k < local_N; ++k)
    {
        const int j = owned_begin() + k;
        double* p = &local[k * NCOL];
        p[0]  = s.states2[j].bssn.chi;
        p[1]  = s.states2[j].bssn.h_zz;
        p[2]  = s.states2[j].bssn.h_ww;
        p[3]  = s.states2[j].bssn.A_zz;
        p[4]  = s.states2[j].bssn.A_ww;
        p[5]  = s.states2[j].bssn.K;
        p[6]  = s.states2[j].bssn.c_chris_Z;
        p[7]  = s.states2[j].csf.phi_re;
        p[8]  = s.states2[j].csf.phi_im;
        p[9]  = s.states2[j].csf.K_phi_re;
        p[10] = s.states2[j].csf.K_phi_im;
        p[11] = s.states2[j].bssn.alpha;
        p[12] = s.states2[j].bssn.beta;
        p[13] = use_CCZ4 ? s.theta[j] : 0.0;
        p[14] = s.states2[j].rsf.psi;
        p[15] = s.states2[j].rsf.K_psi;
    }

    std::vector<double> full;
    gather_owned_to_root(local, full, NCOL);

    if (mpi_rank != 0) return;

    std::ofstream data_file{file_name};
    if (!data_file)
    {
        cerr << file_name << " could not be opened for writing!\n";
        bss_mpi::abort_all(1);
    }

    for (int j = 0; j < n_gridpoints; ++j)
    {
        if (!active_points[j]) continue;
        const double* p = &full[j * NCOL];
        data_file << std::setprecision(16) << (j + grid_offset) * dr
                  << "   " << p[0]  << "    " << p[1]  << "    " << p[2]  << "    " << p[3]
                  << "   " << p[4]  << "    " << p[5]  << "    " << p[6]
                  << "    " << p[7] << "    " << p[8]  << "   " << p[9]  << "    " << p[10]
                  << "    " << p[11]<< "    " << p[12] << "    " << p[13]
                  << "   " << p[14] << "    " << p[15] << endl;
    }
}

//writes diagnostic quantities to output file
// MPI: collective. All ranks pack their owned rows into a local buffer and
// track a local max of ricci_4_matter; the row buffer is gathered to rank 0
// for writing, and ricci_4_grid_max is set to the global max via MPI_Allreduce.
void Spacetime:: write_diagnostics(std::string file_name)
{
    int length = (current_slice_ptr->states2).size();
    vector<double>& aux_test = rho; // knob: which auxiliary variable to dump

    const double dr_local = R / ((double)length - 1.0);

    // 8 doubles per row (skip the `r` column — reconstructed from j on rank 0):
    //   Ham, Mom_Z, det_h, aux_test, A, ricci_4_val, r4m, MS mass (last column written only if output_ms_profile)
    constexpr int NCOL = 8;
    std::vector<double> local(local_N * NCOL, 0.0);
    double local_r4max = 0.0;

    const int j_lo = std::max(owned_begin(), 0);
    const int j_hi = std::min(owned_end(), length - 2);
    for (int j = j_lo; j < j_hi; j++)
    {
        if (!active_points[j]) continue;

        const double& phi_re = current_slice_ptr->states2[j].csf.phi_re;
        const double& phi_im = current_slice_ptr->states2[j].csf.phi_im;
        const double A   = sqrt(phi_re * phi_re + phi_im * phi_im);
        const double r4m = ricci_4_matter(j);
        if (r4m > local_r4max) local_r4max = r4m;

        double* p = &local[(j - owned_begin()) * NCOL];
        p[0] = Ham_adm[j];
        p[1] = Mom_Z[j];
        p[2] = det_h[j];
        p[3] = aux_test[j];
        p[4] = A;
        p[5] = ricci_4_val(j);
        p[6] = r4m;
        p[7] = output_ms_profile ? ms_mass_at_point(current_slice_ptr, j) : 0.;
    }
    ricci_4_grid_max = bss_mpi::max(local_r4max);

    std::vector<double> full;
    gather_owned_to_root(local, full, NCOL);

    if (mpi_rank != 0) return;

    std::ofstream data_file{file_name};
    if (!data_file)
    {
        std::cerr << file_name << " could not be opened for writing!\n";
        bss_mpi::abort_all(1);
    }

    for (int j = 0; j < length - 2; ++j)
    {
        if (!active_points[j]) continue;
        const double* p = &full[j * NCOL];
        data_file << std::setprecision(10) << (j + grid_offset) * dr_local
                  << "   " << p[0] << "    " << p[1] << "    " << p[2]
                  << "    " << p[3] << "    " << p[4] << "    " << p[5]
                  << "    " << p[6];
        if (output_ms_profile) data_file << "    " << p[7];
        data_file << endl;
    }
}

//read additional parameters from BSparams.par
void Spacetime::read_parameters(bool quiet)
{
    std::ifstream params{ "BSParams.par" };

    if (!params)
    {
        std::cerr << "Could not open BSParams.par\n";
        abort();
    }

    string current_line{};

    while (getline(params, current_line))
    {
        fill_parameter(current_line, "min_z = ", min_z, quiet);
        fill_parameter(current_line, "min_chi = ", min_chi, quiet);
        fill_parameter(current_line, "sigma_BSSN = ", sigma_BSSN, quiet);
        fill_parameter(current_line, "max_stored_slices = ", max_stored_slices, quiet);
        fill_parameter(current_line, "eta = ", eta, quiet);
        fill_parameter(current_line, "damping_factor = ", damping_factor, quiet);
        fill_parameter(current_line, "write_interval = ", write_interval, quiet);
        fill_parameter(current_line, "write_CN_interval = ", write_CN_interval, quiet);
        fill_parameter(current_line, "BS_resolution_factor = ", BS_resolution_factor, quiet);
        fill_parameter(current_line, "evolve_shift = ", evolve_shift, quiet);
        fill_parameter(current_line, "adaptive_timestep = ", adaptive_timestep, quiet);
        fill_parameter(current_line, "make_tangherlini = ", make_tangherlini, quiet);
        fill_parameter(current_line, "store_A0 = ", store_A0, quiet);
        fill_parameter(current_line, "run_quietly = ", run_quietly, quiet);
        fill_parameter(current_line, "start_time = ", start_time, quiet);
        fill_parameter(current_line, "checkpoint_time = ", checkpoint_time, quiet);
        fill_parameter(current_line, "read_thinshell = ", read_thinshell, quiet);
        fill_parameter(current_line, "cutoff_frac = ", cutoff_frac, quiet);
        fill_parameter(current_line, "only_BS_violation = ", only_BS_violation, quiet);
        fill_parameter(current_line, "run_spacetime_solver = ", run_spacetime_solver, quiet);
        fill_parameter(current_line, "gamma_fac = ", gamma_fac, quiet);
        fill_parameter(current_line, "spatially_varying_BC = ", spatially_varying_BC, quiet);
        fill_parameter(current_line, "use_CCZ4 = ", use_CCZ4, quiet);
        fill_parameter(current_line, "c1 = ", c1, quiet);
        fill_parameter(current_line, "c2 = ", c2, quiet);
        fill_parameter(current_line, "c3 = ", c3, quiet);
        fill_parameter(current_line, "covariant_Z4 = ", covariant_Z4, quiet);
        fill_parameter(current_line, "theta_off = ", theta_off, quiet);
        fill_parameter(current_line, "shock_gauge = ", shock_gauge, quiet);
        fill_parameter(current_line, "shock_fac = ", shock_fac, quiet);
        fill_parameter(current_line, "shock_gauge = ", shock_gauge, quiet);
        fill_parameter(current_line, "one_log_fac = ", one_log_fac, quiet);
        fill_parameter(current_line, "stop_on_migrate = ", stop_on_migrate, quiet);
        fill_parameter(current_line, "cell_centered = ", cell_centered, quiet);
        fill_parameter(current_line, "add_real_field = ", add_real_field, quiet);
        fill_parameter(current_line, "real_amp = ", real_amp, quiet);
        fill_parameter(current_line, "real_sigma = ", real_sigma, quiet);
        fill_parameter(current_line, "real_center = ", real_center, quiet);
        fill_parameter(current_line, "critical_study = ", critical_study, quiet);
        fill_parameter(current_line, "critical_gaussian_start = ", critical_gaussian_start, quiet);
        fill_parameter(current_line, "past_super_time = ", past_super_time, quiet);
        fill_parameter(current_line, "ms_radius_frac = ", ms_radius_frac, quiet);
        fill_parameter(current_line, "ms_bc_radius_frac = ", ms_bc_radius_frac, quiet);
        fill_parameter(current_line, "output_ms_profile = ", output_ms_profile, quiet);
        fill_parameter(current_line, "gamma_constraint_BC = ", gamma_constraint_BC, quiet);
        fill_parameter(current_line, "constraint_preserving_BCs = ", constraint_preserving_BCs, quiet);
        fill_parameter(current_line, "KO_outer_cutoff = ", KO_outer_cutoff, quiet);
        fill_parameter(current_line, "KO_chi_outer_frac = ", KO_chi_outer_frac, quiet);
        fill_parameter(current_line, "constraints_from_metric = ", constraints_from_metric, quiet);
        fill_parameter(current_line, "k1_theta_inner_factor = ", k1_theta_inner_factor, quiet);
        fill_parameter(current_line, "inner_damp_frac = ", inner_damp_frac, quiet);
        fill_parameter(current_line, "shift_taper_frac = ", shift_taper_frac, quiet);
        fill_parameter(current_line, "critical_state = ", critical_state, quiet);
        fill_parameter(current_line, "critical_eps = ", critical_eps, quiet);
        fill_parameter(current_line, "lapse_thresh = ", lapse_thresh, quiet);
        fill_parameter(current_line, "hi_guess = ", hi_guess, quiet);
        fill_parameter(current_line, "lo_guess = ", lo_guess, quiet);
        fill_parameter(current_line, "sub_min_time = ", sub_min_time, quiet);
        fill_parameter(current_line, "subcritical_time = ", subcritical_time, quiet);
        fill_parameter(current_line, "do_ah_search = ", do_ah_search, quiet);
        fill_parameter(current_line, "AH_search_steps = ", AH_search_steps, quiet);
        fill_parameter(current_line, "fill_subcritical = ", fill_subcritical, quiet);
        fill_parameter(current_line, "A_thresh = ", A_thresh, quiet);
        fill_parameter(current_line, "dump_time = ", dump_time, quiet);

        fill_param_array(current_line, "refinement_points = ", refinement_points, quiet);

    }
    cout << sigma_BSSN << eta << endl;
}

// Build the RK4 update slice with optional preservation of outer refinement levels.
void Spacetime::get_update_slice(const BSSNSlice& s1,
                                 const BSSNSlice& s2,
                                 const BSSNSlice& s3,
                                 const BSSNSlice& s4,
                                 BSSNSlice& update_slice,
                                 double dt,
                                 int skipped_levels)
{
    const int N = static_cast<int>(s1.states2.size());

    // Ensure update_slice structural fields are consistent
    update_slice.R = s1.R;
    update_slice.has_BH = s1.has_BH;
    update_slice.refinement_points = s1.refinement_points;
    update_slice.use_CCZ4 = s1.use_CCZ4;
    update_slice.states2.resize(N);
    if (update_slice.use_CCZ4)
        update_slice.theta.resize(N);

    const bool do_preserve = (skipped_levels > 0 && !refinement_levels.empty());
    const int max_ref = do_preserve ? refinement_levels.back() : 0;

    // MPI: compute RK4 combination only on this rank's owned range. In serial
    // mode owned range covers the full grid, so this matches the previous
    // (combo = (dt/6)*(s1 + 2*s2 + 2*s3 + s4)) semantics exactly.
    const int lo = std::max(owned_begin(), 0);
    const int hi = std::min(owned_end(), N);
    const double c = dt / 6.0;

    for (int j = lo; j < hi; ++j)
    {
        if (do_preserve && refinement_levels[j] > max_ref - skipped_levels)
            continue; // preserve existing update_slice[j]

        update_slice.states2[j] = c * (s1.states2[j] + 2.0 * s2.states2[j] + 2.0 * s3.states2[j] + s4.states2[j]);
        if (update_slice.use_CCZ4)
            update_slice.theta[j] = c * (s1.theta[j] + 2.0 * s2.theta[j] + 2.0 * s3.theta[j] + s4.theta[j]);
    }
}

void Spacetime::slice_axpy(BSSNSlice& result, const BSSNSlice& base, double c, const BSSNSlice& rhs)
{
    const int N = static_cast<int>(base.states2.size());

    // Structural fields follow base; sizes ensured here so first-call case
    // (result starts with empty states2) still works.
    result.R = base.R;
    result.has_BH = base.has_BH;
    result.refinement_points = base.refinement_points;
    result.use_CCZ4 = base.use_CCZ4;
    result.states2.resize(N);
    if (result.use_CCZ4)
        result.theta.resize(N);

    const int lo = std::max(owned_begin(), 0);
    const int hi = std::min(owned_end(), N);

    for (int j = lo; j < hi; ++j)
    {
        result.states2[j] = base.states2[j] + c * rhs.states2[j];
        if (result.use_CCZ4)
            result.theta[j] = base.theta[j] + c * rhs.theta[j];
    }
}

//halves number of gridpoints while keeping R. Should only be called before evolving!
void Spacetime::halve_resolution()
{
    int n_old = n_gridpoints;

    // Downsample the full composite state (geometry + matter)
    vector<State> old_states(n_old);

    //cut off one point at end to get odd # of gridpoints if needed
    if (n_gridpoints % 2 == 0)
        {
            n_gridpoints -= 1;
            R *= (n_gridpoints) / (n_old);
        }

    for (int j = 0; j < n_gridpoints; j++) // fill old array with states
        old_states[j] = slices[0].states2[j];

    n_gridpoints = (n_gridpoints + 1) / 2;

    // Write back every other point and then shrink
    for (int k = 0; k < n_gridpoints; k++)
        slices[0].states2[k] = old_states[2 * k];
    slices[0].states2.resize(n_gridpoints);

    dr = R / (n_gridpoints - 1.0);
}

void Spacetime::fill_active_points()
{
    //first deactivate all points in case of empty refinement_levels
    for (int j = 0; j < n_gridpoints; j++)
        active_points[j] = 0;

    int step_size = 1; //amount we step forward by in activating points
    unsigned int refinement_layers_crossed = 0;

    //now fill in points, halving the number of points filled in every time we pass a refinement layer
    for (int k = 0; k < n_gridpoints; k += step_size)
        {
            active_points[k] = 1;

            //when we cross refinement layer, double step size unless we're already on the last one
            if (refinement_layers_crossed < refinement_points.size() && k > refinement_points[refinement_layers_crossed]  )
            {
                refinement_layers_crossed++;
                step_size *= 2;
            }
        }
    //for (int j = 0; j < n_gridpoints; j++)
        //cout << active_points[j] << endl;
}

void Spacetime::fill_refinement_levels()
{
    refinement_levels.resize(n_gridpoints);
    for (int j = 0; j < n_gridpoints; j++)
        refinement_levels[j] = slices[0].get_refinement_level(j, refinement_points);
}

//corrects the initial field momentum in such a way that constraints are satisfied. Must have filled in initial slice already.
//seems to only work for negative initial perturbation, at least for early tests... now deprecated
void Spacetime::fix_initial_field_mom()
{
    compute_auxiliary_quantities(&slices[0], 0);
    compute_diagnostics(&slices[0]);
    bool failed = 0;

   // vector<double> old_momenta(n_gridpoints);

    for (int j = 0; j < n_gridpoints; j++)
    {
    //old_momenta[j] = slices[0].states2[j].csf.K_phi_im;
    double correction_sq = slices[0].states2[j].csf.K_phi_im * slices[0].states2[j].csf.K_phi_im + Ham[j] / (32. * M_PI);

        if (correction_sq < 0. && failed == 0)
            {
                failed = 1;
                cout << "WARNING: trick to adjust initial field momentum may have failed starting on gridpoint " << j << endl;
            }

       if (correction_sq >= 0.)
           slices[0].states2[j].csf.K_phi_im = - sqrt(correction_sq);
    }
}

void Spacetime::compute_dtK(int time_index)
{
    if (time_index == 0)
    {
        cout<< "WARNING: called compute_dtK on timestep 0" << endl;
        return;
    }

    BSSNSlice& slice_current = slices[time_index];
    BSSNSlice& slice_last = slices[time_index - 1];

    // MPI: partial sum-of-squares over owned range, then combine.
    double dtK_local = 0.0;
    for (int j = owned_begin(); j < owned_end(); j++)
    {
        if (!active_points[j]) //inactive points are never updated
            continue;

        double z = (j + 0.5) * dr;

        if (z < R * (only_BS_violation ? (r_99 / R) : 0.9))
        {
            double local_dtK =  (slice_current.states2[j].bssn.K - slice_last.states2[j].bssn.K ) / dt;
            dtK_local += quad_weights[j] * local_dtK * local_dtK * z * z   * pow(slice_current.states2[j].bssn.chi, -1.5);
        }
    }
    dtK_L2 = sqrt(bss_mpi::sum(dtK_local));
    dtK_L2 = abs(dtK_L2 - 1);
}

void Spacetime::prepare_ham_solve()  //EXPERIMENTAL: try to return to pure isotropic coords + solve Ham constraint mid-run
{

}

//helper function to set all temp array sizes to n_gridpoints and initialize some necessary variables
void Spacetime::resize_temp_arrays()
{
    h_ZZ.resize(n_gridpoints);
    h_WW.resize(n_gridpoints);

    chris_Zww.resize(n_gridpoints);
    chris_Zzz.resize(n_gridpoints);
    aux_constraint.resize(n_gridpoints);
    D_zz_alpha.resize(n_gridpoints);
    D_ww_alpha.resize(n_gridpoints);
    D_zz_alpha_TF.resize(n_gridpoints);
    D_ww_alpha_TF.resize(n_gridpoints);
    R_zz.resize(n_gridpoints);
    R_ww.resize(n_gridpoints);
    R_zz_TF.resize(n_gridpoints);
    R_ww_TF.resize(n_gridpoints);

    rho.resize(n_gridpoints);
    j_z.resize(n_gridpoints);
    S_zz.resize(n_gridpoints);
    S_ww.resize(n_gridpoints);
    S.resize(n_gridpoints);
    S_zz_TF.resize(n_gridpoints);
    S_ww_TF.resize(n_gridpoints);

    Ham.resize(n_gridpoints);
    Ham_adm.resize(n_gridpoints);
    ham_corr.resize(n_gridpoints);
    Mom_Z.resize(n_gridpoints);
    det_h.resize(n_gridpoints);

    if (use_CCZ4)
    {
        theta_Z.resize(n_gridpoints);
    }

    dtK_L2 = 0;

    last_active_j = n_gridpoints - 1; //find last active gridpoint
    while (!active_points[last_active_j])
        last_active_j--;
    BSSNSlice::outer_last_j = last_active_j; // used by the one-sided outer stencils of the constraint-preserving BCs
    BSSNSlice::level_table = &refinement_levels; // per-index refinement level for d_z (avoids recomputing it per call)
    deriv_cache_slice = nullptr; // stencil geometry changed: drop the cached derivative stencil

    fill_quad_weights();
    build_inactive_stencils();
}

//Makes the refinement layout consistent with the final grid. refinement_points are given as indices into the grid before
//initialize() trims it (by two points and then by cutoff_frac), and the level/active-point arrays were built for that
//untrimmed grid. A refinement point is kept only if at least 6 points of the coarser level fit beyond it on the trimmed
//grid (needed by the boundary interpolation and by the outer boundary stencil); it and any later ones are dropped
//otherwise. The active-point and level arrays are then rebuilt for the final grid size.
void Spacetime::finalize_refinement()
{
    std::vector<int> kept;
    for (size_t k = 0; k < refinement_points.size(); k++)
    {
        const int coarse_step = 1 << (k + 1);
        if (refinement_points[k] + 1 + 6 * coarse_step <= n_gridpoints - 1)
            kept.push_back(refinement_points[k]);
        else
            break;
    }

    if (kept.size() < refinement_points.size() && mpi_rank == 0)
    {
        cout << "WARNING: refinement point(s)";
        for (size_t k = kept.size(); k < refinement_points.size(); k++)
            cout << " " << refinement_points[k];
        cout << " do not fit on the trimmed grid of " << n_gridpoints << " points and will be ignored" << endl;
    }
    refinement_points = kept;

    active_points.assign(n_gridpoints, false);
    fill_active_points();
    fill_refinement_levels();

    //refinement boundary points: points of the refinement pattern whose centred 5-point stencil (at their level's spacing)
    //would reach a point that is not part of the pattern. They cannot be evolved consistently (doing so with interpolated
    //neighbours is unstable for the Proca system), so they are removed from the active set and filled by interpolation
    //(see build_inactive_stencils / fill_inactive_points). There is one such point per refinement boundary.
    refinement_boundary_points.clear();
    if (!refinement_points.empty())
    {
        int last = n_gridpoints - 1;
        while (last > 0 && !active_points[last]) last--;

        for (int j = 0; j <= last; j++)
        {
            if (!active_points[j]) continue;
            const int rf = 1 << (refinement_levels[j] - 1);
            for (int m : {-2, -1, 1, 2})
            {
                const int i = std::abs(j + m * rf);
                if (i <= last && !active_points[i])
                {
                    refinement_boundary_points.push_back(j);
                    break;
                }
            }
        }
        for (int j : refinement_boundary_points)
            active_points[j] = false;

        if (!run_quietly && mpi_rank == 0)
        {
            cout << "Refinement boundary points (interpolated, not evolved):";
            for (int j : refinement_boundary_points) cout << " " << j;
            cout << endl;
        }
    }

    slices[0].refinement_points = refinement_points;
    slices[0].active_points = active_points;
}

//Precomputes, for every inactive point strictly inside the last active point, a 4-point Lagrange interpolation stencil
//over the nearest active points (two on each side where possible, otherwise three below and one above); refinement
//boundary points get a symmetric stencil instead.
void Spacetime::build_inactive_stencils()
{
    inactive_stencils.clear();

    //refinement boundary points: symmetric 4-point stencil j - 3s, j - s, j + s, j + 3s over active points, with the
    //smallest such spacing s (weights -1/16, 9/16, 9/16, -1/16)
    std::vector<bool> is_boundary(n_gridpoints, false);
    for (int j : refinement_boundary_points)
    {
        for (int s = 1; s <= (1 << 12); s *= 2)
        {
            const int nodes[4] = {j - 3 * s, j - s, j + s, j + 3 * s};
            bool usable = true;
            for (int i : nodes)
                if (i < 0 || i > last_active_j || !active_points[i]) { usable = false; break; }
            if (!usable) continue;

            InactiveStencil st;
            st.j = j;
            const double w[4] = {-1. / 16., 9. / 16., 9. / 16., -1. / 16.};
            for (int m = 0; m < 4; m++) { st.idx[m] = nodes[m]; st.w[m] = w[m]; }
            inactive_stencils.push_back(st);
            is_boundary[j] = true;
            break;
        }
    }

    std::vector<int> act; //indices of active points, in order
    for (int j = 0; j <= last_active_j; j++)
        if (active_points[j])
            act.push_back(j);

    int a = 0; //index into act of the last active point below j
    for (int j = 0; j < last_active_j; j++)
    {
        if (active_points[j])
        {
            while (a + 1 < (int)act.size() && act[a + 1] <= j) a++;
            continue;
        }

        if (is_boundary[j]) continue; //already has its symmetric stencil

        //act[a] < j < act[a + 1]
        int lo = a - 1;
        if (a + 2 >= (int)act.size()) lo = a - 2; //no second point above: use three below
        if (lo < 0) lo = 0;
        if (lo + 3 >= (int)act.size()) continue; //not enough active points (should not happen)

        InactiveStencil st;
        st.j = j;
        for (int m = 0; m < 4; m++) st.idx[m] = act[lo + m];
        for (int m = 0; m < 4; m++)
        {
            double w = 1.;
            for (int l = 0; l < 4; l++)
                if (l != m) w *= (double)(j - st.idx[l]) / (double)(st.idx[m] - st.idx[l]);
            st.w[m] = w;
        }
        inactive_stencils.push_back(st);
    }
}

//Overwrites every inactive point of the slice by interpolation from the active points (no-op on a uniform grid)
void Spacetime::fill_inactive_points(BSSNSlice& slice)
{
    for (const InactiveStencil& st : inactive_stencils)
    {
        slice.states2[st.j] = st.w[0] * slice.states2[st.idx[0]] + st.w[1] * slice.states2[st.idx[1]]
                            + st.w[2] * slice.states2[st.idx[2]] + st.w[3] * slice.states2[st.idx[3]];

        if (use_CCZ4 && slice.theta.size() == slice.states2.size())
            slice.theta[st.j] = st.w[0] * slice.theta[st.idx[0]] + st.w[1] * slice.theta[st.idx[1]]
                              + st.w[2] * slice.theta[st.idx[2]] + st.w[3] * slice.theta[st.idx[3]];
    }
}


//---------------------------------------------------------------------------------------------------------------------
// Berger-Oliger subcycling for the coarsening-outward refinement (adaptive_timestep = 1)
//---------------------------------------------------------------------------------------------------------------------
//Builds the per-level point sets. Level L owns the active points of its region (rp_{L-1}, rp_L]; its buffer is the
//BO_BUFFER coarse-lattice points just inside the finer region (stepped provisionally with the coarse level so that the
//owned points near the interface see a complete, self-consistent coarse stencil; overwritten by the fine solution at
//the end of the coarse step); its zone is the first few owned points, which the finer level reads through the Hermite
//interpolant. Contamination from the buffer's inner edge travels 3 lattice points per RHS evaluation (KO stencil), i.e.
//15 points for the 4 RK stages plus the extra RHS for the interpolant, so BO_BUFFER = 16 keeps the owned points exact.
void Spacetime::bo_setup()
{
    const int BO_BUFFER = 16, BO_ZONE = 5;
    bo_levels.clear();
    const int n_lev = (int)refinement_points.size() + 1;
    for (int L = 1; L <= n_lev; L++)
    {
        BOLevel lev;
        lev.level = L;
        lev.res_fac = 1 << (L - 1);
        lev.region_lo = (L == 1) ? 0 : refinement_points[L - 2] + 1;
        lev.region_hi = (L == n_lev) ? last_active_j : refinement_points[L - 1];
        for (int j = lev.region_lo; j <= lev.region_hi; j++)
            if (active_points[j]) lev.owned.push_back(j);
        if (lev.owned.empty())
        {
            if (mpi_rank == 0) cout << "ERROR: refinement level " << L << " owns no active points" << endl;
            bss_mpi::abort_all(1);
        }
        if (L > 1)
        {
            const int j0 = lev.owned.front();
            for (int m = 1; m <= BO_BUFFER; m++)
            {
                const int j = j0 - m * lev.res_fac;
                //the innermost buffer point still needs 3 lattice points (KO stencil) below it that are not mirrored
                if (j - 3 * lev.res_fac < 0 || j < bo_levels[L - 2].region_lo || !active_points[j])
                {
                    if (mpi_rank == 0)
                        cout << "ERROR: adaptive_timestep needs refinement region " << L - 1 << " to hold at least "
                             << (BO_BUFFER + 3) * lev.res_fac << " points for the level-" << L << " buffer" << endl;
                    bss_mpi::abort_all(1);
                }
                lev.buffer.push_back(j);
            }
            for (int m = 0; m < BO_ZONE && m < (int)lev.owned.size(); m++)
                lev.zone.push_back(lev.owned[m]);
            const size_t nz = lev.zone.size();
            lev.y_old.resize(nz); lev.f_old.resize(nz); lev.y_new.resize(nz); lev.f_new.resize(nz);
            lev.th_old.assign(nz, 0.); lev.fth_old.assign(nz, 0.); lev.th_new.assign(nz, 0.); lev.fth_new.assign(nz, 0.);
        }
        lev.comp = lev.buffer;
        lev.comp.insert(lev.comp.end(), lev.owned.begin(), lev.owned.end());
        //inactive points that the derivative (+-2) and KO (+-3) stencils of this level's points can read: these must be
        //re-interpolated at every stage (the uniform-step scheme refreshes all inactive points; here only the ones needed)
        std::vector<bool> needed(n_gridpoints, false);
        for (size_t m = 0; m < lev.comp.size(); m++)
        {
            const int j = lev.comp[m];
            const int rf = (m < lev.buffer.size()) ? lev.res_fac : (1 << (refinement_levels[j] - 1));
            for (int d = -3; d <= 3; d++)
            {
                int i = j + d * rf;
                if (i < 0) i = -i - BSSNSlice::origin_mirror_shift;
                if (i >= 0 && i <= last_active_j && !active_points[i]) needed[i] = true;
            }
        }
        for (size_t k = 0; k < inactive_stencils.size(); k++)
            if (needed[inactive_stencils[k].j]) lev.bdry_stencils.push_back((int)k);
        bo_levels.push_back(lev);
    }
    if (!run_quietly && mpi_rank == 0)
    {
        cout << "Berger-Oliger subcycling over " << n_lev << " levels: dt_fine = " << dt_fine << ", dt_coarse = " << dt << endl;
        for (const BOLevel& lev : bo_levels)
            cout << "  level " << lev.level << ": spacing " << lev.res_fac << " dr, points [" << lev.region_lo << ", " << lev.region_hi
                 << "], " << lev.owned.size() << " owned, " << lev.buffer.size() << " buffer, " << lev.zone.size() << " zone, "
                 << lev.bdry_stencils.size() << " interpolated point(s)" << endl;
    }
}

//Per-point weights that are fixed in time: the Gamma-driver shift taper (constraint-preserving BCs) and the CCZ4 outer
//taper of the chi dissipation. Same expressions as were previously evaluated at every RHS call.
void Spacetime::fill_taper_weights()
{
    shift_taper_w.assign(n_gridpoints, 1.);
    chi_ko_w.assign(n_gridpoints, 1.);
    const double z_out = (n_gridpoints - 1 + grid_offset) * dr;
    for (int j = 0; j < n_gridpoints; j++)
    {
        const double zj = (j + grid_offset) * dr;
        shift_taper_w[j] = 0.5 * (1. - tanh((zj - shift_taper_frac * z_out) / (0.04 * z_out)));
        chi_ko_w[j] = 0.5 * (1. - tanh((zj - KO_chi_outer_frac * z_out) / (0.05 * z_out)));
    }
}

//Fills W at the zone points of (coarse) level L with the cubic Hermite interpolant of that level's step in time
void Spacetime::bo_fill_zone(int L, double tau, BSSNSlice& W)
{
    const BOLevel& lev = bo_levels[L - 1];
    double th = (tau - lev.t_old) / lev.dt;
    th = std::max(0., std::min(1., th));
    const double h00 = 2. * th * th * th - 3. * th * th + 1., h10 = th * th * th - 2. * th * th + th;
    const double h01 = -2. * th * th * th + 3. * th * th,     h11 = th * th * th - th * th;
    for (size_t m = 0; m < lev.zone.size(); m++)
    {
        const int j = lev.zone[m];
        W.states2[j] = h00 * lev.y_old[m] + (h10 * lev.dt) * lev.f_old[m] + h01 * lev.y_new[m] + (h11 * lev.dt) * lev.f_new[m];
        if (use_CCZ4)
            W.theta[j] = h00 * lev.th_old[m] + h10 * lev.dt * lev.fth_old[m] + h01 * lev.th_new[m] + h11 * lev.dt * lev.fth_new[m];
    }
}

void Spacetime::bo_fill_bdry(const BOLevel& lev, BSSNSlice& W)
{
    for (int k : lev.bdry_stencils)
    {
        const InactiveStencil& st = inactive_stencils[k];
        W.states2[st.j] = st.w[0] * W.states2[st.idx[0]] + st.w[1] * W.states2[st.idx[1]]
                        + st.w[2] * W.states2[st.idx[2]] + st.w[3] * W.states2[st.idx[3]];
        if (use_CCZ4)
            W.theta[st.j] = st.w[0] * W.theta[st.idx[0]] + st.w[1] * W.theta[st.idx[1]]
                          + st.w[2] * W.theta[st.idx[2]] + st.w[3] * W.theta[st.idx[3]];
    }
}

//Advances level L (and, recursively, all finer levels) from t_L to t_L + dt_L in place on W. On entry W holds the
//level's points at t_L and the finer levels at t_L; coarser levels are already at later times with their Hermite data set.
void Spacetime::bo_advance_level(int L, double t_L, double dt_L, BSSNSlice& W)
{
    BOLevel& lev = bo_levels[L - 1];
    const int n_lev = (int)bo_levels.size();
    const bool coarsest = (L == n_lev);
    const std::vector<int>& P = lev.comp;
    const size_t nP = P.size();
    const double n = D - 2.;
    current_slice_ptr = &W; //Spacetime::d_z / d_zz (used by slice_rhs) differentiate the slice current_slice_ptr points to

    //fine values at the buffer points (restored before the finer level steps from them)
    std::vector<State> saved_buf(lev.buffer.size());
    std::vector<double> saved_buf_th(lev.buffer.size(), 0.);
    for (size_t m = 0; m < lev.buffer.size(); m++)
    {
        saved_buf[m] = W.states2[lev.buffer[m]];
        if (use_CCZ4) saved_buf_th[m] = W.theta[lev.buffer[m]];
    }

    std::vector<State> y0(nP);
    std::vector<double> th0(nP, 0.);
    for (size_t m = 0; m < nP; m++)
    {
        y0[m] = W.states2[P[m]];
        if (use_CCZ4) th0[m] = W.theta[P[m]];
    }

    //buffer points use the coarse stencil spacing
    if (L > 1)
    {
        BSSNSlice::forced_res_fac = lev.res_fac;
        BSSNSlice::forced_below_j = lev.region_lo - 1;
    }

    auto set_stage = [&](const BSSNSlice& k, double c, double tau)
    {
        for (size_t m = 0; m < nP; m++)
        {
            W.states2[P[m]] = y0[m] + c * k.states2[P[m]];
            if (use_CCZ4) W.theta[P[m]] = th0[m] + c * k.theta[P[m]];
        }
        if (!coarsest) bo_fill_zone(L + 1, tau, W);
        bo_fill_bdry(lev, W);
    };

    if (!coarsest) bo_fill_zone(L + 1, t_L, W);
    bo_fill_bdry(lev, W);
    BSSNSlice &k1 = lev.k1, &k2 = lev.k2, &k3 = lev.k3, &k4 = lev.k4; //preallocated; only the level's points are (re)written
    slice_rhs_into(&W, k1, 0, &P, coarsest);
    set_stage(k1, 0.5 * dt_L, t_L + 0.5 * dt_L);
    slice_rhs_into(&W, k2, 0, &P, coarsest);
    set_stage(k2, 0.5 * dt_L, t_L + 0.5 * dt_L);
    slice_rhs_into(&W, k3, 0, &P, coarsest);
    set_stage(k3, dt_L, t_L + dt_L);
    slice_rhs_into(&W, k4, 0, &P, coarsest);

    //RK4 update, traceless A and chi floor on the level's points
    for (size_t m = 0; m < nP; m++)
    {
        const int j = P[m];
        W.states2[j] = y0[m] + (dt_L / 6.) * (k1.states2[j] + 2.0 * k2.states2[j] + 2.0 * k3.states2[j] + k4.states2[j]);
        if (use_CCZ4) W.theta[j] = th0[m] + (dt_L / 6.) * (k1.theta[j] + 2.0 * k2.theta[j] + 2.0 * k3.theta[j] + k4.theta[j]);
        State& st = W.states2[j];
        const double A = st.bssn.A_zz / st.bssn.h_zz + n * st.bssn.A_ww / st.bssn.h_ww;
        st.bssn.A_zz -= st.bssn.h_zz * A / (D - 1.);
        st.bssn.A_ww -= st.bssn.h_ww * A / (D - 1.);
        if (st.bssn.chi < min_chi) st.bssn.chi = min_chi;
    }
    bo_fill_bdry(lev, W);

    if (L > 1)
    {
        //Hermite data for the finer level: one extra RHS evaluation at the new time, only at the zone points
        BSSNSlice& fN = lev.fN;
        slice_rhs_into(&W, fN, 0, &lev.zone, false);
        lev.t_old = t_L;
        lev.dt = dt_L;
        for (size_t m = 0; m < lev.zone.size(); m++)
        {
            const int j = lev.zone[m];
            lev.y_new[m] = W.states2[j];
            lev.f_new[m] = fN.states2[j];
            lev.f_old[m] = k1.states2[j];
            if (use_CCZ4) { lev.th_new[m] = W.theta[j]; lev.fth_new[m] = fN.theta[j]; lev.fth_old[m] = k1.theta[j]; }
        }
        for (size_t m = 0; m < nP; m++) //y_old at the zone points
        {
            auto it = std::find(lev.zone.begin(), lev.zone.end(), P[m]);
            if (it != lev.zone.end())
            {
                const size_t mz = it - lev.zone.begin();
                lev.y_old[mz] = y0[m];
                if (use_CCZ4) lev.th_old[mz] = th0[m];
            }
        }
        BSSNSlice::forced_res_fac = 0;
        BSSNSlice::forced_below_j = -1;

        //restore the fine solution at the buffer and advance the finer level by two half-steps
        for (size_t m = 0; m < lev.buffer.size(); m++)
        {
            W.states2[lev.buffer[m]] = saved_buf[m];
            if (use_CCZ4) W.theta[lev.buffer[m]] = saved_buf_th[m];
        }
        bo_advance_level(L - 1, t_L, 0.5 * dt_L, W);
        bo_advance_level(L - 1, t_L + 0.5 * dt_L, 0.5 * dt_L, W);

        //the finer level has left the zone at the interpolant's end value (= y_new); restore exactly
        for (size_t m = 0; m < lev.zone.size(); m++)
        {
            W.states2[lev.zone[m]] = lev.y_new[m];
            if (use_CCZ4) W.theta[lev.zone[m]] = lev.th_new[m];
        }
    }
}

//radial quadrature weights for integrals over the slice (masses, charges, energies, constraint norms). With mesh
//refinement only every 2^(level-1)-th point is active, so each active point j is weighted by half the distance
//between its active neighbours, dr * (j_next - j_prev) / 2 (the origin is mirrored, j_prev = -1, and the last active
//point uses its inner spacing). Inactive points get weight 0. On a uniform grid every weight is exactly dr.
void Spacetime::fill_quad_weights()
{
    quad_weights.assign(n_gridpoints, 0.);

    for (int j = 0; j < n_gridpoints; j++)
    {
        if (!active_points[j])
            continue;

        int prev = j - 1;
        while (prev >= 0 && !active_points[prev])
            prev--;

        int next = j + 1;
        while (next < n_gridpoints && !active_points[next])
            next++;

        if (next >= n_gridpoints) //last active point
            next = j + (j - prev);

        quad_weights[j] = 0.5 * (next - prev) * dr;
    }
}

//add perturbation directly to spacetime field profile, thus skipping BS step
void Spacetime::add_spacetime_pert(double a, double k, double center)
{
    BSSNSlice& s = *current_slice_ptr;
    double k2 = k * k;
    double dr = R / ((double)n_gridpoints - 1.0);

    for (int j = 0; j < n_gridpoints; j++)
    {
        double phase = std::arg( std::complex<double>(s.states2[j].csf.phi_re, s.states2[j].csf.phi_im));
        double r = j * dr;

        s.states2[j].csf.phi_re += cos(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states2[j].csf.phi_im += sin(phase) * a * exp ( -pow (r - center, 2.) / k2);
        s.states2[j].csf.K_phi_re += 0.5 * sin(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states2[j].bssn.alpha);
        s.states2[j].csf.K_phi_im -= 0.5 * cos(phase) * omega * a * exp ( -pow (r - center, 2.) / k2) / (2. * s.states2[j].bssn.alpha);
    }
    solve_initial_ham(run_quietly);
}

// Add a massless real scalar Gaussian to the provided slice.
// psi(z) = real_amp * exp(-0.5 * ((z - real_center)/real_sigma)^2), K_psi(z) = 0
void Spacetime::add_real_gaussian(BSSNSlice& slice)
{
    // Guard: if disabled or degenerate sigma, do nothing
    if (!add_real_field || real_sigma <= 0.0 || real_amp == 0.0)
        return;

    const int N = n_gridpoints; // use current spacetime grid size
    const double dr_local = slice.R / ((double)N - 1.0);
    const double offset = 0.5 * (double)cell_centered; // consistent with grid_offset

    const double inv_sigma2 = 1.0 / (real_sigma * real_sigma);

    for (int j = 0; j < N; ++j)
    {
        const double z = (j + offset) * dr_local;
        const double dz = z - real_center;
        const double gauss = real_amp * exp(-0.5 * dz * dz * inv_sigma2);

        slice.states2[j].rsf.psi = gauss;
        slice.states2[j].rsf.K_psi= 0.0;
    }
}

//read data from BosonStar to spacetime and construct initial time slice
void Spacetime::initialize(BosonStar& boson_star, bool skip_read)
{
    wave_mode = 0; //make 1 for MR testing purposes only
    ah_radius = -1.;
    D = boson_star.D;

    //inherit parameters from BS
    n_gridpoints = boson_star.n_gridpoints;
    R = boson_star.R;
    courant_factor = boson_star.courant_factor;
    stop_time = boson_star.stop_time;
    mu = boson_star.mu;
    lambda = boson_star.lambda;
    sigma = boson_star.sigma;
    solitonic = boson_star.solitonic;
    omega = boson_star.omega;
    isotropic = boson_star.isotropic;
    M = boson_star.M;
    r_99 = boson_star.r_99;
    BS_perturbed = boson_star.perturb;

    if (!skip_read)
    {
        // Initialize defaults prior to reading params for backwards compatibility
        gamma_fac = 0.75; // 3/4 by default
        spatially_varying_BC = 1;
        one_log_fac = 2.0;
    }
    critical_state = 2;

    if (!skip_read)
    {   
        // Default subcritical_time if not provided
        subcritical_time = stop_time;
        refinement_points = {};
        read_parameters();
    }

    if (critical_study)
        do_ah_search = 1;
    


    //cout << refinement_points.size() << endl;
    if (!refinement_points.empty() && refinement_points[0] <= 0) // disable refinement if sentinel present
        refinement_points.clear();

    // MVP MPI scope: multi-rank runs only support a uniform grid with no
    // adaptive timestep. The ghost-zone width and load balance change with
    // refinement, and outer refinement levels get skipped on the last rank
    // when adaptive_timestep is on; both are out of scope for the first pass.
    if (mpi_size > 1 && (!refinement_points.empty() || adaptive_timestep))
    {
        if (mpi_rank == 0)
            cerr << "ERROR: MPI multi-rank mode requires a uniform grid: "
                 << "refinement_points must be empty (or its first entry <= 0) "
                 << "and adaptive_timestep must be 0." << endl;
        bss_mpi::abort_all(1);
    }

    if (read_thinshell)
        BS_resolution_factor = 1;

    if ((BS_resolution_factor & (BS_resolution_factor - 1)) != 0 || BS_resolution_factor <= 0)
        {
            cerr << "ERROR: BS_resolution_factor must be a power of 2" << endl;
            exit(1);
        }

    active_points.resize(n_gridpoints);
    fill_active_points();
    fill_refinement_levels();

    // Cache original BosonStar grid so we can restore after temporary high-res solve
    const int bs_n_gridpoints_orig = boson_star.n_gridpoints;
    bool bumped_bs_resolution = false;

    //solve BS at higher resolution and read in data to first slice, if not starting in other mode (checkpoint / thinshell read / gaussian start)
    if (BS_resolution_factor > 1 && start_time == 0 && !read_thinshell && !boson_star.gaussian_start)
    {
        if (mpi_rank == 0)
        {
            boson_star.n_gridpoints = boson_star.n_gridpoints * BS_resolution_factor - BS_resolution_factor + 1;

            if (!run_quietly)
                cout << " \nRe-solving BS at higher resolution with " << boson_star.n_gridpoints <<  endl;

            boson_star.solve(1);
            boson_star.write_field();
            boson_star.fill_isotropic_arrays();
            boson_star.write_isotropic();
        }
        boson_star.broadcast();

        //using higher-res omega seems like an improvement
        omega = boson_star.omega;
        bumped_bs_resolution = true;
    }

    if (read_thinshell && start_time == 0)
    {
        boson_star.read_thinshell();
        boson_star.write_field();

        if (boson_star.isotropic)
        {
            boson_star.fill_isotropic_arrays();
            boson_star.write_isotropic();
        }

        R = boson_star.R;
        n_gridpoints = boson_star.n_gridpoints;
        omega = boson_star.omega;

        cout << "R = " << boson_star.R << endl;
    }

    //add perturbation to "standard" BS
    if (!read_thinshell && !boson_star.gaussian_start && boson_star.perturb)
    {
        if (mpi_rank == 0)
        {
            boson_star.add_perturbation(boson_star.perturb_amp, boson_star.perturb_spread, boson_star.perturb_center);
            boson_star.fill_given_A(omega);
        }
        boson_star.broadcast();
    }

    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    int num_timesteps = ceil(stop_time / dt);

    slices.resize(std::min(num_timesteps + 1, max_stored_slices));
    slices[0].active_points = active_points;
    slices[0].refinement_points = refinement_points;
    slices[0].use_CCZ4 = use_CCZ4;

    //cout << "About to read" << endl;

    if (start_time == 0)
        slices[0].read_BS_data(boson_star, BS_resolution_factor, isotropic, cell_centered);

    // Restore BosonStar gridpoints so repeated initialize() calls (e.g., during critical study)
    // don't keep multiplying the resolution.
    if (bumped_bs_resolution)
        boson_star.n_gridpoints = bs_n_gridpoints_orig;

    grid_offset = 0.5 * (double)cell_centered;
    BSSNSlice::origin_mirror_shift = cell_centered ? 1 : 0;

    BSSNSlice::one_sided_outer = constraint_preserving_BCs;
    deriv_cache_slice = nullptr;

    if (!run_quietly)
        cout << "Read BS data" << endl;

    // Optionally add a real scalar field Gaussian before evolution
    //if (add_real_field && start_time == 0)
    //    add_real_gaussian(slices[0]);

    //solve Hamiltonian constraint directly in isotropic coords for initial slice if requested
    if (run_spacetime_solver && start_time == 0)
        solve_initial_ham(run_quietly);

    //cut off outermost 2 gripoints, where the christoffel symbols will be generally polluted by garbage due to not having data to take derivatives there. Might be better to just extrapolate long term.
    n_gridpoints -= 2 ;

    R *= (n_gridpoints - 1.0) / (n_gridpoints + 1.0); //also need to rescale R to avoid stretching solution
    slices[0].R = R;

    int n_old = n_gridpoints;
    n_gridpoints = round(cutoff_frac * n_gridpoints); //shrink domain by cutoff_frac, ideally to remove detritus in H
    R = (R * (n_gridpoints - 1.0)) / (n_old - 1.0); //rescale R again to keep dr the same

    // Keep composite state vector in sync with the trimmed grid size
    slices[0].states2.resize(n_gridpoints);
    slices[0].R = R;

    if (use_CCZ4)
        slices[0].theta.resize(n_gridpoints);

    if (start_time > 0)
        slices[0].read_checkpoint(start_time, n_gridpoints);

    //NOTE: n_gridpoints and R should not change after this point!!!

    //the grid has been trimmed, so make the refinement layout consistent with its final size
    finalize_refinement();

    // MPI domain decomposition is computed now, once the grid size is final
    // (initialize() has already trimmed n_gridpoints twice above).
    init_domain_decomposition();

    //resize all auxiliary/diagnostic arrays as appropriate
    resize_temp_arrays();

    //compute auxiliary/diagnostic quantities on initial slice
    current_slice_ptr = &slices[0];

    if (make_tangherlini)
        M = slices[0].make_tangherlini(1., min_chi, D);

    if (make_tangherlini || slices[0].states2[0].bssn.chi < 10 * min_chi)
        slices[0].has_BH = 1;
    else slices[0].has_BH = 0;

    compute_auxiliary_quantities(current_slice_ptr);
    rho0_init = make_tangherlini ? 1. : rho[0];
    compute_diagnostics(current_slice_ptr);
    dtK_L2 = 0;

    double mass0 = slice_mass(&slices[0]);
    double charge0 = slice_charge(&slices[0]);
    M = mass0;
    M_bc = slice_ms_mass(&slices[0], ms_bc_radius_frac); // mass used in the Sommerfeld asymptotics

    if (mpi_rank == 0)
    {
        std::ofstream dynamical_file{"dynamical_constants.dat"};
        dynamical_file << mass0 << "    " << charge0 << "    " << mass0 - mu * charge0;

        if (!run_quietly)
            cout << "Dynamical N =" << charge0 << ", M = "  << mass0 << " and binding energy is " << mass0 - mu * charge0  << endl;
    }
}


//Time evolution! Must have initialized using a constructed BosonStar first.
void Spacetime::evolve()
{
    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    if (critical_study) ricci_4_ctr_max = 0.0;

    //Berger-Oliger subcycling: the loop below then advances by the coarsest level's step and every level L takes
    //2^(n_levels - L) sub-steps of 2^(L-1) dt_fine inside bo_advance_level
    fill_taper_weights();

    const bool use_bo = adaptive_timestep && !refinement_points.empty();
    if (use_bo)
    {
        dt_fine = dt;
        dt *= (1 << refinement_points.size());
        bo_setup();
    }

    if (!run_quietly && mpi_rank == 0)
        cout << "dr = " << dr << ", dt = " << dt <<  "   " << endl;

    // MPI: initial slice was populated identically on all ranks (each rank ran
    // the BosonStar ODE), so slices[0].states2 has valid data at every index and
    // no exchange_ghosts is strictly needed at t=0. Call it anyway as a
    // low-cost safety belt in case initial data ever diverges across ranks.
    exchange_ghosts(slices[0]);

    // AH search is MPI-aware: each rank evaluates ah_indicator on its owned range
    // (using ghost data at the owned edges) and MPI_Allgathers into a full array.

    int num_timesteps = ceil(stop_time / dt);
    int last_checkpoint_time = 0;
    int last_dump_time = 0;
    double super_detect_time = -1.0; //time at which supercriticality was first detected; -1 means not yet detected
    super_ah_radius = 0.0;
    super_M = 0.0;
    super_ah_mass = 0.0;
    super_ms_mass = 0.0;

    //write constraint norms at each timestep to file (rank 0 only)
    std::ofstream constraints_file;
    if (mpi_rank == 0)
    {
        constraints_file.open("constraint_norms.dat");
        if (!constraints_file)
        {
            std::cerr << "constraint_norms.dat could not be opened for writing!\n";
            bss_mpi::abort_all(1);
        }

        constraints_file << "#time"<< "   " << "Ham_L2"
                << "   " << "Mom_L2"<<  "   " << "chi_central" << "   "
                << "A_central" << "   "  << "AH radius" << "   " <<  "M_ADM" << "   "  << "Noether charge"
                << "   " << "Central dK/dt" << "   " << "freq_central" << "   " << "lapse_central" <<
                "   " << "complex field energy" << "   " << "real field energy" << "   " << "central 4D ricci scalar" << "   " << "phi_re_central" << "   " << "AH mass" << "   " << "MS mass" << endl;
    }

    // A0 is the initial central field amplitude, used later for migration
    // detection and the subcritical check. Only rank 0 has states2[0] valid
    // in the general case; broadcast so all ranks share the same value.
    double A0 = 0.0;
    if (owns_origin())
        A0 = sqrt (pow(slices[0].states2[0].csf.phi_re,2) + pow(slices[0].states2[0].csf.phi_im,2));
    bss_mpi::bcast_doubles(&A0, 1);

    if (!run_quietly && mpi_rank == 0)
        cout <<" \n Will evolve with up to " << num_timesteps << " time steps \n" << endl;
        
    //s_i are returned RHS's, t represents temporary RHS + current_slice quantities that must be stored so derivatives can be accessed
    BSSNSlice s1, s2, s3, s4, t1, t2, t3, update_slice;

    // Pure-evolution timing: measure wall time inside the main stepping loop
    // only. All ranks have just completed a collective (bcast of A0) so they
    // enter the loop in sync; final cout is rank-0-only.
    bss_mpi::barrier();
    const auto evol_t0 = std::chrono::steady_clock::now();
    int steps_executed = 0;

    for (int time_step = 0; time_step < num_timesteps; time_step++)
    {
        ++steps_executed;
        double t = start_time + time_step * dt;

        const int skipped_levels = 0; //(the former frozen-RHS level skipping is superseded by the Berger-Oliger path below)

        //fill out array until we've reached maximum number of stored slices, then update last element + rotate at end.
        int n = (time_step > max_stored_slices - 2) ? (max_stored_slices - 2) : time_step;
        double phase_ctr = 0.0, phase_last = 0.0;
        if (time_step > 0) phase_last = phase_ctr;
        // MPI: states2[0] is only valid on rank 0. Compute central-point diagnostics
        // on rank 0 and broadcast the few scalars other ranks need for control flow.
        if (owns_origin())
            phase_ctr = std::arg( std::complex<double>(slices[n].states2[0].csf.phi_re, slices[n].states2[0].csf.phi_im));

        double omega_approx = (phase_ctr - phase_last) / dt;

        if (omega_approx < -M_PI / dt) omega_approx += 2 * M_PI / dt; //corrects jumps due to branch cuts

        if (n > 0) compute_dtK(n);
        M = slice_mass(current_slice_ptr); // volume-integral mass, output as M_ADM
        {
            const double ms_now = slice_ms_mass(current_slice_ptr, ms_bc_radius_frac); // collective
            if (std::isfinite(ms_now)) M_bc = ms_now;           // keep last good value if the slice has gone bad
        }

        // Central-point scalars needed on every rank for control flow.
        // Layout: [0]=A_ctr, [1]=has_BH_now, [2]=chi_ctr, [3]=alpha_ctr,
        //         [4]=phi_re_ctr, [5]=phase_ctr, [6]=test_ctr.
        // (Rank 0 also uses [2..5] to write constraint_norms.dat; the bundle keeps it in one bcast.)
        double central[7] = {0., 0., 0., 0., 0., 0., 0.};
        if (owns_origin())
        {
            const double Ap0 = sqrt (pow(slices[n].states2[0].csf.phi_re,2) + pow(slices[n].states2[0].csf.phi_im,2));
            const double Ap1 = sqrt (pow(slices[n].states2[1].csf.phi_re,2) + pow(slices[n].states2[1].csf.phi_im,2));
            central[0] = (cell_centered) ? (Ap0 - 0.5 * dr * (Ap1 - Ap0)) : Ap0;
            central[1] = (make_tangherlini || slices[n].states2[0].bssn.chi < 10 * min_chi) ? 1.0 : 0.0;
            central[2] = slices[n].states2[0].bssn.chi;
            central[3] = slices[n].states2[0].bssn.alpha;
            central[4] = slices[n].states2[0].csf.phi_re;
            central[5] = phase_ctr;
            central[6] = test_ctr;
        }
        bss_mpi::bcast_doubles(central, 7);
        const double A_ctr        = central[0];
        const bool   has_BH_now   = (central[1] != 0.0);
        const double chi_central  = central[2];
        const double alpha_central= central[3];
        const double phi_re_ctr   = central[4];
        phase_ctr                 = central[5];
        test_ctr                  = central[6];

        double ricci_4 = (critical_study ? ricci_4_grid_max : 0.0);

        if (critical_study && ricci_4 > ricci_4_ctr_max && t > 10)
            ricci_4_ctr_max = ricci_4;

        if (time_step % write_CN_interval == 0) //write time-dependent diagnostics to constraint_norms.dat
        {
            compute_scalar_energies(current_slice_ptr);
            const double charge_now = slice_charge(current_slice_ptr);
            const double ah_mass_now = get_ah_mass(); // 0 on non-owners since ah_radius<0 when AH search is off
            const double ms_mass_now = slice_ms_mass(current_slice_ptr, ms_radius_frac); // collective
            if (mpi_rank == 0)
            {
                constraints_file << std::setprecision (10) << start_time + dt * time_step << "   " << Ham_L2
                << "   " << Mom_L2 <<  "   " << chi_central << "   "
                << A_ctr << "   "  << ah_radius << "   " <<  M << "   "  << charge_now
                << "   " << dtK_L2 << "   " << omega_approx << "   " << alpha_central <<
                "   " << E_phi << "   " << E_psi << "   " << ricci_4 << "   " << phi_re_ctr << "   " << ah_mass_now << "   " << ms_mass_now << endl;
            }
        }

        slices[n + 1].states2.resize(n_gridpoints);
        slices[n + 1].R = R;
        slices[n + 1].refinement_points = refinement_points;
        slices[n + 1].active_points = active_points;

        if (has_BH_now)
            slices[n + 1].has_BH = 1;

        if (use_CCZ4)
            slices[n + 1].theta.resize(n_gridpoints);

        //evaluate intermediate RK4 quantities. exchange_ghosts refreshes each
        //temporary's ghost slots with the corresponding time-shifted values
        //from neighbors before the next RHS reads them across stencils.
        current_slice_ptr = &slices[n];

        if (use_bo)
        {
            //diagnostics on the synchronised slice, then advance all levels in place on a copy
            if (time_step % write_CN_interval == 0)
            {
                compute_auxiliary_quantities(current_slice_ptr, 0);
                compute_diagnostics(current_slice_ptr);
            }
            slices[n + 1].states2 = slices[n].states2;
            if (use_CCZ4) slices[n + 1].theta = slices[n].theta;
            bo_advance_level((int)bo_levels.size(), t, dt, slices[n + 1]); //sets current_slice_ptr = &slices[n + 1]
            current_slice_ptr = &slices[n + 1];
        }
        else
        {
        s1 = slice_rhs(current_slice_ptr, skipped_levels);
        slice_axpy(t1, slices[n], 0.5 * dt, s1);
        fill_inactive_points(t1);
        exchange_ghosts(t1);

        if (time_step % write_CN_interval == 0)
            compute_diagnostics(current_slice_ptr);

        current_slice_ptr = &t1;
        s2 = slice_rhs(current_slice_ptr, skipped_levels);
        slice_axpy(t2, slices[n], 0.5 * dt, s2);
        fill_inactive_points(t2);
        exchange_ghosts(t2);

        current_slice_ptr = &t2;
        s3 = slice_rhs(current_slice_ptr, skipped_levels);
        slice_axpy(t3, slices[n], dt, s3);
        fill_inactive_points(t3);
        exchange_ghosts(t3);

        current_slice_ptr = &t3;
        s4 = slice_rhs(current_slice_ptr, skipped_levels);

        //update slice
        get_update_slice(s1, s2, s3, s4, update_slice, dt, skipped_levels);
        slice_axpy(slices[n + 1], slices[n], 1.0, update_slice);

        //enforce that A is traceless
        current_slice_ptr = &slices[n + 1];
        make_A_traceless(current_slice_ptr);

        //enforce minimum chi on owned points (ghosts get refreshed by exchange_ghosts below)
        for (int j = owned_begin(); j < owned_end(); j++)
        {
            State& st = slices[n + 1].states2[j];
            if (st.bssn.chi < min_chi) st.bssn.chi = min_chi;
        }
        }

        //refresh inactive (refinement) points from the updated active points
        fill_inactive_points(slices[n + 1]);

        // Refresh ghost zones of the freshly-updated slice so next step's RHS
        // sees consistent-time values across the domain boundary.
        exchange_ghosts(slices[n + 1]);

        if (do_ah_search && time_step % AH_search_steps == 0)
        {
            double ah_radius_old = ah_radius;
            ah_radius = apparent_horizon_search(current_slice_ptr);
            if ( ah_radius_old < 0. && ah_radius > 0. && !critical_study && mpi_rank == 0)
                cout << "Apparent horizon found at radius " << ah_radius << " at time " << t + dt << endl;
        }

        if (time_step % write_interval == 0)
        {
            // Both writes are collective (they gather from every rank; rank 0 writes).
            write_current_slice();
            write_diagnostics();
        }

        //write checkpoint files (collective: gather to rank 0, rank 0 writes)
        if ((int)std::floor(t) % checkpoint_time == 0 && (int)std::floor(t) > last_checkpoint_time)
        {
            write_checkpoint(*current_slice_ptr, "checkpoint" + std::to_string((int)std::floor(t)) + ".dat");
            if (mpi_rank == 0)
                cout << "Wrote checkpoint at t = " << (int)std::floor(t) << endl;
            last_checkpoint_time = (int)std::floor(t);
        }

        //write dump files at multiples of dump_time
        if (dump_time > 0)
        {
            int current_dump = (int)(std::floor(t) / dump_time) * dump_time;
            if (current_dump > last_dump_time)
            {
                std::string tag = std::to_string(current_dump);
                write_current_slice("SliceData_" + tag + ".dat");
                write_diagnostics("Diagnostics_" + tag + ".dat");
                last_dump_time = current_dump;
            }
        }

        //cycles slice array back by one so that last entry can be overwritten
        if (time_step >= max_stored_slices - 2)
            rotate(slices.begin(), slices.begin() + 1, slices.end());

        if ((time_step + 1) % 10 == 0 && !run_quietly && mpi_rank == 0)
            cout << "Time step " << time_step + 1 << " complete! t = " << t << endl;

        // Collective NaN check: any rank sees a NaN in its owned chi -> abort together.
        bool local_nan = false;
        for (int j = owned_begin(); j < owned_end() && !local_nan; j++)
            if (std::isnan(slices[n + 1].states2[j].bssn.chi)) local_nan = true;
        if (bss_mpi::any(local_nan))
        {
            if (mpi_rank == 0) cout << "chi became nan on step " << time_step << endl;
            if (critical_study && (super_detect_time >= 0.0 || ah_radius > 0.))
            {
                if (mpi_rank == 0) cout << "Horizon was detected before nan appeared; recording run as supercritical." << endl;
                if (super_detect_time < 0.0)
                {
                    critical_state = 1;
                    critical_time = t;
                    super_detect_time = t;
                    if (ah_radius > 0.) super_ah_radius = ah_radius;
                    super_ah_mass = get_ah_mass();
                    super_M = M;
                    const double ms_now = slice_ms_mass(current_slice_ptr, ms_radius_frac); // slice may already hold nans here
                    if (std::isfinite(ms_now)) super_ms_mass = ms_now;
                }
                break;
            }
            bss_mpi::abort_all(1);
        }

        //if stop_on_migrate enabled, exit when central amplitude changes by more than 10%
        if (stop_on_migrate && abs(A0 - test_ctr) > 0.1 * A0)
        {
            if (mpi_rank == 0)
                cout << "Migration occurred on time step  " << time_step <<" at time t = " << time_step * dt << endl;
            bss_mpi::abort_all(1);
        }

        if (critical_study && super_detect_time < 0.0 && (ah_radius > 0. || A_ctr > A_thresh)) //use AH formation or central amplitude as supercritical indicator
        //(critical_study && slices[n + 1].states2[0].bssn.alpha < lapse_thresh)
        {
            if (mpi_rank == 0) cout << "Supercritical at time " << t << endl;
            critical_state = 1;
            critical_time = t;
            super_detect_time = t;
            if (ah_radius > 0.) { super_ah_radius = ah_radius; super_ah_mass = get_ah_mass(); }
            super_M = M;
            super_ms_mass = slice_ms_mass(current_slice_ptr, ms_radius_frac);
        }

        //if past_super_time > 0, keep the last known good AH radius / mass while extending, then break when the window elapses
        if (critical_study && super_detect_time >= 0.0)
        {
            if (ah_radius > 0.) { super_ah_radius = ah_radius; super_ah_mass = get_ah_mass(); }
            super_M = M;
            super_ms_mass = slice_ms_mass(current_slice_ptr, ms_radius_frac);
            if (past_super_time <= 0.0 || t >= super_detect_time + past_super_time)
                break;
        }
        //wait some time to avoid initial transients, then use net central amplitude decrease as subcritical indicator
        if (critical_study && super_detect_time < 0.0 && A_ctr < A0 && t > sub_min_time)
        {
            if (mpi_rank == 0) cout << "Subcritical at time " << t << endl;
            critical_state = 0;
            critical_time = t;
            break;
        }

        // Additional subcritical criterion: if no AH by subcritical_time, assume subcritical
        if (critical_study && super_detect_time < 0.0 && t >= subcritical_time)
        {
            if (ah_radius < 0.0)
            {
                if (mpi_rank == 0) cout << "Subcritical by time criterion at time " << t << endl;
                critical_state = 0;
                critical_time = t;
                break;
            }
        }
    }

    bss_mpi::barrier();
    const auto evol_t1 = std::chrono::steady_clock::now();
    const double evol_elapsed = std::chrono::duration<double>(evol_t1 - evol_t0).count();
    if (mpi_rank == 0)
    {
        const double us_per_step = (steps_executed > 0)
            ? (evol_elapsed * 1.0e6 / steps_executed) : 0.0;
        cout << "[timing] evolve_loop: " << evol_elapsed << " s over "
             << steps_executed << " timesteps ("
             << us_per_step << " us/step) on " << mpi_size << " ranks" << endl;
    }
}

// Tune a parameter by bisection to bracket the critical threshold.
// hi_guess must produce supercritical (critical_state = 1), lo_guess subcritical (critical_state = 0).
// Will repeatedly initialize and evolve, restoring tuning_param after initialize in (likely) case it's a Spacetime member overwritten by params.
void Spacetime::tune_to_critical(double& tuning_param, double hi_guess, double lo_guess, BosonStar* bs)
{
    std::ofstream subcritical_file{"subcritical.dat"};
    std::ofstream supercritical_file{"supercritical.dat"};

    subcritical_file << "#tuning_param   critical_time     ricci_4_max" << endl;
    supercritical_file << "#tuning_param   critical_time   ah_radius_final   M_final   AH_mass_final   MS_mass_final" << endl;

    if (bs == nullptr)
    {
        cerr << "ERROR: tune_to_critical called with null BosonStar pointer" << endl;
        exit(1);
    }

    // Ensure we record critical outcomes during evolve()
    critical_study = true;

    auto run_with_param = [&](double val) -> int {
        // Assign target value, initialize (which may overwrite), then restore
        tuning_param = val;
        double tmp = tuning_param;
        cout << "Running with tuning_param = " << std::setprecision(16) << tuning_param << endl;

        if (critical_gaussian_start)
        {
            // Gaussian scalar cloud initial data path: perturb_amp is the tuning param,
            // real scalar field is disabled.
            bs->add_real_field = false;
            bs->clear_BS();
            bs->omega = bs->enforced_freq;
            bs->add_perturbation(bs->perturb_amp, bs->perturb_spread, 0.);

            //if (!bs->fill_given_A(bs->omega, 0))
                bs->default_metric_vars();

            bs->fill_isotropic_arrays();
        }
        else
        {
            if (!bs->solve()) //re-solve BS now that real perturbation is added at initial phase.
            {
                cout << "solve failed; assuming supercritical for tuning_param = " << std::setprecision(16) << tuning_param << endl;
                return 1;
            }
            //cout << "Real amplitude: " << bs->real_amp << endl;

            // Ensure isotropic arrays (and rsf_iso) reflect the newly solved BS state so initialize/read_BS_data sees updated profiles
            bs->fill_isotropic_arrays();
        }

        // Skip re-reading parameters so tuned value persists into initial data seeding
        initialize(*bs, /*skip_read=*/true);
        if (run_spacetime_solver)
            solve_initial_ham(run_quietly);
        tuning_param = tmp;
        evolve();
        if (critical_state != 0 && critical_state != 1)
        {
            cerr << "ERROR: evolve() ended with critical_state = " << critical_state
                 << ", expected 0 (subcritical) or 1 (supercritical)." << endl;
            exit(1);
        }
        return critical_state;
    };

    // Verify high guess is supercritical
    int hi_state = run_with_param(hi_guess);
    if (hi_state != 1)
    {
        cerr << "ERROR: hi_guess = " << std::setprecision(16) << hi_guess
             << " did not produce supercritical outcome (critical_state = 1). Got critical_state = "
             << hi_state << "." << endl;
        exit(1);
    }
    supercritical_file << std::setprecision(16)
                               << hi_guess<< "   " << critical_time << "   " << super_ah_radius << "   " << super_M << "   " << super_ah_mass << "   " << super_ms_mass << endl;

    // Verify low guess is subcritical
    int lo_state = run_with_param(lo_guess);
    if (lo_state != 0)
    {
        cerr << "ERROR: lo_guess = " << std::setprecision(16) << lo_guess
             << " did not produce subcritical outcome (critical_state = 0). Got critical_state = "
             << lo_state << "." << endl;
        exit(1);
    }

    subcritical_file << std::setprecision(16)
                             << lo_guess << "   " << critical_time << "   " << ricci_4_ctr_max << endl;
    double hi = hi_guess;
    double lo = lo_guess;

    // Bisection until interval smaller than tolerance
    while (fabs(hi - lo) / lo > critical_eps)
    {
        double mid = 0.5 * (hi + lo);
        int state = run_with_param(mid);
        if (state == 1)
        {
            hi = mid;   // supercritical at mid -> reduce upper bound
            supercritical_file << std::setprecision(16)
                               << mid << "   " << critical_time << "   " << super_ah_radius << "   " << super_M << "   " << super_ah_mass << "   " << super_ms_mass << endl;
        }
        else if (state == 0)
        {
            lo = mid;   // subcritical at mid -> raise lower bound
            subcritical_file << std::setprecision(16)
                             << mid << "   " << critical_time << "   " << ricci_4_ctr_max << endl;
        }
        else
        {
            cerr << "ERROR: Unexpected critical_state during bisection: " << state << endl;
            exit(1);
        }
    }

    tuning_param = 0.5 * (hi + lo);
    cout << std::setprecision(16)
         << "Critical tuning complete. tuning_param = " << tuning_param
         << ", bracket [" << lo << ", " << hi << "] with tol = " << critical_eps << endl;


    if (fill_subcritical)
    {
        double step = critical_eps;
        // Fill in subcritical data points for plotting
        while (fabs(step) < fabs(0.01 * tuning_param))
        {
            subcritical_file <<"#Critical tuning finished here!" << endl;
            double sub_val = tuning_param - step;
            run_with_param(sub_val);
            subcritical_file << std::setprecision(16)
                             << sub_val << "   " << critical_time << "   " << ricci_4_ctr_max << endl;
            step *= 2.0;
        }
    }

    exit(0);
}

//Expansion for AH search
double Spacetime::ah_indicator_at(const BSSNSlice* slice_ptr, int j) 
{
    // Access BSSN variables to avoid unused-variable warnings and as a cue for future implementation.
    const BSSNState& s = slice_ptr->states2[j].bssn;
    double z = (j + grid_offset) * dr;
    if (z == 0.0)
        z = 1e-12; // avoid division by zero at origin

    // outgoing null expansion / (D-2): sqrt(chi/h_zz) d_z ln R_areal - K^theta_theta, with
    // d_z ln R_areal = 1/z + d_z h_ww / (2 h_ww) - d_z chi / (2 chi) (the D = 4 form had the prefactor 2 hard-coded)
    double ah_indicator = sqrt(s.chi / s.h_zz) *
                          (1. / z + 0.5 * (d_z(v_h_ww, j) - s.h_ww * d_z(v_chi, j) / s.chi) / s.h_ww)
                          - (s.A_ww + s.K * s.h_ww / (D - 1.)) / s.h_ww;

    return ah_indicator;
}

// Scan slice for sign changes of the AH indicator and return outermost root radius. Returns -1.0 if no AH is found.
// MPI: each rank evaluates the indicator on its owned range (using ghost data
// for derivative stencils at the owned edges), Allgatherv assembles a full
// indicator array on every rank, and every rank walks the same array so the
// returned radius is identical across ranks with no follow-up broadcast needed.
double Spacetime::apparent_horizon_search(const BSSNSlice* slice_ptr)
{
    if (slice_ptr == nullptr || slice_ptr->states2.empty())
        return -1.0;

    const int N = static_cast<int>(slice_ptr->states2.size());
    const double dr_local = slice_ptr->R / (static_cast<double>(N) - 1.0);
    const double offset = 0.5 * static_cast<double>(cell_centered); // match grid_offset convention

    // Build a full-grid indicator array. In serial we compute it inline via the
    // usual walk; under MPI each rank computes its owned slice and Allgathers.
    std::vector<double> full;
    if (mpi_size > 1)
    {
        std::vector<double> local(local_N);
        for (int k = 0; k < local_N; ++k)
            local[k] = ah_indicator_at(slice_ptr, owned_begin() + k);
        allgather_owned(local, full, 1);
    }

    auto f_at = [&](int j) -> double {
        return (mpi_size > 1) ? full[j] : ah_indicator_at(slice_ptr, j);
    };

    double outermost_r = -1.0;

    double f_prev = f_at(0);
    double r_prev = (0 + offset) * dr_local;

    //only active points carry evolved data, so look for sign changes between consecutive active points
    for (int j = 1; j < N; ++j)
    {
        if (j < (int)active_points.size() && !active_points[j])
            continue;

        double f_curr = f_at(j);
        double r_curr = (j + offset) * dr_local;

        // Exact zero at a node counts as a root
        if (f_curr == 0.0)
            outermost_r = r_curr;

        // Look for sign changes indicating a crossing
        if ((f_prev > 0.0 && f_curr < 0.0) || (f_prev < 0.0 && f_curr > 0.0))
        {
            // Linear interpolation to estimate root location within [r_prev, r_curr]
            double r_zero = r_prev - f_prev * (r_curr - r_prev) / (f_curr - f_prev);
            if (r_zero > outermost_r)
                outermost_r = r_zero;
        }

        f_prev = f_curr;
        r_prev = r_curr;
    }
    return outermost_r;
}

#endif /*EVOLUTIONVARIABLES*/