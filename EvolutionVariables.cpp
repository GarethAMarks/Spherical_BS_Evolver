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

using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::vector;
using std::min;
using std::max;
using std::rotate;
// Spacetime constructor: initialize reusable CSF model with this context
Spacetime::Spacetime() : csf_model(this), rsf_model(this), cvf_model(this), spf_model(this) {

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

// Compute (counts, displs) for MPI_Gatherv/Allgatherv given a per-point double
// stride. Uses the same base/rem split as init_domain_decomposition so the
// layout matches what each rank actually owns.
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
    // 28 doubles per point (theta as 0 when not CCZ4 — but we handle theta
    // as a separate gather below when CCZ4 is on); then 8 Proca and 4 Dirac variables.
    constexpr int NCOL = 28;
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
        const CVF& v = slice.states2[j].cvf;
        p[16] = v.X_re;    p[17] = v.X_im;    p[18] = v.E_re; p[19] = v.E_im;
        p[20] = v.Xphi_re; p[21] = v.Xphi_im; p[22] = v.Z_re; p[23] = v.Z_im;
        const CSpF& sp = slice.states2[j].spf;
        p[24] = sp.u_re; p[25] = sp.u_im; p[26] = sp.v_re; p[27] = sp.v_im;
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
        scratch.states2[j].cvf = CVF{p[16], p[17], p[18], p[19], p[20], p[21], p[22], p[23]};
        scratch.states2[j].spf = CSpF{p[24], p[25], p[26], p[27]};
    }
    scratch.write_slice(file_name, proca, dirac);
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
//   [BSSNState (9)] [CSF (4)] [RSF (2)] [CVF (8)] [CSpF (4)]
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
            p[15] = s.cvf.X_re;    p[16] = s.cvf.X_im;    p[17] = s.cvf.E_re; p[18] = s.cvf.E_im;
            p[19] = s.cvf.Xphi_re; p[20] = s.cvf.Xphi_im; p[21] = s.cvf.Z_re; p[22] = s.cvf.Z_im;
            p[23] = s.spf.u_re; p[24] = s.spf.u_im; p[25] = s.spf.v_re; p[26] = s.spf.v_im;
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
            s.cvf = CVF{p[15], p[16], p[17], p[18], p[19], p[20], p[21], p[22]};
            s.spf = CSpF{p[23], p[24], p[25], p[26]};
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

// Composite State elementwise operations (applies to all of csf, bssn, rsf, cvf and spf)
State operator+(const State& a, const State& b) {
    return State{a.csf + b.csf, a.bssn + b.bssn, a.rsf + b.rsf, a.cvf + b.cvf, a.spf + b.spf};
}

State operator-(const State& a, const State& b) {
    return State{a.csf - b.csf, a.bssn - b.bssn, a.rsf - b.rsf, a.cvf - b.cvf, a.spf - b.spf};
}

State operator*(const State& a, const State& b) {
    return State{a.csf * b.csf, a.bssn * b.bssn, a.rsf * b.rsf, a.cvf * b.cvf, a.spf * b.spf};
}

State operator*(double c, const State& a) {
    return State{c * a.csf, c * a.bssn, c * a.rsf, c * a.cvf, c * a.spf};
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
double BSSNSlice::d_z(bssn_var var, int index, int order = 1 )
{
    // Use composite state vector size for grid length to ensure compatibility with RK temp slices
    int n_gridpoints = states2.size();
    // Use cell-centred grid spacing to match read_BS_data (innermost point at z = 0.5*dr)
    double dr = R / ((double)n_gridpoints - 1.0);

    // stencil spacing 2^(ref_level-1), where the refinement level starts at 1 for no refinement; subcycling may override it
    int res_fac;
    if (forced_res_fac && index >= 0 && index < (int)forced_res_fac->size() && (*forced_res_fac)[index] > 0)
        res_fac = (*forced_res_fac)[index];
    else
        res_fac = 1 << (get_refinement_level(index, refinement_points) - 1);

    //check index is valid and error out if not
   if (index < 0 || index >= n_gridpoints )
    {
        cerr << "ERROR: invalid index requested in radial derivative" << endl;
        exit(1);
    }

    //set of indices using which derivative will be evaluated
    vector<int> J{index - 2 * res_fac, index - res_fac, index, index + res_fac, index + 2 * res_fac};

    //enforce BC at z = 0 by using symmetry
    if (index <= 1)
        J[0] = -J[0];
    if (index == 0)
        J[1] = -J[1];


    //At outer edge just use central value for rightmost two; this is cut off in initialize()/ overwritten with BCs anyway
    if (index >= n_gridpoints - res_fac)
        J[3] = J[2];
    if (index >= n_gridpoints - 2 * res_fac)
        J[4] = J[3];

    //set of values used to compute derivative
    vector<double> F{0., 0., 0., 0., 0.};

    //Fill out five values; note that at outer boundary this gives unreliable results
    if (index < n_gridpoints)
    {
        for (int j  = 0; j < 5; j++)
        {
            switch (var)
            {
                case v_chi:
                    F[j] = states2[J[j]].bssn.chi;
                    break;
                case v_h_zz:
                    F[j] = states2[J[j]].bssn.h_zz;
                    break;
                case v_h_ww:
                    F[j] = states2[J[j]].bssn.h_ww;
                    break;
                case v_A_zz:
                    F[j] = states2[J[j]].bssn.A_zz;
                    break;
                case v_A_ww:
                    F[j] = states2[J[j]].bssn.A_ww;
                    break;
                case v_K:
                    F[j] = states2[J[j]].bssn.K;
                    break;
                case v_c_chris_Z:
                    F[j] = states2[J[j]].bssn.c_chris_Z;
                    break;
                case v_phi_re:
                    F[j] = states2[J[j]].csf.phi_re;
                    break;
                case v_phi_im:
                    F[j] = states2[J[j]].csf.phi_im;
                    break;
                case v_K_phi_re:
                    F[j] = states2[J[j]].csf.K_phi_re;
                    break;
                case v_K_phi_im:
                    F[j] = states2[J[j]].csf.K_phi_im;
                    break;
                case v_alpha:
                    F[j] = states2[J[j]].bssn.alpha;
                    break;
                case v_beta:
                    F[j] = states2[J[j]].bssn.beta;
                    break;
                case v_theta:
                    F[j] = theta[J[j]];
                    break;
                case v_psi:
                    F[j] = states2[J[j]].rsf.psi;
                    break;
                case v_K_psi:
                    F[j] = states2[J[j]].rsf.K_psi;
                    break;
                case v_X_re:
                    F[j] = states2[J[j]].cvf.X_re;
                    break;
                case v_X_im:
                    F[j] = states2[J[j]].cvf.X_im;
                    break;
                case v_E_re:
                    F[j] = states2[J[j]].cvf.E_re;
                    break;
                case v_E_im:
                    F[j] = states2[J[j]].cvf.E_im;
                    break;
                case v_Xphi_re:
                    F[j] = states2[J[j]].cvf.Xphi_re;
                    break;
                case v_Xphi_im:
                    F[j] = states2[J[j]].cvf.Xphi_im;
                    break;
                case v_Z_re:
                    F[j] = states2[J[j]].cvf.Z_re;
                    break;
                case v_Z_im:
                    F[j] = states2[J[j]].cvf.Z_im;
                    break;
                case v_u_re:
                    F[j] = states2[J[j]].spf.u_re;
                    break;
                case v_u_im:
                    F[j] = states2[J[j]].spf.u_im;
                    break;
                case v_v_re:
                    F[j] = states2[J[j]].spf.v_re;
                    break;
                case v_v_im:
                    F[j] = states2[J[j]].spf.v_im;
                    break;
                default:
                    cerr << "ERROR: invalid variable requested for differentiation" << endl;
                    exit(0);
            }
        }

        bool parity_is_odd = ((var == v_c_chris_Z || var == v_beta || var == v_X_re || var == v_X_im || var == v_E_re || var == v_E_im
                               || var == v_v_re || var == v_v_im)/*|| (has_BH && (var == v_chi || var == v_alpha))*/);
       //if (var == v_alpha) cout << parity_is_odd<< endl;

        //account for odd parity of contracted christoffel symbols, beta, the radial Proca components and the Dirac v across z = 0
            if (index <= 1 && parity_is_odd )
                F[0] = -F[0];
            if (index == 0 && parity_is_odd)
                F[1] = -F[1];
    }
    return fivePointDeriv(res_fac * dr, order, F[0],F[1],F[2],F[3],F[4]);
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
    states2[j].cvf = CVF{};
    states2[j].spf = CSpF{};

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
    int res_fac = (rhs_res_fac && (*rhs_res_fac)[j] > 0) ? (*rhs_res_fac)[j] : (1 << (slice_ptr->get_refinement_level(j, refinement_points) - 1));

    // indices for 7-point centered stencil with symmetry fix near the origin
    std::array<int, 7> J = {j - 3 * res_fac, j - 2 * res_fac, j - 1 * res_fac, j, j + 1 * res_fac, j + 2 * res_fac, j + 3 * res_fac};
    if (j < 3) {
        for (int k = 0; k < 7; ++k) J[k] = std::max(-J[k], J[k]);
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
    // Christoffel symbol and shift, the Proca X_z, E^z and the Dirac v) with the correct sign flip for reflected points
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
        if (proca)
        {
            damping_corr.cvf.X_re = odd_corr([](const State& st) { return st.cvf.X_re; });
            damping_corr.cvf.X_im = odd_corr([](const State& st) { return st.cvf.X_im; });
            damping_corr.cvf.E_re = odd_corr([](const State& st) { return st.cvf.E_re; });
            damping_corr.cvf.E_im = odd_corr([](const State& st) { return st.cvf.E_im; });
        }
        if (dirac)
        {
            damping_corr.spf.v_re = odd_corr([](const State& st) { return st.spf.v_re; });
            damping_corr.spf.v_im = odd_corr([](const State& st) { return st.spf.v_im; });
        }
    }

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

    // asymptotic states and their derivatives where relevant
    State asymp_state{}; // matter part zero by default
    asymp_state.bssn = (BSSNState){chi_asymp(R - dr * res_fac), 1., 1., 0., 0., 0., 0., alpha_asymp(R - dr * res_fac), 0.};
    State asymp_deriv{}; // r-derivative of asymptotic expansion (matter zero)
    asymp_deriv.bssn = (BSSNState){d_chi_asymp(R - dr * res_fac), 0., 0., 0., 0., 0., 0., d_alpha_asymp(R - dr * res_fac), 0.};

    //version without spatially varying asymptotics
    if (!spatially_varying_BC)
    {
        asymp_state.bssn = (BSSNState){1., 1., 1., 0., 0., 0., 0., 1., 0.};
        asymp_deriv.bssn = (BSSNState){0., 0., 0., 0., 0., 0., 0., 0., 0.};
    }

    //maybe adjust to account for purported 1/r^2 decay in K!
    State N{};
    N.bssn = (BSSNState){1., 1., 1., 1., 1., 1., 1., 1., 1.};
    N.csf = (CSF){1., 1., 1., 1.};
    N.cvf = CVF{1., 1., 1., 1., 1., 1., 1., 1.};
    N.spf = CSpF{1., 1., 1., 1.};
    
    if (add_real_field)
        N.rsf = (RSF){1., 1.};

    N = pow(R, (4. - D))  * N;

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
    char_speeds.cvf = CVF{1., 1., 1., 1., 1., 1., 1., 1.};
    char_speeds.spf = CSpF{1., 1., 1., 1.};

    if (add_real_field)
        char_speeds.rsf = (RSF){1., 1.};

    rhs.states2[last_active_j - res_fac] = (-1.) * char_speeds * ( p4_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                         + (0.5 * D - 1.) * N * (s4 - asymp_state) /  (dr * (last_active_j - res_fac)) - asymp_deriv);

    //update variable asymptotic states to outermost edge
    if (spatially_varying_BC)
    {asymp_state.bssn.chi = chi_asymp(R);  asymp_state.bssn.alpha = alpha_asymp(R);
    asymp_deriv.bssn.chi = d_chi_asymp(R);  asymp_deriv.bssn.alpha = d_alpha_asymp(R);}

    rhs.states2[last_active_j] = (-1) * char_speeds * ( p5_stencil(dr * res_fac, s1, s2, s3, s4, s5) 
                                 + (0.5 * D - 1.) * N * (s5 - asymp_state) / (dr * last_active_j) - asymp_deriv);

    if (use_CCZ4)
    {
        double& t1 = slice_ptr->theta[last_active_j - 4 * res_fac];
        double& t2 = slice_ptr->theta[last_active_j - 3 * res_fac];
        double& t3 = slice_ptr->theta[last_active_j - 2 * res_fac];
        double& t4 = slice_ptr->theta[last_active_j - 1 * res_fac];
        double& t5 = slice_ptr->theta[last_active_j];

        rhs.theta[last_active_j - res_fac] = -(p4_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                             + (0.5 * D - 1.) * t4 / (dr * (last_active_j - res_fac)));

        rhs.theta[last_active_j] = -(p5_stencil(dr * res_fac, t1, t2, t3, t4, t5) 
                                   + (0.5 * D - 1.) * t5 / (dr * (last_active_j)));
    }
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
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j) * dr;

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

//reads data from a solved ProcaStar into the initial slice (Proca star analogue of read_BS_data).
//At t = 0, A = f dt + i g dr, so with zero shift: X_phi = -f / alpha (real), X_z = i g (imaginary) and
//E_z = F_{zt} / alpha = (d_z f - omega g) / alpha (real), raised with gamma^zz to give the evolved E^z.
void BSSNSlice::read_PS_data(ProcaStar& proca_star, int BS_resolution_factor, bool isotropic, bool cell_centered)
{
    if ((proca_star.n_gridpoints + BS_resolution_factor - 1) % BS_resolution_factor != 0 )
        cout << "WARNING: Incompatible resolution factor used in read_PS_data" << endl;

    int n_gridpoints = (proca_star.n_gridpoints + BS_resolution_factor - 1) / BS_resolution_factor;
    states2.resize(n_gridpoints);

    if (use_CCZ4)
       theta.resize(n_gridpoints);

    R = proca_star.R;
    double dr = R / ((double)n_gridpoints - 1.0);
    double D = proca_star.D;

    int N_ps = proca_star.n_gridpoints;
    double dr_ps = proca_star.R / (double)(max(1, N_ps - 1));

    std::vector<double> X_nodes(N_ps);
    for (int k = 0; k < N_ps; ++k) X_nodes[k] = k * dr_ps;

    if (proca_star.g_array.size() != proca_star.state.size())
        proca_star.fill_g_array();

    //source arrays: conformal factor (psi or X), log lapse, f, and the radial components g, E (odd in r)
    std::vector<double> c_src(N_ps), phi_src(N_ps), f_src(N_ps), g_src(N_ps), E_src(N_ps);
    for (int k = 0; k < N_ps; ++k)
    {
        if (isotropic)
        {
            c_src[k] = proca_star.psi_iso_array[k];  phi_src[k] = proca_star.phi_iso_array[k];
            f_src[k] = proca_star.f_iso_array[k];    g_src[k] = proca_star.g_iso_array[k];   E_src[k] = proca_star.E_iso_array[k];
        }
        else
        {
            c_src[k] = proca_star.state[k].X;  phi_src[k] = proca_star.state[k].phi;
            f_src[k] = proca_star.state[k].f;  g_src[k] = proca_star.g_array[k];   E_src[k] = proca_star.state[k].E;
        }
    }

    //even functions have zero first derivative at the origin, odd ones zero second derivative
    auto make_spline = [&](tk::spline& sp, const std::vector<double>& y, bool odd)
    {
        if (odd) sp.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0);
        else     sp.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        sp.set_points(X_nodes, y, tk::spline::cspline);
    };
    tk::spline s_c, s_phi, s_f, s_g, s_E;
    make_spline(s_c, c_src, 0); make_spline(s_phi, phi_src, 0); make_spline(s_f, f_src, 0);
    make_spline(s_g, g_src, 1); make_spline(s_E, E_src, 1);

    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j + 0.5 * (double)cell_centered) * dr;

        if (isotropic)
        {
            states2[j].bssn.chi = pow(s_c(z), -4.);
            states2[j].bssn.h_zz = 1.;
            states2[j].bssn.h_ww = 1.;
        }
        else
        {
            double X_interp = s_c(z);
            states2[j].bssn.chi = pow(X_interp, -2. / (D - 1.));
            states2[j].bssn.h_zz = pow(X_interp, (2. * D - 4.) / (D - 1.));
            states2[j].bssn.h_ww = pow(X_interp, -2. / (D - 1.));
        }

        //static star: K_ij = 0
        states2[j].bssn.A_zz = 0.;
        states2[j].bssn.A_ww = 0.;
        states2[j].bssn.K = 0.;

        double alpha = exp(s_phi(z));
        states2[j].bssn.alpha = alpha;
        states2[j].bssn.beta = 0.;

        double g_ZZ = states2[j].bssn.chi / states2[j].bssn.h_zz; //physical gamma^zz

        states2[j].csf = CSF{0., 0., 0., 0.};
        states2[j].rsf = RSF{0., 0.};

        CVF& v = states2[j].cvf;
        v.Xphi_re = -s_f(z) / alpha;
        v.Xphi_im = 0.;
        v.X_re = 0.;
        v.X_im = s_g(z);
        v.E_re = g_ZZ * s_E(z) / alpha;
        v.E_im = 0.;
        v.Z_re = 0.;
        v.Z_im = 0.;

        if (use_CCZ4)
            theta[j] = 0.;
    }

    //contracted Christoffel symbols, as in read_BS_data
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j) * dr;
        double h_ZZ = 1 / states2[j].bssn.h_zz;
        double h_WW = 1 / states2[j].bssn.h_ww;

        states2[j].bssn.c_chris_Z = 0.5 * h_ZZ * h_ZZ * d_z(v_h_zz,j) + (D - 2.) * (-0.5 * h_ZZ * h_WW * d_z(v_h_ww,j));

        if (z > 1e-8)
            states2[j].bssn.c_chris_Z += (D - 2.) * (h_WW * (1 - h_ZZ * states2[j].bssn.h_ww) / z);

        if (isotropic)
            states2[j].bssn.c_chris_Z = 0.;
    }
}

//reads data from a solved DiracStar into the initial slice (Dirac star analogue of read_BS_data). At t = 0 the spinor
//amplitudes are real: u = g (even in r), v = f (odd), the physical radial spinor functions, which are orthonormal-frame
//components and so are the same in areal and isotropic coordinates. The star is static, so K_ij = 0 and beta = 0.
void BSSNSlice::read_DS_data(DiracStar& dirac_star, int BS_resolution_factor, bool isotropic, bool cell_centered)
{
    if ((dirac_star.n_gridpoints + BS_resolution_factor - 1) % BS_resolution_factor != 0 )
        cout << "WARNING: Incompatible resolution factor used in read_DS_data" << endl;

    int n_gridpoints = (dirac_star.n_gridpoints + BS_resolution_factor - 1) / BS_resolution_factor;
    states2.resize(n_gridpoints);

    if (use_CCZ4)
       theta.resize(n_gridpoints);

    R = dirac_star.R;
    double dr = R / ((double)n_gridpoints - 1.0);
    double D = dirac_star.D;

    int N_ds = dirac_star.n_gridpoints;
    double dr_ds = dirac_star.R / (double)(max(1, N_ds - 1));

    std::vector<double> X_nodes(N_ds);
    for (int k = 0; k < N_ds; ++k) X_nodes[k] = k * dr_ds;

    //source arrays: conformal factor (psi or X), log lapse, and the physical spinor functions g (even) and f (odd)
    std::vector<double> c_src(N_ds), phi_src(N_ds), g_src(N_ds), f_src(N_ds);
    for (int k = 0; k < N_ds; ++k)
    {
        if (isotropic)
        {
            c_src[k] = dirac_star.psi_iso_array[k];  phi_src[k] = dirac_star.phi_iso_array[k];
            g_src[k] = dirac_star.g_iso_array[k];    f_src[k] = dirac_star.f_iso_array[k];
        }
        else
        {
            const DiracFieldState& st = dirac_star.state[k];
            const double r = k * dr_ds, sa = exp(0.5 * st.phi);
            c_src[k] = st.X;  phi_src[k] = st.phi;
            g_src[k] = (k == 0) ? dirac_star.A_central : st.G / (sa * r);
            f_src[k] = (k == 0) ? 0. : st.F / (sa * r);
        }
    }

    auto make_spline = [&](tk::spline& sp, const std::vector<double>& y, bool odd)
    {
        if (odd) sp.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0);
        else     sp.set_boundary(tk::spline::first_deriv, 0.0, tk::spline::second_deriv, 0.0);
        sp.set_points(X_nodes, y, tk::spline::cspline);
    };
    tk::spline s_c, s_phi, s_g, s_f;
    make_spline(s_c, c_src, 0); make_spline(s_phi, phi_src, 0); make_spline(s_g, g_src, 0); make_spline(s_f, f_src, 1);

    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j + 0.5 * (double)cell_centered) * dr;

        if (isotropic)
        {
            states2[j].bssn.chi = pow(s_c(z), -4.);
            states2[j].bssn.h_zz = 1.;
            states2[j].bssn.h_ww = 1.;
        }
        else
        {
            double X_interp = s_c(z);
            states2[j].bssn.chi = pow(X_interp, -2. / (D - 1.));
            states2[j].bssn.h_zz = pow(X_interp, (2. * D - 4.) / (D - 1.));
            states2[j].bssn.h_ww = pow(X_interp, -2. / (D - 1.));
        }

        states2[j].bssn.A_zz = 0.;
        states2[j].bssn.A_ww = 0.;
        states2[j].bssn.K = 0.;
        states2[j].bssn.alpha = exp(s_phi(z));
        states2[j].bssn.beta = 0.;

        states2[j].csf = CSF{0., 0., 0., 0.};
        states2[j].rsf = RSF{0., 0.};
        states2[j].cvf = CVF{};
        states2[j].spf = CSpF{s_g(z), 0., s_f(z), 0.};

        if (use_CCZ4)
            theta[j] = 0.;
    }

    //contracted Christoffel symbols, as in read_BS_data
    for (int j = 0; j < n_gridpoints; j++)
    {
        double z = (j) * dr;
        double h_ZZ = 1 / states2[j].bssn.h_zz;
        double h_WW = 1 / states2[j].bssn.h_ww;

        states2[j].bssn.c_chris_Z = 0.5 * h_ZZ * h_ZZ * d_z(v_h_zz,j) + (D - 2.) * (-0.5 * h_ZZ * h_WW * d_z(v_h_ww,j));

        if (z > 1e-8)
            states2[j].bssn.c_chris_Z += (D - 2.) * (h_WW * (1 - h_ZZ * states2[j].bssn.h_ww) / z);

        if (isotropic)
            states2[j].bssn.c_chris_Z = 0.;
    }
}

//writes BSSN evolution variables from a particular slice to text file; if write_cvf, the 8 Proca variables are appended,
//and if write_spf, the 4 Dirac variables
void BSSNSlice::write_slice(std::string file_name, bool write_cvf, bool write_spf)
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
        << "   " << states2[j].rsf.psi << "    " << states2[j].rsf.K_psi;

        if (write_cvf)
        {
            const CVF& v = states2[j].cvf;
            data_file << "    " << v.X_re << "    " << v.X_im << "    " << v.E_re << "    " << v.E_im
                      << "    " << v.Xphi_re << "    " << v.Xphi_im << "    " << v.Z_re << "    " << v.Z_im;
        }
        if (write_spf)
        {
            const CSpF& sp = states2[j].spf;
            data_file << "    " << sp.u_re << "    " << sp.u_im << "    " << sp.v_re << "    " << sp.v_im;
        }
        data_file << endl;
    }
}

//reads slice data from checkpoint file of form checkpoint(time).dat w/o parentheses. If read_cvf (read_spf), the file is
//expected to contain all columns written by write_slice with write_cvf (write_spf), including the 8 Proca (4 Dirac) variables.
void BSSNSlice::read_checkpoint(int time, int n_gridpoints, bool read_cvf, bool read_spf)
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

        if (read_spf)
        {
            double t0, psi, K_psi;
            CSpF sp{};
            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta >> t0
                   >> psi >> K_psi >> sp.u_re >> sp.u_im >> sp.v_re >> sp.v_im)
            {
                if (j >= n_gridpoints) break;
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                states2[j].rsf.psi = psi; states2[j].rsf.K_psi = K_psi;
                states2[j].spf = sp;
                if (use_CCZ4) theta[j] = t0;
                j++;
            }
        }
        else if (read_cvf)
        {
            double t0, psi, K_psi;
            CVF v{};
            if(iss >> z >> chi >> h_zz >> h_ww >> A_zz >> A_ww >> K >> c_chris_Z >> phi_re >> phi_im >> K_phi_re >> K_phi_im >> alpha >> beta >> t0
                   >> psi >> K_psi >> v.X_re >> v.X_im >> v.E_re >> v.E_im >> v.Xphi_re >> v.Xphi_im >> v.Z_re >> v.Z_im)
            {
                if (j >= n_gridpoints) break;
                states2[j].bssn.chi = chi; states2[j].bssn.h_zz = h_zz; states2[j].bssn.h_ww = h_ww;
                states2[j].bssn.A_zz = A_zz; states2[j].bssn.A_ww = A_ww; states2[j].bssn.K = K; states2[j].bssn.c_chris_Z = c_chris_Z;
                states2[j].csf.phi_re = phi_re; states2[j].csf.phi_im = phi_im; states2[j].csf.K_phi_re = K_phi_re; states2[j].csf.K_phi_im = K_phi_im;
                states2[j].bssn.alpha = alpha; states2[j].bssn.beta = beta;
                states2[j].rsf.psi = psi; states2[j].rsf.K_psi = K_psi;
                states2[j].cvf = v;
                if (use_CCZ4) theta[j] = t0;
                j++;
            }
        }
        else if (use_CCZ4)
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

    while ( k < n_refinements && j >= refinement_points[k] - pow(2, k))//check this difference -- meant to ensure we can use a stencil at points spaced by 2^(k + 1) safely at j
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
    return current_slice_ptr->d_z(var, index, order);
}

double Spacetime::d_zz(bssn_var var, int index)
{
    return current_slice_ptr->d_zz(var, index);
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

//asymptotic values of chi and alpha and their derivatives. TODO: generalize D > 4
double Spacetime::chi_asymp(double r)
{
    return(isotropic) ? (pow(1. +  0.5 * M / r, -4.)) : pow(1 - 2. * M / r, -1./3.); //use isotropic/areal Schwarzchilld value for chi as appropriate
}

double Spacetime::alpha_asymp(double r)
{
        return sqrt(1 - 2. * M / R);
}

double Spacetime::d_chi_asymp(double r)
{
    return(isotropic) ? (2 * M * pow(1. +  0.5 * M / r, -5.) / (r * r)) : -(2./3.) * M * pow(1 - 2. * M / r, -4./3.) / (r * r);
}

double Spacetime::d_alpha_asymp(double r)
{
        return (M / (R * R)) / sqrt(1 - 2. * M / r);
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
    if (proca || dirac)
    {
        cout << "WARNING: Spacetime Hamiltonian solver is not implemented for Proca or Dirac fields; aborting request" << endl;
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

    double R_zz_c = n * inv_hww * (R_zz_c_t1 + chris_Wwz * chris_wwz + 2. * chris_Wwz * chris_zww)
                    - 0.5 * inv_hzz * dzz_h_zz + h_zz * d_z(v_c_chris_Z, j)
                    + c_chris_Z * chris_zzz + 3. * inv_hzz * chris_zzz * cZzz;

    double R_ww_c = -0.5 * inv_hzz * dzz_h_ww + 3. * inv_hzz * chris_Wwz * chris_wwz - 0.5 * n * inv_hww * d_h_ww_z
                    + inv_hww * (cZww * chris_zww + 2. * cZww * chris_wwz)
                    + h_ww * c_chris_Z_overz + c_chris_Z * chris_wwz - inv_hww * h_zw_diff;

    R_zz[j] = R_zz_c + R_zz_chi;
    R_ww[j] = R_ww_c + R_ww_chi;

    if (use_CCZ4) //if CCZ4 is used, add supplementary terms to the (extended) Ricci tensor and fill out theta_Z
    {
        double theta = 0.5 * chi * (c_chris_Z - inv_hzz * cZzz - n * inv_hww * cZww);
        theta_Z[j] = theta;
        R_zz[j] += theta * (chi * d_z_h_zz + h_zz * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_zzz * inv_chi;
        R_ww[j] += theta * (chi * d_z_h_ww - h_ww * d_z_chi) * (inv_chi * inv_chi) - 2. * theta * chris_wwz * inv_chi;
    }

    // traceless part of Ricci components
    R_zz_TF[j] = n * (R_zz[j] - h_zz * inv_hww * R_ww[j]) * invDm1;
    R_ww_TF[j] = (R_ww[j] - h_ww * inv_hzz * R_zz[j]) * invDm1;

    const BSSNState& bssn = slice_ptr->states2[j].bssn;

    if (proca)
    {
        // matter quantities via ComplexVectorField helpers; |X_phi| at the centre plays the role of |phi|
        const CVF& cvf = slice_ptr->states2[j].cvf;
        if (j == 0) test_ctr = sqrt(cvf.Xphi_re * cvf.Xphi_re + cvf.Xphi_im * cvf.Xphi_im);

        rho[j] = cvf_model.rho(cvf, bssn);
        j_z[j] = cvf_model.j_z(cvf, bssn);
        S_zz[j] = cvf_model.S_zz(cvf, bssn);
        S_ww[j] = cvf_model.S_ww(cvf, bssn);
    }
    else if (dirac)
    {
        // matter quantities via ComplexSpinorField helpers; |u| at the centre plays the role of |phi|
        const CSpF& spf = slice_ptr->states2[j].spf;
        if (j == 0) test_ctr = sqrt(spf.u_re * spf.u_re + spf.u_im * spf.u_im);
        const CSpF d_z_spf{slice_ptr->d_z(v_u_re, j), slice_ptr->d_z(v_u_im, j), slice_ptr->d_z(v_v_re, j), slice_ptr->d_z(v_v_im, j)};

        rho[j] = spf_model.rho(spf, bssn, d_z_spf, z);
        j_z[j] = spf_model.j_z(spf, bssn, d_z_spf);
        S_zz[j] = spf_model.S_zz(spf, bssn, d_z_spf);
        S_ww[j] = spf_model.S_ww(spf, bssn, d_z_spf, z);
    }
    else
    {
        // |phi|
        double mod_phi = sqrt(phi_re * phi_re + phi_im * phi_im);
        if (j == 0) test_ctr = mod_phi;

        // matter quantities via ComplexScalarField helpers (reuse member to avoid per-point construction)
        CSF csf = slice_ptr->states2[j].csf;
        // Only phi derivatives are used; K_phi derivatives and second derivatives not used in stress-energy
        CSF d_z_csf{d_z_phi_re, d_z_phi_im, 0.0, 0.0};
        CSF d_zz_csf{0.0, 0.0, 0.0, 0.0};

        rho[j] = csf_model.rho(csf, bssn, d_z_csf, d_zz_csf);
        j_z[j] = csf_model.j_z(csf, bssn, d_z_csf, d_zz_csf);
        S_zz[j] = csf_model.S_zz(csf, bssn, d_z_csf, d_zz_csf);
        S_ww[j] = csf_model.S_ww(csf, bssn, d_z_csf, d_zz_csf);
    }

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
BSSNSlice Spacetime::slice_rhs(BSSNSlice* slice_ptr, int skipped_levels)
{
    double n = D - 2.;

    //define return slice and size its states array appropriately
    BSSNSlice rhs;
    rhs.R = R;
    rhs.has_BH = slice_ptr->has_BH;
    rhs.refinement_points = slice_ptr->refinement_points;

    rhs.use_CCZ4 = use_CCZ4;
    rhs.states2.resize(n_gridpoints); // composite State vector (initialized to zeros)

    int max_ref = (refinement_levels.size() == 0 ) ? 1 : refinement_levels[refinement_levels.size() - 1]; //refinement level at outermost bdry
    const int res_fac_outer = (max_ref > 0) ? (1 << (max_ref - 1)) : 1;

    // MPI: only compute RHS at points this rank owns. Under replicated storage,
    // rhs.states2 is still sized n_gridpoints; non-owned slots are left zero
    // (initialized above) and are never read by anyone downstream because the
    // corresponding rank owns them. The upper bound n_gridpoints - 2 stays
    // because the last two active gridpoints are set by sommerfeld_BC on the
    // rank that owns the outer edge.
    const int j_lo = std::max(owned_begin(), 0);
    const int j_hi = std::min(owned_end(), n_gridpoints - 2);
    for (int j = j_lo; j < j_hi; j++)
    {
        if (rhs_mask ? !(*rhs_mask)[j] : (!active_points[j] || refinement_levels[j] > max_ref - skipped_levels))//perform no calculations and skip on inactive points
            continue;
        point_updates++;

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
        rhs.states2[j].bssn.beta =  (evolve_shift)? (beta * d_z_beta + gamma_fac * c_chris_Z - eta * beta): 0.;
       
        if (proca)
        {
            // ComplexVectorField RHS for the Proca variables
            CVF d_z_cvf{d_z(v_X_re, j), d_z(v_X_im, j), d_z(v_E_re, j), d_z(v_E_im, j),
                        d_z(v_Xphi_re, j), d_z(v_Xphi_im, j), d_z(v_Z_re, j), d_z(v_Z_im, j)};

            rhs.states2[j].cvf = cvf_model.rhs_CVF(slice_ptr->states2[j].cvf, slice_ptr->states2[j].bssn, d_z_cvf,
                                                   d_z_alpha, d_z_beta, d_z_chi, d_z_h_zz, d_z_h_ww, z);
            rhs.states2[j].csf = CSF{0.0, 0.0, 0.0, 0.0};
        }
        else if (dirac)
        {
            // ComplexSpinorField RHS for the Dirac variables
            CSpF d_z_spf{d_z(v_u_re, j), d_z(v_u_im, j), d_z(v_v_re, j), d_z(v_v_im, j)};

            rhs.states2[j].spf = spf_model.rhs_CSpF(slice_ptr->states2[j].spf, slice_ptr->states2[j].bssn, d_z_spf,
                                                    d_z_alpha, d_z_chi, d_z_h_ww, z);
            rhs.states2[j].csf = CSF{0.0, 0.0, 0.0, 0.0};
            rhs.states2[j].cvf = CVF{};
        }
        else
        {
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
            rhs.states2[j].cvf = CVF{};
        }

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
            rhs.theta.resize(n_gridpoints);

            rhs.states2[j].bssn.K = beta * d_z_K - chi * inv_hzz * D_zz_alpha[j]  - n * inv_hww * chi * D_ww_alpha[j]
                            + alpha * (chi * inv_hzz * R_zz[j] + n * chi * inv_hww * R_ww[j] + K * (K - 2. * theta)
                            - (D - 1.) * c1 * (1. + c2) * theta + 8. * M_PI * (S[j] - (D - 1.) * rho[j]) / (D - 2.));

            rhs.states2[j].bssn.A_zz += -2. * alpha * theta * A_zz;
            rhs.states2[j].bssn.A_ww += -2. * alpha * theta * A_ww;

            //note this is now \hat{\Gamma}^z, not \tilde{\Gamma}^z. So we first remove BSSN specific damping then add new terms.
            rhs.states2[j].bssn.c_chris_Z +=   sigma_BSSN * aux_constraint[j] - 2. * d_z_alpha * theta * inv_hzz + 2. * alpha * inv_hzz * d_z(v_theta, j)
                                        - 2. * c1 * alpha * theta_Z[j] / chi - 4. * alpha * K * theta_Z[j] * invDm1 / chi
                                        + sigma_BSSN * theta_Z[j]  * ((3. - D) * d_z_beta + 2. * n * beta_z) / chi;

            rhs.states2[j].bssn.alpha +=  2. * alpha * alpha * f * theta; //4. * alpha * theta; //corrects for K -> K - 2 * theta in 1+log slicing
            //rhs.states2[j].bssn.beta +=  (evolve_shift)? ( -2. * gamma_fac * theta_Z[j] / chi): 0.; //replace theta hat -> theta tilde in CCZ4 moving puncture gauge; seems to produce instability...

            if (!theta_off)
                rhs.theta[j] = beta * d_z(v_theta, j) + 0.5 * alpha * (chi * inv_hzz * R_zz[j] - inv_hzz * inv_hzz * A_zz * A_zz + n * K * K * invDm1
                            + n * (chi * inv_hww * R_ww[j] - A_ww * A_ww / (h_ww * h_ww))  - 2. * K * theta - 2. * theta_Z[j] * d_z_alpha / alpha
                            - c1 * (D + c2 * n) * theta  - 16. * M_PI * rho[j]);
            else
                rhs.theta[j] = 0.;
        }

        // add damping away from edges
        if (damping_factor != 0. && j < n_gridpoints - 3 * res_fac_outer)
            apply_damping_at_point(j, slice_ptr, rhs.states2, rhs.theta);
    }

    //apply Sommerfeld boundary conditions at outer boundary
    if (rhs_do_sommerfeld)
        sommerfeld_BC(rhs, res_fac_outer, slice_ptr);

    return rhs;
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

        if (D == 4.)
            mass += -0.25 * quad_weights[j] * (Ham[j]- h_ZZ[j] * chi * R_zz[j] - n * h_WW[j] * chi * R_ww[j]) * z * z * pow(chi, -1.25) ; //-1.25 is spurious! Still haven't figured out why this works...
        else
             mass += -0.125 * pow(M_PI, 0.5 * (D - 3.))  * quad_weights[j] * (Ham[j]- h_ZZ[j] * chi * R_zz[j] - n * h_WW[j] * chi * R_ww[j]) * pow(z, D - 2.) * pow(chi, -1.25 - 0.35 * (D - 4.)) / tgamma(0.5 * (D - 1.));
    }
    return bss_mpi::sum(mass);
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

        if (dirac) // 4 pi z^2 n sqrt(gamma), n = (|u|^2 + |v|^2) / (2 pi) the number density (D = 4)
            charge += 4. * M_PI * quad_weights[j] * spf_model.number_density(slice_ptr->states2[j].spf) * z * z * pow(chi, -1.5);
        else if (proca) // Omega_{D-2} * charge density * sqrt(gamma) z^{D-2}, with sqrt(gamma) = chi^{(1-D)/2} as below
        {
            const double sphere_area = 2. * pow(M_PI, 0.5 * (D - 1.)) / tgamma(0.5 * (D - 1.));
            charge += sphere_area * quad_weights[j] * cvf_model.charge_density(slice_ptr->states2[j].cvf) * pow(z, n) * pow(chi, 0.5 * (1. - D));
        }
        else if (D == 4.)
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

        // Dirac field contribution (stored in E_phi, the "complex field energy")
        if (dirac)
        {
            const CSpF& spf = slice_ptr->states2[j].spf;
            const CSpF d_z_spf{slice_ptr->d_z(v_u_re, j), slice_ptr->d_z(v_u_im, j), slice_ptr->d_z(v_v_re, j), slice_ptr->d_z(v_v_im, j)};
            const double rho_D = spf_model.rho(spf, bssn, d_z_spf, z);
            e_phi += quad_weights[j] * sphere_area * rho_D * std::pow(z, D - 2.) * std::pow(chi, 0.5 * (1.0 - D));
        }
        // Proca field contribution (stored in E_phi, the "complex field energy")
        else if (proca)
        {
            const double rho_A = cvf_model.rho(slice_ptr->states2[j].cvf, bssn);
            e_phi += quad_weights[j] * sphere_area * rho_A * std::pow(z, D - 2.) * std::pow(chi, 0.5 * (1.0 - D));
        }
        // Complex scalar field contribution
        else
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

    const BSSNSlice& s = *current_slice_ptr;
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
        contribution = 0.5 * std::sqrt(h_ww0 / chi0) * ah_radius;
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
    double gauss_sq = 0., gauss_ref_sq = 0.; //Proca Gauss constraint and its reference scale mu^2 |X_phi|

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
            Ham_L2 += quad_weights[j] * Ham[j] * Ham[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
            Mom_L2 += quad_weights[j] * Mom_Z[j] * Mom_Z[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
        }

        if (proca)
        {
            const CVF& v = slice_ptr->states2[j].cvf;
            CVF d_z_cvf{};
            d_z_cvf.E_re = slice_ptr->d_z(v_E_re, j);
            d_z_cvf.E_im = slice_ptr->d_z(v_E_im, j);
            double C_re, C_im;
            cvf_model.gauss_constraint(v, slice_ptr->states2[j].bssn, d_z_cvf, slice_ptr->d_z(v_chi, j),
                                       slice_ptr->d_z(v_h_zz, j), slice_ptr->d_z(v_h_ww, j), z, C_re, C_im);
            Gauss[j] = sqrt(C_re * C_re + C_im * C_im);

            if (z > start_val && z < end_val)
            {
                const double w = quad_weights[j] * pow(z, D - 2.) * pow(chi, (1. - D) / 2.);
                gauss_sq += w * Gauss[j] * Gauss[j];
                gauss_ref_sq += w * pow(mu, 4) * (v.Xphi_re * v.Xphi_re + v.Xphi_im * v.Xphi_im);
            }
        }
    }

    //take square root to get L^2 norms; normalize by 16 * pi * central energy density.
    // MPI: combine per-rank sums-of-squares before sqrt.
    Ham_L2 = sqrt(bss_mpi::sum(Ham_L2)) / (16. * M_PI * rho0_init);
    Mom_L2 = sqrt(bss_mpi::sum(Mom_L2)) / (16. * M_PI * rho0_init);

    //Proca Gauss constraint, relative to the L2 norm of mu^2 X_phi over the same region
    if (proca)
    {
        const double num = bss_mpi::sum(gauss_sq), den = bss_mpi::sum(gauss_ref_sq);
        Gauss_L2 = (den > 0.) ? sqrt(num / den) : sqrt(num);
    }
}

//writes diagnostic quantities to output file. maybe just make appropriate call to write_slice() instead?
// MPI: collective. All ranks pack their owned rows; rank 0 gathers and writes.
void Spacetime:: write_current_slice(std::string file_name)
{
    const BSSNSlice& s = *current_slice_ptr;
    dr = R / ((double)n_gridpoints - 1.0);

    // 28 doubles per point (all state fields; theta packed as 0 if not CCZ4); then 8 Proca variables, written only for
    // Proca runs, and 4 Dirac variables, written only for Dirac runs, so that boson star output is unchanged.
    constexpr int NCOL = 28;
    std::vector<double> local(local_N * NCOL);
    for (int k = 0; k < local_N; ++k)
    {
        const int j = owned_begin() + k;
        double* p = &local[k * NCOL];
        const CVF& v = s.states2[j].cvf;
        p[16] = v.X_re;    p[17] = v.X_im;    p[18] = v.E_re; p[19] = v.E_im;
        p[20] = v.Xphi_re; p[21] = v.Xphi_im; p[22] = v.Z_re; p[23] = v.Z_im;
        const CSpF& sp = s.states2[j].spf;
        p[24] = sp.u_re; p[25] = sp.u_im; p[26] = sp.v_re; p[27] = sp.v_im;
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
                  << "   " << p[14] << "    " << p[15];
        if (proca)
            for (int c = 16; c < 24; ++c) data_file << "    " << p[c];
        if (dirac)
            for (int c = 24; c < 28; ++c) data_file << "    " << p[c];
        data_file << endl;
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
    //   Ham, Mom_Z, det_h, aux_test, A, ricci_4_val, r4m, Gauss (Proca only, written only for Proca runs)
    constexpr int NCOL = 8;
    std::vector<double> local(local_N * NCOL, 0.0);
    double local_r4max = 0.0;

    const int j_lo = std::max(owned_begin(), 0);
    const int j_hi = std::min(owned_end(), length - 2);
    for (int j = j_lo; j < j_hi; j++)
    {
        if (!active_points[j]) continue;

        const double A   = central_amplitude(*current_slice_ptr, j); //|phi|, or |X_phi| for Proca
        const double r4m = ricci_4_matter(j);
        if (r4m > local_r4max) local_r4max = r4m;

        double* p = &local[(j - owned_begin()) * NCOL];
        p[0] = Ham[j];
        p[1] = Mom_Z[j];
        p[2] = det_h[j];
        p[3] = aux_test[j];
        p[4] = A;
        p[5] = ricci_4_val(j);
        p[6] = r4m;
        p[7] = proca ? Gauss[j] : 0.;
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
        if (proca) data_file << "    " << p[7];
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
        fill_parameter(current_line, "proca_damping = ", proca_damping, quiet);

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
    Mom_Z.resize(n_gridpoints);
    det_h.resize(n_gridpoints);
    Gauss.assign(n_gridpoints, 0.);

    if (use_CCZ4)
    {
        theta_Z.resize(n_gridpoints);
    }

    dtK_L2 = 0;

    last_active_j = n_gridpoints - 1; //find last active gridpoint
    while (!active_points[last_active_j])
        last_active_j--;

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
    proca = false;
    dirac = false;
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

    setup_before_reading(skip_read);

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

    setup_after_reading();
}


//parameter reading, consistency checks and grid setup shared by boson star and Proca star initialisation.
//Requires the parameters inherited from the star (D, n_gridpoints, R, stop_time, ...) to have been set.
void Spacetime::setup_before_reading(bool skip_read)
{
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
    

    if (D != 4.)
        spatially_varying_BC = 0; //only try this for D = 4 for now

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
}

//trims the grid, sets up domain decomposition/auxiliary arrays and computes initial diagnostics; shared by boson star
//and Proca star initialisation. Must be called after the initial slice has been read in.
void Spacetime::setup_after_reading()
{
    grid_offset = 0.5 * (double)cell_centered;

    if (!run_quietly)
        cout << (dirac ? "Read DS data" : (proca ? "Read PS data" : "Read BS data")) << endl;

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
        slices[0].read_checkpoint(start_time, n_gridpoints, proca, dirac);

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

    if (mpi_rank == 0)
    {
        std::ofstream dynamical_file{"dynamical_constants.dat"};
        dynamical_file << mass0 << "    " << charge0 << "    " << mass0 - mu * charge0;

        if (!run_quietly)
            cout << "Dynamical N =" << charge0 << ", M = "  << mass0 << " and binding energy is " << mass0 - mu * charge0  << endl;
    }
}

//|phi| for the complex scalar field, or |X_phi| = |n^a A_a| for the Proca field (which equals f(0)/alpha(0) = A_central for
//a stationary Proca star), at index j of a slice
double Spacetime::central_amplitude(const BSSNSlice& slice, int j) const
{
    if (dirac) //|u| = |upper spinor amplitude|, which equals A_central = g0 for a stationary Dirac star
        return sqrt(pow(slice.states2[j].spf.u_re, 2) + pow(slice.states2[j].spf.u_im, 2));
    if (proca)
        return sqrt(pow(slice.states2[j].cvf.Xphi_re, 2) + pow(slice.states2[j].cvf.Xphi_im, 2));
    return sqrt(pow(slice.states2[j].csf.phi_re, 2) + pow(slice.states2[j].csf.phi_im, 2));
}

double Spacetime::central_phase(const BSSNSlice& slice, int j) const
{
    if (dirac)
        return std::arg(std::complex<double>(slice.states2[j].spf.u_re, slice.states2[j].spf.u_im));
    if (proca)
        return std::arg(std::complex<double>(slice.states2[j].cvf.Xphi_re, slice.states2[j].cvf.Xphi_im));
    return std::arg(std::complex<double>(slice.states2[j].csf.phi_re, slice.states2[j].csf.phi_im));
}

//read data from a solved ProcaStar to spacetime and construct initial time slice (Proca analogue of initialize(BosonStar&)).
//Only the free (mini) Proca potential is supported in evolution.
void Spacetime::initialize(ProcaStar& proca_star, bool skip_read)
{
    proca = true;
    dirac = false;
    wave_mode = 0;
    ah_radius = -1.;
    D = proca_star.D;

    //inherit parameters from PS
    n_gridpoints = proca_star.n_gridpoints;
    R = proca_star.R;
    courant_factor = proca_star.courant_factor;
    stop_time = proca_star.stop_time;
    mu = proca_star.mu;
    lambda = proca_star.lambda;
    sigma = proca_star.sigma;
    solitonic = proca_star.solitonic;
    omega = proca_star.omega;
    isotropic = proca_star.isotropic;
    M = proca_star.M;
    r_99 = proca_star.r_99;
    BS_perturbed = 0;

    setup_before_reading(skip_read);

    if (proca_star.self_interacting())
    {
        if (mpi_rank == 0)
            cerr << "ERROR: evolution of self-interacting Proca stars (lambda != 0 or solitonic = 1) is not implemented;"
                 << " only the free (mini) Proca field is supported." << endl;
        bss_mpi::abort_all(1);
    }
    if (read_thinshell || critical_study || make_tangherlini)
    {
        if (mpi_rank == 0)
            cerr << "ERROR: read_thinshell, critical_study and make_tangherlini are not supported for Proca stars." << endl;
        bss_mpi::abort_all(1);
    }
    if (damping_factor <= 0. && mpi_rank == 0)
        cout << "WARNING: the first-order Proca system is unstable at the origin without Kreiss-Oliger dissipation;"
             << " damping_factor of about 1 is recommended for Proca evolutions." << endl;
    if (add_real_field || run_spacetime_solver)
    {
        if (mpi_rank == 0)
            cout << "WARNING: add_real_field and run_spacetime_solver are not supported for Proca stars; ignoring them." << endl;
        add_real_field = 0;
        run_spacetime_solver = 0;
    }

    const int ps_n_gridpoints_orig = proca_star.n_gridpoints;
    bool bumped_resolution = false;

    //solve PS at higher resolution and read in data to first slice, as for the BS
    if (BS_resolution_factor > 1 && start_time == 0)
    {
        if (mpi_rank == 0)
        {
            proca_star.n_gridpoints = proca_star.n_gridpoints * BS_resolution_factor - BS_resolution_factor + 1;

            if (!run_quietly)
                cout << " \nRe-solving PS at higher resolution with " << proca_star.n_gridpoints <<  endl;

            proca_star.solve(1);
            proca_star.write_field();
            if (isotropic)
            {
                proca_star.fill_isotropic_arrays();
                proca_star.write_isotropic();
            }
        }
        proca_star.broadcast();

        omega = proca_star.omega;
        bumped_resolution = true;
    }

    //optional perturbation of the star (scales the matter fields and re-solves the metric)
    if (proca_star.perturb && start_time == 0)
    {
        if (mpi_rank == 0)
        {
            proca_star.add_perturbation(proca_star.perturb_amp, proca_star.perturb_spread, proca_star.perturb_center, run_quietly);
            proca_star.write_field("PSdata_perturbed.dat");
            if (isotropic)
            {
                proca_star.fill_isotropic_arrays();
                proca_star.write_isotropic();
            }
        }
        proca_star.broadcast();
        M = proca_star.M;
        r_99 = proca_star.r_99;
    }

    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    int num_timesteps = ceil(stop_time / dt);

    slices.resize(std::min(num_timesteps + 1, max_stored_slices));
    slices[0].active_points = active_points;
    slices[0].refinement_points = refinement_points;
    slices[0].use_CCZ4 = use_CCZ4;

    if (start_time == 0)
        slices[0].read_PS_data(proca_star, BS_resolution_factor, isotropic, cell_centered);

    if (bumped_resolution)
        proca_star.n_gridpoints = ps_n_gridpoints_orig;

    setup_after_reading();
}

//read data from a solved DiracStar to spacetime and construct initial time slice (Dirac analogue of initialize(BosonStar&)).
//Only the free field in D = 4 is supported.
void Spacetime::initialize(DiracStar& dirac_star, bool skip_read)
{
    proca = false;
    dirac = true;
    wave_mode = 0;
    ah_radius = -1.;
    D = dirac_star.D;

    //inherit parameters from the Dirac star
    n_gridpoints = dirac_star.n_gridpoints;
    R = dirac_star.R;
    courant_factor = dirac_star.courant_factor;
    stop_time = dirac_star.stop_time;
    mu = dirac_star.mu;
    lambda = dirac_star.lambda;
    sigma = dirac_star.sigma;
    solitonic = dirac_star.solitonic;
    omega = dirac_star.omega;
    isotropic = dirac_star.isotropic;
    M = dirac_star.M;
    r_99 = dirac_star.r_99;
    BS_perturbed = 0;

    setup_before_reading(skip_read);

    if (D != 4.)
    {
        if (mpi_rank == 0)
            cerr << "ERROR: Dirac star evolution is only implemented for D = 4." << endl;
        bss_mpi::abort_all(1);
    }
    if (read_thinshell || critical_study || make_tangherlini)
    {
        if (mpi_rank == 0)
            cerr << "ERROR: read_thinshell, critical_study and make_tangherlini are not supported for Dirac stars." << endl;
        bss_mpi::abort_all(1);
    }
    if (damping_factor <= 0. && mpi_rank == 0)
        cout << "WARNING: the Dirac system is first order in space, like the Proca system; Kreiss-Oliger dissipation"
             << " (damping_factor of about 1) is recommended for Dirac evolutions." << endl;
    if (add_real_field || run_spacetime_solver)
    {
        if (mpi_rank == 0)
            cout << "WARNING: add_real_field and run_spacetime_solver are not supported for Dirac stars; ignoring them." << endl;
        add_real_field = 0;
        run_spacetime_solver = 0;
    }

    const int ds_n_gridpoints_orig = dirac_star.n_gridpoints;
    bool bumped_resolution = false;
    BS_perturbed = dirac_star.perturb;

    //solve the star at higher resolution and read in data to first slice, as for the BS
    if (BS_resolution_factor > 1 && start_time == 0)
    {
        if (mpi_rank == 0)
        {
            dirac_star.n_gridpoints = dirac_star.n_gridpoints * BS_resolution_factor - BS_resolution_factor + 1;

            if (!run_quietly)
                cout << " \nRe-solving DS at higher resolution with " << dirac_star.n_gridpoints <<  endl;

            dirac_star.solve(1);
            dirac_star.write_field();
            if (isotropic)
            {
                dirac_star.fill_isotropic_arrays();
                dirac_star.write_isotropic();
            }
        }
        dirac_star.broadcast();

        omega = dirac_star.omega;
        bumped_resolution = true;
    }

    //optional perturbation of the star (scales the spinor amplitudes and re-solves the metric)
    if (dirac_star.perturb && start_time == 0)
    {
        if (mpi_rank == 0)
        {
            dirac_star.add_perturbation(dirac_star.perturb_amp, dirac_star.perturb_spread, dirac_star.perturb_center, run_quietly);
            dirac_star.write_field("DSdata_perturbed.dat");
            if (isotropic)
            {
                dirac_star.fill_isotropic_arrays();
                dirac_star.write_isotropic();
            }
        }
        dirac_star.broadcast();
        M = dirac_star.M;
        r_99 = dirac_star.r_99;
    }

    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    int num_timesteps = ceil(stop_time / dt);

    slices.resize(std::min(num_timesteps + 1, max_stored_slices));
    slices[0].active_points = active_points;
    slices[0].refinement_points = refinement_points;
    slices[0].use_CCZ4 = use_CCZ4;

    if (start_time == 0)
        slices[0].read_DS_data(dirac_star, BS_resolution_factor, isotropic, cell_centered);

    if (bumped_resolution)
        dirac_star.n_gridpoints = ds_n_gridpoints_orig;

    setup_after_reading();
}

// ================= Berger-Oliger subcycling of refinement levels (adaptive_timestep = 1) =================
void Spacetime::setup_subcycling()
{
    const int B = subcycle_buffer;
    n_levels_sub = refinement_levels.empty() ? 1 : refinement_levels[last_active_j];
    sublevels.assign(n_levels_sub, SubLevel());
    sub_start.assign(n_levels_sub, BSSNSlice());
    sub_stages.assign(n_levels_sub, std::vector<BSSNSlice>());
    sub_step.assign(n_levels_sub, 0.);

    for (int k = 1; k <= n_levels_sub; k++)
    {
        SubLevel& Lv = sublevels[k - 1];
        Lv.k = k; Lv.s = 1 << (k - 1);
        Lv.a = -1; Lv.b = -1;
        for (int j = 0; j <= last_active_j; j++)
            if (refinement_levels[j] == k) { if (Lv.a < 0) Lv.a = j; Lv.b = j; }
        int q0 = -1;
        for (int j = Lv.a; j >= 0 && j <= Lv.b; j++) if (active_points[j]) { q0 = j; break; }
        if (Lv.a < 0 || q0 < 0) { cerr << "ERROR: subcycling setup failed for refinement level " << k << endl; exit(1); }
        auto aligned = [&](int j) { return ((j - q0) % Lv.s + Lv.s) % Lv.s == 0; };

        Lv.mask.assign(n_gridpoints, 0); Lv.rf.assign(n_gridpoints, 0);
        for (int j = Lv.a; j <= Lv.b; j++)
            if (aligned(j))
            {
                if (!active_points[j]) { cerr << "ERROR: subcycling: point " << j << " of refinement level " << k << " is inactive" << endl; exit(1); }
                Lv.writeback.push_back(j);
            }
        // inner buffer: points of the finer level (at this level's spacing), stopping short of the origin so that the
        // stencils stay well defined; outer buffer: points of the coarser level k+1, filled by interpolation
        std::vector<int> inner;
        if (k > 1)
        {
            if (Lv.a - B * Lv.s < 3 * Lv.s)
            {
                cerr << "ERROR: subcycling (adaptive_timestep = 1): refinement level " << k << " starts too close to the origin to hold its "
                     << subcycle_buffer << "-point inner buffer; move the first refinement points outwards" << endl;
                exit(1);
            }
            for (int j = Lv.a - 1; j >= Lv.a - B * Lv.s; j--)
                if (aligned(j)) inner.push_back(j);
        }
        if (k < n_levels_sub)
        {
            for (int j = Lv.b + 1; j <= std::min(Lv.b + B * Lv.s, last_active_j); j++)
                if (aligned(j)) Lv.outer_buffer.push_back(j);
            for (int j : Lv.outer_buffer)
                if (refinement_levels[j] != k + 1)
                {
                    cerr << "ERROR: subcycling (adaptive_timestep = 1): refinement level " << k + 1 << " is too narrow to hold the "
                         << subcycle_buffer << "-point buffer of level " << k << "; move the refinement points further apart" << endl;
                    exit(1);
                }
        }
        for (int j : inner) Lv.evolve.push_back(j);
        for (int j : Lv.writeback) Lv.evolve.push_back(j);
        for (int j : Lv.outer_buffer) Lv.evolve.push_back(j);
        for (int j : Lv.evolve) { Lv.mask[j] = 1; Lv.rf[j] = Lv.s; }

        if (!run_quietly && mpi_rank == 0)
            cout << "Subcycling level " << k << ": spacing " << Lv.s << ", points " << Lv.a << ".." << Lv.b << " (" << Lv.writeback.size()
                 << " evolved), inner buffer " << inner.size() << ", outer buffer " << Lv.outer_buffer.size() << endl;
    }
}

// Fills the outer buffer of level k in W for its substep m (0 or 1) within the current step of level k+1. Values of level
// k+1 at the intermediate time come from RK4 dense output over its last step; points between level k+1's points are
// filled by 6-point Lagrange interpolation over the merged points of both levels (the last points of level k, already at
// the substep's start time in cur, followed by those of level k+1), centred on the target so nothing is extrapolated.
void Spacetime::fill_outer_buffer(BSSNSlice& W, const BSSNSlice& cur, int k, int m)
{
    const SubLevel& Lv = sublevels[k - 1];
    const SubLevel& Lc = sublevels[k]; // level k+1
    const int ic = k;                  // storage index of level k+1
    const bool ccz4 = use_CCZ4 && !cur.theta.empty();

    // dense output weights for y(T + th D) = y(T) + D (b1 k1 + b2 k2 + b3 k3 + b4 k4), b2 = b3
    const double th = 0.5 * m;
    const double b1 = th - 1.5 * th * th + 2. * th * th * th / 3.;
    const double b23 = th * th - 2. * th * th * th / 3.;
    const double b4 = -0.5 * th * th + 2. * th * th * th / 3.;
    const std::vector<BSSNSlice>& kk = sub_stages[ic];
    const BSSNSlice& y0 = sub_start[ic];
    const double D = sub_step[ic];
    if (kk.size() != 4) { cerr << "ERROR: subcycling: level " << k + 1 << " has no stored step" << endl; exit(1); }

    struct Node { int j; State s; double theta; };
    std::vector<Node> nodes;
    const int n_own = 8;
    for (size_t q = (Lv.writeback.size() > (size_t)n_own ? Lv.writeback.size() - n_own : 0); q < Lv.writeback.size(); q++)
    {
        const int j = Lv.writeback[q];
        nodes.push_back({j, cur.states2[j], ccz4 ? cur.theta[j] : 0.});
    }
    const int jmax = Lv.outer_buffer.empty() ? Lv.b : Lv.outer_buffer.back();
    int extra = 0;
    for (int j : Lc.writeback)
    {
        if (j <= Lv.b) continue;
        Node nd{j, y0.states2[j] + D * (b1 * kk[0].states2[j] + b23 * (kk[1].states2[j] + kk[2].states2[j]) + b4 * kk[3].states2[j]), 0.};
        if (ccz4) nd.theta = y0.theta[j] + D * (b1 * kk[0].theta[j] + b23 * (kk[1].theta[j] + kk[2].theta[j]) + b4 * kk[3].theta[j]);
        nodes.push_back(nd);
        if (j > jmax && ++extra >= 3) break;
    }
    if (nodes.size() < 6) { cerr << "ERROR: subcycling: too few interpolation points for the buffer of level " << k << endl; exit(1); }

    for (int j : Lv.outer_buffer)
    {
        auto it = std::lower_bound(nodes.begin(), nodes.end(), j, [](const Node& nd, int jj) { return nd.j < jj; });
        const int hi = it - nodes.begin();
        if (it != nodes.end() && it->j == j) { W.states2[j] = it->s; if (ccz4) W.theta[j] = it->theta; continue; }
        const int lo = std::max(0, std::min(hi - 3, (int)nodes.size() - 6));
        State v{};
        double tv = 0.;
        for (int a = 0; a < 6; a++)
        {
            double w = 1.;
            for (int b = 0; b < 6; b++)
                if (b != a) w *= (double)(j - nodes[lo + b].j) / (double)(nodes[lo + a].j - nodes[lo + b].j);
            v = v + w * nodes[lo + a].s;
            tv += w * nodes[lo + a].theta;
        }
        W.states2[j] = v;
        if (ccz4) W.theta[j] = tv;
    }
}

// one RK4 step of size Dk restricted to the points of level Lv, using that level's stencil spacing; optionally stores the stages
void Spacetime::rk4_masked(BSSNSlice& W, const SubLevel& Lv, double Dk, std::vector<BSSNSlice>* stages)
{
    rhs_mask = &Lv.mask; rhs_res_fac = &Lv.rf; rhs_do_sommerfeld = (Lv.k == n_levels_sub);
    W.forced_res_fac = &Lv.rf;

    current_slice_ptr = &W;
    BSSNSlice s1 = slice_rhs(&W);
    BSSNSlice t1 = W + (0.5 * Dk) * s1; t1.forced_res_fac = &Lv.rf;
    current_slice_ptr = &t1;
    BSSNSlice s2 = slice_rhs(&t1);
    BSSNSlice t2 = W + (0.5 * Dk) * s2; t2.forced_res_fac = &Lv.rf;
    current_slice_ptr = &t2;
    BSSNSlice s3 = slice_rhs(&t2);
    BSSNSlice t3 = W + Dk * s3; t3.forced_res_fac = &Lv.rf;
    current_slice_ptr = &t3;
    BSSNSlice s4 = slice_rhs(&t3);

    W = W + (Dk / 6.) * (s1 + 2. * s2 + 2. * s3 + s4);
    W.forced_res_fac = nullptr;
    if (stages) { stages->assign(4, BSSNSlice()); (*stages)[0] = s1; (*stages)[1] = s2; (*stages)[2] = s3; (*stages)[3] = s4; }
    rhs_mask = nullptr; rhs_res_fac = nullptr; rhs_do_sommerfeld = true;
}

// advances level k by Dk and then, recursively, all finer levels by two steps of Dk / 2; m = index of this step within
// the step of level k+1
void Spacetime::subcycle_level(int k, BSSNSlice& cur, double Dk, int m)
{
    const SubLevel& Lv = sublevels[k - 1];

    BSSNSlice W = cur;
    if (k < n_levels_sub)
        fill_outer_buffer(W, cur, k, m);

    if (k > 1)
    {   // keep what level k-1 needs to fill its outer buffer
        sub_start[k - 1] = cur;
        sub_step[k - 1] = Dk;
        rk4_masked(W, Lv, Dk, &sub_stages[k - 1]);
    }
    else
        rk4_masked(W, Lv, Dk);

    for (int j : Lv.writeback)
    {
        cur.states2[j] = W.states2[j];
        if (use_CCZ4 && W.theta.size() == cur.theta.size()) cur.theta[j] = W.theta[j];
        if (cur.states2[j].bssn.chi < min_chi) cur.states2[j].bssn.chi = min_chi;
    }
    make_A_traceless(&cur);

    if (k > 1)
    {
        subcycle_level(k - 1, cur, 0.5 * Dk, 0);
        subcycle_level(k - 1, cur, 0.5 * Dk, 1);
    }
}

//Time evolution! Must have initialized using a constructed BosonStar first.
void Spacetime::evolve()
{
    dr = R / (n_gridpoints - 1.0);
    dt = courant_factor * dr;
    if (critical_study) ricci_4_ctr_max = 0.0;

    const bool subcycling = adaptive_timestep && !refinement_points.empty();
    if (subcycling)
    {
        setup_subcycling();
        dt = courant_factor * dr * (1 << (n_levels_sub - 1)); //step of the coarsest level, i.e. of the time loop
    }
    point_updates = 0;

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
                "   " << "complex field energy" << "   " << "real field energy" << "   " << "central 4D ricci scalar" << "   " << "phi_re_central" << "   " << "AH mass";
        if (proca) constraints_file << "   " << "Gauss_L2"; //for Proca runs, A_central/phi_re_central/field energy refer to X_phi and the Proca field
        constraints_file << endl;
    }

    // A0 is the initial central field amplitude, used later for migration
    // detection and the subcritical check. Only rank 0 has states2[0] valid
    // in the general case; broadcast so all ranks share the same value.
    double A0 = 0.0;
    if (owns_origin())
        A0 = central_amplitude(slices[0], 0);
    bss_mpi::bcast_doubles(&A0, 1);

    if (!run_quietly && mpi_rank == 0)
        cout <<" \n Will evolve with up to " << num_timesteps << " time steps \n" << endl;
        
    //s_i are returned RHS's, t represents temporary RHS + current_slice quantities that must be stored so derivatives can be accessed
    BSSNSlice s1, s2, s3, s4, t1, t2, t3, update_slice;

    double phase_last = 0.0; //central phase on the previous step, for the central frequency diagnostic

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

        int skipped_levels = 0; //skipping of outer levels is superseded by subcycling (adaptive_timestep = 1)

        //fill out array until we've reached maximum number of stored slices, then update last element + rotate at end.
        int n = (time_step > max_stored_slices - 2) ? (max_stored_slices - 2) : time_step;
        double phase_ctr = 0.0;
        // MPI: states2[0] is only valid on rank 0. Compute central-point diagnostics
        // on rank 0 and broadcast the few scalars other ranks need for control flow.
        if (owns_origin())
            phase_ctr = central_phase(slices[n], 0);

        if (n > 0) compute_dtK(n);
        if (subcycling)
        {   // the auxiliary arrays (R_zz, ...) used by slice_mass were last written by the masked RHS evaluations of
            // the finest level, which leave stale or buffer-contaminated values elsewhere; refresh them on this slice
            current_slice_ptr = &slices[n];
            compute_auxiliary_quantities(current_slice_ptr);
        }
        M = slice_mass(current_slice_ptr); // maybe remove if causes bad bdry oscillations?

        // Central-point scalars needed on every rank for control flow.
        // Layout: [0]=A_ctr, [1]=has_BH_now, [2]=chi_ctr, [3]=alpha_ctr,
        //         [4]=phi_re_ctr, [5]=phase_ctr, [6]=test_ctr.
        // (Rank 0 also uses [2..5] to write constraint_norms.dat; the bundle keeps it in one bcast.)
        double central[7] = {0., 0., 0., 0., 0., 0., 0.};
        if (owns_origin())
        {
            const double Ap0 = central_amplitude(slices[n], 0);
            const double Ap1 = central_amplitude(slices[n], 1);
            central[0] = (cell_centered) ? (Ap0 - 0.5 * dr * (Ap1 - Ap0)) : Ap0;
            central[1] = (make_tangherlini || slices[n].states2[0].bssn.chi < 10 * min_chi) ? 1.0 : 0.0;
            central[2] = slices[n].states2[0].bssn.chi;
            central[3] = slices[n].states2[0].bssn.alpha;
            central[4] = dirac ? slices[n].states2[0].spf.u_re : (proca ? slices[n].states2[0].cvf.Xphi_re : slices[n].states2[0].csf.phi_re);
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

        //central frequency from the phase change over the last step, with the phase difference wrapped into (-pi, pi].
        //Boson star initial data rotate as e^{+i omega t} and Proca/Dirac star data as e^{-i omega t}; the sign is chosen so
        //that a stationary star of either kind reports omega > 0.
        double dphase = phase_ctr - phase_last;
        while (dphase > M_PI) dphase -= 2. * M_PI;
        while (dphase <= -M_PI) dphase += 2. * M_PI;
        double omega_approx = (time_step > 0) ? ((proca || dirac) ? -1. : 1.) * dphase / dt : 0.;
        phase_last = phase_ctr;

        double ricci_4 = (critical_study ? ricci_4_grid_max : 0.0);

        if (critical_study && ricci_4 > ricci_4_ctr_max && t > 10)
            ricci_4_ctr_max = ricci_4;

        if (time_step % write_CN_interval == 0) //write time-dependent diagnostics to constraint_norms.dat
        {
            compute_scalar_energies(current_slice_ptr);
            const double charge_now = slice_charge(current_slice_ptr);
            const double ah_mass_now = get_ah_mass(); // 0 on non-owners since ah_radius<0 when AH search is off
            if (mpi_rank == 0)
            {
                constraints_file << std::setprecision (10) << start_time + dt * time_step << "   " << Ham_L2
                << "   " << Mom_L2 <<  "   " << chi_central << "   "
                << A_ctr << "   "  << ah_radius << "   " <<  M << "   "  << charge_now
                << "   " << dtK_L2 << "   " << omega_approx << "   " << alpha_central <<
                "   " << E_phi << "   " << E_psi << "   " << ricci_4 << "   " << phi_re_ctr << "   " << ah_mass_now;
                if (proca) constraints_file << "   " << Gauss_L2;
                constraints_file << endl;
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

        if (subcycling)
        {
            if (time_step % write_CN_interval == 0)
                compute_diagnostics(current_slice_ptr); //auxiliary quantities on this slice were refreshed above
            BSSNSlice next = slices[n];
            subcycle_level(n_levels_sub, next, dt, 0);
            slices[n + 1] = next;
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
        }

        //enforce that A is traceless
        current_slice_ptr = &slices[n + 1];
        make_A_traceless(current_slice_ptr);

        //enforce minimum chi on owned points (ghosts get refreshed by exchange_ghosts below)
        for (int j = owned_begin(); j < owned_end(); j++)
        {
            State& st = slices[n + 1].states2[j];
            if (st.bssn.chi < min_chi) st.bssn.chi = min_chi;
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
                    super_M = M;
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
            if (ah_radius > 0.) super_ah_radius = ah_radius;
            super_M = M;
        }

        //if past_super_time > 0, keep the last known good AH radius / mass while extending, then break when the window elapses
        if (critical_study && super_detect_time >= 0.0)
        {
            if (ah_radius > 0.) super_ah_radius = ah_radius;
            super_M = M;
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
    if (subcycling && !run_quietly && mpi_rank == 0)
        cout << "Subcycling: " << point_updates << " RHS point evaluations" << endl;

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
    supercritical_file << "#tuning_param   critical_time   ah_radius_final   M_final" << endl;

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
                               << hi_guess<< "   " << critical_time << "   " << super_ah_radius << "   " << super_M << endl;

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
                               << mid << "   " << critical_time << "   " << super_ah_radius << "   " << super_M << endl;
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

    double ah_indicator = sqrt(s.chi / s.h_zz) *
                          (2. / z + (d_z(v_h_ww, j) - s.h_ww * d_z(v_chi, j) / s.chi) / s.h_ww)
                          - (D - 2.) * (s.A_ww + s.K * s.h_ww / (D - 1.)) / s.h_ww;

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