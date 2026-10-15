#ifndef MPI_HELPERS_H_
#define MPI_HELPERS_H_

// Thin wrapper around MPI so the codebase builds either serially (plain g++)
// or with MPI (mpicxx -DUSE_MPI). In serial builds, mpi_rank()==0 and
// mpi_size()==1, and every collective is a no-op.

#include <cstdlib>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace bss_mpi
{
    inline void init(int* argc, char*** argv)
    {
#ifdef USE_MPI
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized)
            MPI_Init(argc, argv);
#else
        (void)argc; (void)argv;
#endif
    }

    inline void finalize()
    {
#ifdef USE_MPI
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized)
            MPI_Finalize();
#endif
    }

    inline int rank()
    {
#ifdef USE_MPI
        int r = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &r);
        return r;
#else
        return 0;
#endif
    }

    inline int size()
    {
#ifdef USE_MPI
        int s = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &s);
        return s;
#else
        return 1;
#endif
    }

    inline void barrier()
    {
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    // Abort all ranks (used for fatal-error paths so MPI runs don't deadlock).
    [[noreturn]] inline void abort_all(int code = 1)
    {
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD, code);
#endif
        std::exit(code);
    }

    // Collective reductions. In serial builds these are identity operations.
    inline double sum(double local)
    {
#ifdef USE_MPI
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global;
#else
        return local;
#endif
    }

    inline double max(double local)
    {
#ifdef USE_MPI
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global;
#else
        return local;
#endif
    }

    inline bool any(bool local)
    {
#ifdef USE_MPI
        int lo = local ? 1 : 0, go = 0;
        MPI_Allreduce(&lo, &go, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        return go != 0;
#else
        return local;
#endif
    }

    // Broadcast a small fixed-size buffer of doubles from rank 0 to all ranks.
    inline void bcast_doubles(double* buf, int n, int root = 0)
    {
#ifdef USE_MPI
        MPI_Bcast(buf, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
#else
        (void)buf; (void)n; (void)root;
#endif
    }

    inline void bcast_ints(int* buf, int n, int root = 0)
    {
#ifdef USE_MPI
        MPI_Bcast(buf, n, MPI_INT, root, MPI_COMM_WORLD);
#else
        (void)buf; (void)n; (void)root;
#endif
    }
}

#endif // MPI_HELPERS_H_
