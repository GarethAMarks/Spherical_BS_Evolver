#ifndef CCZ4_CPP_
#define CCZ4_CPP_

/*#include "EvolutionVariables.h"
#include "DimensionMacros.h"
#include "mathutils.h"
#include "CCZ4.h"
#include <iomanip>
#include<algorithm>
#include<filesystem>
#include <complex.h>

//additional data (just theta) and functions needed for CCZ4 on a given slice.
class CCZ4Slice
{
    //public vars are BSSN evolution/gauge variables; privates are auxiliary quantities used in evolution/diagnostic equations
    private:
    double refinement_level;

    public:
        std::vector<BSSNState> states; //array of states on a BSSN time slice

        std::vector<int> refinement_points; //points at which refinement should be halved; to be copied from spacetime version at start of each step

        void read_BS_data(BosonStar& boson_star,int BS_resolution_factor = 1., bool isotropic = 1);
        void read_checkpoint(int time, int n_gridpoints);
        void write_CCZ4( std::string file_name = "CCZ4Data.dat");



        double R;
        double d_z(bssn_var var, int index, int order); //takes z derivative of a BSSN variable on the current slice centred on given index
        double d_zz(bssn_var var, int index);  //second z-derivative, basically syntactic sugar for d_z(var, index, 2)


    friend class Spacetime; //Spacetime should be able to access private CCZ4Slice members
    friend CCZ4Slice operator+(const CCZ4Slice& slice1, const CCZ4Slice& slice2);
    friend CCZ4Slice operator*(double c, const CCZ4Slice slice);

*/
#endif /*CCZ4*/
