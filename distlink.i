%module distlink
%feature("autodoc", "0");

%{
#define SWIG_FILE_WITH_INIT
#include "distlink.h"
%}

// function detect_suitable_options
template<typename realfp>
void detect_suitable_options(realfp& max_root_error,
                             realfp& min_root_error,
                             realfp& max_anom_error);
%template(detect_suitable_options) detect_suitable_options<double>;

// class COrbitData
template<typename realfp>
class COrbitData
{
  public:
    COrbitData();
    COrbitData(realfp a_, realfp e_, realfp i_, realfp w_, realfp Om_);

    void set_data(realfp a_, realfp e_, realfp i_, realfp w_, realfp Om_);
    void get_data(realfp& a_, realfp& e_, realfp& i_, realfp& w_, realfp& Om_);

    realfp get_a();
    realfp get_e();
    realfp get_i();
    realfp get_w();
    realfp get_Om();

    const realfp* vectorP();
    const realfp* vectorQ();
    void get_vectors(realfp P_[3], realfp Q_[3]);
};
%template(COrbitData) COrbitData<double>;

// function test_peri_apo
template<typename realfp>
bool test_peri_apo(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2, realfp limit);
%template(test_peri_apo) test_peri_apo<double>;

// struct SMOIDResult
template<typename realfp>
struct SMOIDResult
{
 bool good;                   //true if the result is reliable
 realfp distance;             //minimum distance between orbits
 realfp distance_error;       //numeric uncertainty of the distance
 realfp u1;                   //eccentric anomaly on the first orbit
 realfp u1_error;             //its numeric uncertainty
 realfp u2;                   //eccentric anomaly on the second orbit
 realfp u2_error;             //its numeric uncertainty
 unsigned short root_count;   //number of real roots of the function g(u1)
 realfp min_delta;            //the minimum quantity delta among all non-real roots - see description in the paper
 unsigned long iter_count;    //number of Newtonian iterations of g(u), sum for all 16 roots
 unsigned long iter_count_2D; //number of Newtonian 2D iterations of rho(u,u')
 realfp time;                 //CPU time used (in seconds), actually appears unreliable
                              //in modern hardware (millisecond accuracy appears not enough)
 SMOIDResult();
};
%template(SMOIDResult) SMOIDResult<double>;

// struct SLCResult
template<typename realfp>
struct SLCResult
{
 realfp I;                    //mutual inclination
 realfp l;                    //the first or the third linking coefficient; it depends on I
 realfp lmod;                 //modified first linking coefficient, |lmod| is always smaller than |l|, and is a good upper limit on MOID^2
 realfp l2;                   //the continuous second linking coefficient; it is calculated for every pair of orbits,
                              //regardless of whether the mutual inclination I is small (however small) or not
 SLCResult();
};
%template(SLCResult) SLCResult<double>;

// function LC
template<typename realfp>
SLCResult<realfp> LC(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2, realfp min_mut_incl);
%template(LC) LC<double>;

// function MOID_fast
template<typename realfp>
SMOIDResult<realfp> MOID_fast(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2,
                              realfp max_root_error, realfp min_root_error,
                              realfp nu=static_cast<realfp>(1));
%template(MOID_fast) MOID_fast<double>;

// function MOID_direct_search
template<typename realfp>
SMOIDResult<realfp> MOID_direct_search(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2,
                                       const unsigned int densities[],
                                       realfp max_dist_error, realfp max_anom_error);
%template(MOID_direct_search) MOID_direct_search<double>;

// function restrict_search_range
template<typename realfp>
int restrict_search_range(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2,
                          std::pair<realfp,realfp>& u1, std::pair<realfp,realfp>& u1_);
%template(restrict_search_range) restrict_search_range<double>;
