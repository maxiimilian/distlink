/*
Copyright (c) 2018-2020 R.V. Baluev and D.V. Mikryukov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "distlink.h"

#include <cstdlib>
#include <cmath>
#include <limits>
#include <complex>
#include <ctime>
#include <algorithm>
#include <utility>

#define DIM 21      // must be odd greater than DEG
#define DEG 16

using namespace std;

// square of x
template<typename T> inline T sqr(T x) {return x*x;}
template<typename T> inline short sign(T x) {return (x>0) - (x<0);}

// asin with range check
template<typename realfp> inline realfp safe_asin(realfp x) {return asin(max<realfp>(min<realfp>(x,1),-1));}
// sqrt with range check
template<typename realfp> inline realfp safe_sqrt(realfp x) {return sqrt(max<realfp>(x,0));}
// acosh with range check
template<typename realfp> inline realfp safe_acosh(realfp x) {return acosh(max<realfp>(x,1));}

template<typename realfp> inline realfp atan2h(realfp y, realfp x) {return (log(fabs(x+y))-log(fabs(x-y)))/2;} // real part

template<typename realfp> inline realfp atan2_smart(realfp y, realfp x, bool trig) {return trig ? atan2(y,x) : atan2h(y,x);}

template<typename realfp>
inline void cross_product(const realfp a[3], const realfp b[3], realfp res[3])
{
 res[0] = a[1]*b[2]-a[2]*b[1];
 res[1] = a[2]*b[0]-a[0]*b[2];
 res[2] = a[0]*b[1]-a[1]*b[0];
}

template<typename realfp>
inline realfp dot_product(const realfp a[3], const realfp b[3]) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

template<typename realfp>
inline realfp norm(const realfp a[3]) {return dot_product(a,a);}

//Some extra definitions for the template type "complex" from <complex>
template<typename T> inline complex<T> operator+(const complex<T>& x, int y) {return x+static_cast<T>(y);}
template<typename T> inline complex<T> operator+(int x, const complex<T>& y) {return static_cast<T>(x)+y;}
template<typename T> inline complex<T> operator-(const complex<T>& x, int y) {return x-static_cast<T>(y);}
template<typename T> inline complex<T> operator-(int x, const complex<T>& y) {return static_cast<T>(x)-y;}
template<typename T> inline complex<T> operator*(const complex<T>& x, int y) {return x*static_cast<T>(y);}
template<typename T> inline complex<T> operator*(int x, const complex<T>& y) {return static_cast<T>(x)*y;}
template<typename T> inline complex<T> operator/(const complex<T>& x, int y) {return x/static_cast<T>(y);}
template<typename T> inline complex<T> operator/(int x, const complex<T>& y) {return static_cast<T>(x)/y;}
template<typename T> inline complex<T> inverse(const complex<T>& x) {return conj(x)/norm(x);}

template<typename T>
inline T intpow(const T& x, int n)
{
 if(n==0) return 1;
 if(n==1) return x;
 if(n<0) return intpow(T(1)/x,-n);
 if(n>1) return n%2==0 ? sqr(intpow(x,n/2)) : x*intpow(x,n-1);
}

template<typename realfp>
inline realfp pi()
{
 static const realfp _pi = acos(static_cast<realfp>(-1));
 return _pi;
}

template<typename realfp>
inline realfp angle_wrap(realfp x)
{
 static const realfp circ = 2*pi<realfp>();
 x = fmod(x+pi<realfp>(),circ);
 return x<0 ? x+pi<realfp>() : x-pi<realfp>();
}

//This class is for initial setup.
template<typename realfp, unsigned short dimension>
class CExps
{
 public:
  complex<realfp> expvals[dimension];

  CExps(){
   for(int i=0; i<dimension; i++)
    expvals[i] = polar<realfp>(1,2*i*pi<realfp>()/dimension);}
};

template<int dimension>
class CInitializer
{
 public:
  CExps<float,dimension> expsf;
  CExps<double,dimension> expsd;
  CExps<long double,dimension> expsld;

  CInitializer() {srand(time(0));}
};

const CInitializer<DIM> Init;

template<typename realfp> inline const complex<realfp>* exps();

template<> inline const complex<float>* exps<float>() {return Init.expsf.expvals;}
template<> inline const complex<double>* exps<double>() {return Init.expsd.expvals;}
template<> inline const complex<long double>* exps<long double>() {return Init.expsld.expvals;}

template<typename realfp> inline realfp ierr() {return numeric_limits<realfp>::epsilon();}

template<typename realfp> inline realfp random() {return static_cast<realfp>(rand())/RAND_MAX;}

template<typename realfp>
void detect_suitable_options(realfp& max_root_error,
                             realfp& min_root_error,
                             realfp& max_anom_error)
{
 max_root_error = sqrt(ierr<realfp>());
 min_root_error = ierr<realfp>()*2;
 max_anom_error = ierr<realfp>()*1000;
}

// Implementation of COrbitData.
template<typename realfp>
COrbitData<realfp>::COrbitData(): a(1.), e(0.), i(0.), w(0.), Om(0.) {
    P[0]=1.; P[1]=0.; P[2]=0.;
    Q[0]=0.; Q[1]=1.; Q[2]=0.;}

template<typename realfp>
COrbitData<realfp>::COrbitData(realfp a_, realfp e_, realfp i_, realfp w_, realfp Om_):
    a(a_), e(e_), i(i_), w(w_), Om(Om_) {set_vectors();}

template<typename realfp>
template<typename T>
COrbitData<realfp>::COrbitData(const COrbitData<T>& other):
    a(other.a), e(other.e), i(other.i), w(other.w), Om(other.Om) {set_vectors();}

//Components of P and Q vectors.
template<typename realfp>
void COrbitData<realfp>::set_vectors()
{
 P[0]= cos(w)*cos(Om)-cos(i)*sin(w)*sin(Om);
 P[1]= cos(w)*sin(Om)+cos(i)*sin(w)*cos(Om);
 P[2]= sin(i)*sin(w);
 Q[0]=-sin(w)*cos(Om)-cos(i)*cos(w)*sin(Om);
 Q[1]= cos(i)*cos(w)*cos(Om)-sin(w)*sin(Om);
 Q[2]= sin(i)*cos(w);
}

//This function evaluates the value of the algebraic polynomial.
template<typename realfp>
complex<realfp> polynomial_val(int n, const complex<realfp> c[], const complex<realfp>& z, bool forward=true)
{
 complex<realfp> R;
 if(forward)
 {
  R=c[n];
  for(int i=n-1; i>=0; i--)
   R=R*z+c[i];
 }
 else
 {
  R=c[0];
  for(int i=1; i<=n; i++)
   R=R*z+c[i];
 }
 return R;
}

//This function evaluates simultaneously the value of the algebraic polynomial and of its derivative.
template<typename realfp>
void polynomial_valder(int n, const complex<realfp> c[], const complex<realfp>& z,
                       complex<realfp>& val, complex<realfp>& der, bool forward)
{
 if(forward)
 {
  val=c[n];
  der=c[n]*n;
  for(int i=n-1;i>1;i--)
  {
   val=val*z+c[i];
   der=der*z+c[i]*i;
  }
  if(n>1)
  {
   val=val*z+c[1];
   der=der*z+c[1];
  }
  if(n>0)
   val=val*z+c[0];
 }
 else
 {
  val=c[0];
  der=c[0]*n;
  for(int i=1;i<n-1;i++)
  {
   val=val*z+c[i];
   der=der*z+c[i]*(n-i);
  }
  if(n>1)
  {
   val=val*z+c[n-1];
   der=der*z+c[n-1];
  }
  if(n>0)
   val=val*z+c[n];
 }
}

//This function evaluates simultaneously the value of the algebraic polynomial and of its first and second derivatives.
template<typename realfp>
void polynomial_valder(int n, const complex<realfp> c[], const complex<realfp>& z,
                       complex<realfp>& val, complex<realfp>& der, complex<realfp>& der2, bool forward)
{
 if(forward)
 {
  val =c[n];
  der =c[n]*n;
  der2=c[n]*(n*(n-1));
  for(int i=n-1; i>2; i--)
  {
   val =val*z+c[i];
   der =der*z+c[i]*i;
   der2=der2*z+c[i]*(i*(i-1));
  }
  if(n>2)
  {
   val=val*z+c[2];
   der=der*z+c[2]*2;
   der2=der2*z+c[2]*2;
  }
  if(n>1)
  {
   val=val*z+c[1];
   der=der*z+c[1];
  }
  if(n>0)
   val=val*z+c[0];
 }
 else
 {
  val =c[0];
  der =c[0]*n;
  der2=c[0]*(n*(n-1));
  for(int i=1; i<n-2; i++)
  {
   val =val*z+c[i];
   der =der*z+c[i]*(n-i);
   der2=der2*z+c[i]*((n-i)*(n-i-1));
  }
  if(n>2)
  {
   val=val*z+c[n-2];
   der=der*z+c[n-2]*2;
   der2=der2*z+c[n-2]*2;
  }
  if(n>1)
  {
   val=val*z+c[n-1];
   der=der*z+c[n-1];
  }
  if(n>0)
   val=val*z+c[n];
 }
}

template<typename realfp>
inline complex<realfp> ratio(int n, const complex<realfp> c[], const complex<realfp>& z, bool forward)
{
 complex<realfp> f, fd;
 polynomial_valder(n, c, z, f, fd, forward);
 return f/fd;
}

template<typename realfp>
inline complex<realfp> ratio2(const complex<realfp> c[4], const complex<realfp>& z, bool forward)
{
 return forward ? ((c[2]*z+c[1])*z+c[0])/(2*c[2]*z+c[1]) :
                  ((c[0]*z+c[1])*z+c[2])/(2*c[0]*z+c[1]);
}

template<typename realfp>
inline complex<realfp> ratio4(const complex<realfp> c[4], const complex<realfp>& z, bool forward)
{
 return forward ? ((((c[4]*z+c[3])*z+c[2])*z+c[1])*z+c[0])/(((4*c[4]*z+3*c[3])*z+2*c[2])*z+c[1]) :
                  ((((c[0]*z+c[1])*z+c[2])*z+c[3])*z+c[4])/(((4*c[0]*z+3*c[1])*z+2*c[2])*z+c[3]);
}

/*
  This function evaluates one complex root of the algebraic polynomial with complex coefficients.
  Input:
    n       - degree of the polynomial;
    c       - array of complex coefficients;
    z       - initial approximation of the root;
    maxeps  - maximum relative error allowed for determination of the root;
    mineps  - desirable relative error of the root; it may be set to zero;
    maxcount- maximum number of iterations in a single run;
    restartcount - maximum number of runs.
  Output:
    z       - final approximation of the root.
  Return value: the full number of iterations.
*/
template<typename realfp>
unsigned long Newton(int n, complex<realfp> c[], complex<realfp> &z, realfp maxeps, realfp mineps,
                     unsigned long maxcount, unsigned long restartcount)
{
 complex<realfp> dz;
 maxeps=fabs(maxeps);
 mineps=fabs(mineps);
 complex<realfp> z1=z,z2=z1;

 //correcting inadequate values of the least precision required
 if(!(maxeps>=2*mineps)) maxeps=2*mineps;
 if(!(maxeps>=2*ierr<realfp>())) maxeps=2*ierr<realfp>();

 unsigned long iter_cnt=0;
 realfp err2=0,nrm;
 bool inner=true;
 do
 {
  iter_cnt++;
  //if the Newtonian algorithm hanged in maxcount steps,
  //change the starting approximation and restart it:
  const bool cycled = (iter_cnt>2 && norm(z-z2)<norm(z-z1)*1e-6) || iter_cnt%maxcount == 0;
  complex<realfp> rand_mlt;
  if(cycled) rand_mlt = polar<realfp>((2*random<realfp>()+1)/3,random<realfp>()*2*pi<realfp>());

  nrm=norm(z);
  if(inner) {if(nrm>5) inner=false;}
  else {if(2*nrm<1) inner=true;}
  z2=z1;z1=z;
  if(inner)
  {
   dz=ratio(n,c,z,true);
   z-=cycled?dz*rand_mlt:dz;
  }
  else
  {
   nrm=1/nrm;
   dz=ratio(n,c,conj(z)*nrm,false); // actually, dw
   z/=(1-z*(cycled?dz*rand_mlt:dz));
  }
  err2=norm(dz);
 }
 //require relative precision of at least maxeps...
 while(err2>sqr(maxeps)*nrm && iter_cnt<maxcount*restartcount);

 if(err2>sqr(mineps)*nrm)
 {
  unsigned long _iter_cnt=0;
  realfp derr2;
  do
  {
   nrm=norm(z);
   if(inner)
   {
    dz=ratio(n,c,z,true);
    z-=dz;
   }
   else
   {
    nrm=1/nrm;
    dz=ratio(n,c,conj(z)*nrm,false);
    z/=(1-z*dz);
   }
   const realfp err2_=norm(dz);
   derr2=err2-err2_;
   err2=err2_;
   _iter_cnt++;
  }
  // ...and even more, until the desirable or the maximum possible precision is reached
  while(derr2>ierr<realfp>()*err2 && err2>sqr(mineps)*nrm && _iter_cnt<=maxcount);
  iter_cnt += _iter_cnt;
 }

 return iter_cnt;
}

/*
  This function uses Ruffini-Horner algorithm for decrementing the degree of the algebraic
  polynomial (if one complex root was obtained before calling this function).

  For numerical stability, it extracts:
  factor (z-root)     out of P(z)   = sum_{k=0..n}{c_k z^k} if |root|<1
  factor (1/z-1/root) out of Q(1/z) = sum_{k=0..n}{c_{n-k}*(1/z)^k} otherwise.

  In the first case, it does not spend time for shifting the array of coefficients
  according to c[k]:=c[k+1]. The new set of coefficients starts from c[1] and ends by
  c[n-1]. The element c[0] should be disregarded.

  In the second case, the new set of coefficients starts from c[0] and ends by c[n-2]. The
  element c[n-1] should be disregarded.

  Input:
    n    - degree of the polynomial;
    c    - array of complex coefficients of the polynomial;
    root - the root of the polynomial.
  Output:
    c    - array of the complex coefficients of the new polynomial.
  Return value: pointer to the first coefficient of the new polynomial.
*/
template<typename realfp>
inline complex<realfp>* extract_linear_factor(int n, complex<realfp> c[], const complex<realfp> &root)
{
 if(norm(root)<=1)
 {
  for(int i=n-1;i>0;i--) c[i]+=c[i+1]*root;
  return c+1;
 }
 else
 {
  const complex<realfp> iroot=inverse(root);
  for(int i=1; i<n; i++) c[i]+=c[i-1]*iroot;
  return c;
 }
}

/*
  This function solves the quadratic polynomial equation
*/
template<typename realfp>
unsigned int solve_quadratic(const complex<realfp> c[3], complex<realfp> x[2])
{
 const complex<realfp> D=sqrt(c[1]*c[1]-4*c[2]*c[0]);
 const complex<realfp> tmp=(real(c[1])*real(D)+imag(c[1])*imag(D)>=0 ? -c[1]-D : -c[1]+D);
 x[0] = tmp/(c[2]*2);
 x[1] = (c[0]*2)/tmp;
 if(norm(x[0])<=1) x[0]-=ratio2(c,x[0],true);
 else              x[0]/=(1-x[0]*ratio2(c,inverse(x[0]),false));
 if(norm(x[1])<=1) x[1]-=ratio2(c,x[1],true);
 else              x[1]/=(1-x[1]*ratio2(c,inverse(x[1]),false));
 return 2;
}

/*
  This function solves the quartic polynomial equation
*/
template<typename realfp>
unsigned int solve_quartic(const complex<realfp> c[5], complex<realfp> x[4])
{
 const complex<realfp> r  = 4*c[4];
 const complex<realfp> p  = 2*r*c[2]-3*sqr(c[3]);
 const complex<realfp> q  = 2*c[3]*(sqr(c[3])-r*c[2])+sqr(r)*c[1];
 const complex<realfp> D0 = sqr(c[2])-3*c[3]*c[1]+3*c[0]*r;
 const complex<realfp> D1 = 2*c[2]*sqr(c[2])-9*c[3]*c[2]*c[1]+27*sqr(c[3])*c[0]+27*sqr(c[1])*c[4]-18*r*c[2]*c[0];
 complex<realfp> tmp= sqrt(sqr(D1)-4*D0*sqr(D0));
 const complex<realfp> Qc = (real(D1)*real(tmp)+imag(D1)*imag(tmp)>=0 ? D1+tmp : D1-tmp)/2;
 const realfp Qabs=cbrt(abs(Qc));
 const realfp Qarg=arg(Qc)/3;
 const complex<realfp> Q1=polar<realfp>(Qabs,Qarg);
 static const realfp twopi3 = 2*pi<realfp>()/3;
 const complex<realfp> Q2=polar<realfp>(Qabs,Qarg+twopi3);
 const complex<realfp> Q3=polar<realfp>(Qabs,Qarg-twopi3);
 const complex<realfp> vals[3] = {(Q1+D0/Q1)*r-p, (Q2+D0/Q2)*r-p, (Q3+D0/Q3)*r-p};
 const realfp nrms[3] = {norm(vals[0]), norm(vals[1]), norm(vals[2])};
 tmp=sqrt(vals[nrms[0]>max(nrms[1],nrms[2])?0:(nrms[1]>nrms[2]?1:2)]/3);
 const complex<realfp> S = (real(c[4])>=0 ? tmp : -tmp);
 tmp = sqrt( q/S-sqr(S)-p);
 x[0] = (-c[3]-S+tmp)/r;
 x[1] = (-c[3]-S-tmp)/r;
 tmp = sqrt(-q/S-sqr(S)-p);
 x[2] = (-c[3]+S+tmp)/r;
 x[3] = (-c[3]+S-tmp)/r;
 if(norm(x[0])<=1) x[0]-=ratio4(c,x[0],true);
 else              x[0]/=(1-x[0]*ratio4(c,inverse(x[0]),false));
 if(norm(x[1])<=1) x[1]-=ratio4(c,x[1],true);
 else              x[1]/=(1-x[1]*ratio4(c,inverse(x[1]),false));
 if(norm(x[2])<=1) x[2]-=ratio4(c,x[2],true);
 else              x[2]/=(1-x[2]*ratio4(c,inverse(x[2]),false));
 if(norm(x[3])<=1) x[3]-=ratio4(c,x[3],true);
 else              x[3]/=(1-x[3]*ratio4(c,inverse(x[3]),false));
 return 4;
}

/*
  This function evaluates all roots of algebraic polinomials of degree n.

  P(z) = sum_{k=0..n}{c_k*z^k} = 0 <===> Q(1/z) = sum_{k=0..n}{c_{n-k}*(1/z)^k} = 0.

  The algebraic polynomial is assumed to be originated from a trigonometric polynomial
  with real coefficients of degree n/2.

  Input:
    n      - degree of the algebraic polynomial;
    c      - array of complex coefficients of the polynomial;
    minerr - desirable precision of the roots;
    maxerr - maximum error of the roots allowed.
    trigonometric - boolean indicator

  Output:
    roots  - array of the roots.
*/
template<typename realfp>
unsigned long polynomial_roots(int n, complex<realfp> c[], complex<realfp> roots[],
                               realfp minerr, realfp maxerr, bool trigonometric)
{
 complex<realfp>* c_ = c; // pointer to the first coefficient
 int cnt = 0;
 unsigned long iter_cnt = 0;
 //extracting roots until a square polynomial is obtained
 for(int i=0, degree=n; degree>4; i++)
 {
  iter_cnt += Newton(degree, c_, roots[i], maxerr, minerr, 300, 10); //Newtonian root search
  c_ = extract_linear_factor(degree, c_, roots[i]);                  //extracting the root
  degree--;
  if(degree>4 && roots[i+1]==complex<realfp>(0,0))
   if(trigonometric)
   //if z is a root of an algebraic polynomial
   //produced by a trigonometric polynomial with real coefficients
   //then 1/z^* is also its root;
   //therefore, the starting approximation for the next root is 1/z^*
    roots[i+1]=roots[i]/norm(roots[i]);
   else
   //if z is a root of an algebraic polynomial with real coefficients
   //then z^* is also its root;
   //therefore, the starting approximation for the next root is z^*
    if(fabs(imag(roots[i]))>fabs(real(roots[i]))*1e-3)  // additionally check that it is not a real-valued root
     roots[i+1]=conj(roots[i]);
    else roots[i+1]=complex<realfp>(imag(roots[i]),real(roots[i]));
 }
 //solving the final square equation
 iter_cnt += solve_quartic(c_,roots+n-4);
 return iter_cnt;
}

/*
  This function evaluates the uncertainty of the root of the algebraic polynomial.

  Input:
    n     - degree of the polynomial;
    c     - array of the complex coefficients of the polynomial;
    cerr  - estimated error of the coefficients of the polynomial;
    z     - value of the root.
  Return value:
    relative uncertainty of the root.
*/
template<typename realfp>
realfp root_error(int n, const complex<realfp> c[], realfp cerr, const complex<realfp>& z)
{
 realfp m=norm(z);
 const bool forward = m<=1;
 if(!forward) m=1/m;
 //evaluating expected error of calculating the value of the polynomial:
 realfp ferr2=1;
 for(int i=1; i<=n; i++)
  {ferr2*=m;ferr2++;}
 ferr2*=sqr(cerr);
 //using quadratic approximation to estimate the uncertainty of the root:
 complex<realfp> f, fd, fd2;
 polynomial_valder(n, c, forward ? z : inverse(z), f, fd, fd2, forward);
 const complex<realfp> D=sqrt(fd*fd-2*f*fd2);
 const realfp d2=4*norm(f)/max(norm(fd+D),norm(fd-D));
 return sqrt((d2+ferr2/norm(D))/m);
}

/*
  This function evaluates the maximum root bound.

  Input:
    n     - degree of the polynomial;
    c     - array of the complex coefficients of the polynomial;
  Return value:
    maximum bound on |z_k|.
*/
template<typename realfp>
realfp max_root_bound(int n, const complex<realfp> c[], bool forward)
{
 const realfp nch=norm(c[forward?n:0]);
 if(n<1) return numeric_limits<realfp>::infinity();
 realfp res=norm(c[forward?n-1:1])/nch;
 for(int i=2; i<=n; i++)
 {
  const realfp tmp=pow(norm(c[forward?n-i:i])/nch,static_cast<realfp>(1)/i);
  if(res<tmp) res=tmp;
 }
 return sqrt(res);
}

//This structure contains necessary pre-calculated data for a pair of orbits.
template<typename realfp>
struct SAuxData
{
 //these data are for calculating MOID
 realfp e1, e2, a1, a2;
 realfp alpha1, alpha2;
 realfp K;
 realfp Pp, Ps, Sp, Ss;

 //these data are for calculating linking coefficients
 realfp p1, p2, w1, w2, I, abs_w;
 realfp P1w, P2w, Q1w, Q2w;

 const realfp* P1;
 const realfp* P2;
 const realfp* Q1;
 const realfp* Q2;

 SAuxData(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2)
 {
  e1 = O1.get_e();   e2 = O2.get_e();
  a1 = O1.get_a();   a2 = O2.get_a();
  w1 = O1.get_w();   w2 = O2.get_w();
  if(e1<0) {e1=-e1;w1=w1+pi<realfp>();}
  if(e2<0) {e2=-e2;w2=w2+pi<realfp>();}
  a1 = fabs(a1);       a2 = fabs(a2);
  if(!(e1<=1)) a1=-a1;
  if(!(e2<=1)) a2=-a2;
  p1 = a1*(1-e1*e1); p2 = a2*(1-e2*e2);
  alpha1 = a1/a2;    alpha2 = a2/a1;
  K = alpha2*e2*e2;
  const realfp eta1 = sqrt(fabs(1-e1*e1)); const realfp eta2 = sqrt(fabs(1-e2*e2));
  const realfp   i1 = O1.get_i();    const realfp   i2 = O2.get_i();
  const realfp  Om1 = O1.get_Om();   const realfp  Om2 = O2.get_Om();
  const realfp c1=cos(i1);  const realfp s1=sin(i1);
  const realfp c2=cos(i2);  const realfp s2=sin(i2);
  const realfp w[3] = {c1*s2*cos(Om2)-s1*c2*cos(Om1),
                       c1*s2*sin(Om2)-s1*c2*sin(Om1),
                       s1*s2*sin(Om2-Om1)};
  abs_w = sqrt(norm(w));
  const realfp cosI = c1*c2+s1*s2*cos(Om2-Om1);
  I = cosI>0 ? safe_asin(abs_w) : pi<realfp>()-safe_asin(abs_w);  //radians

  //second small letter in Pp, Ps, Sp, Ss refers to the orbit O2
  P1=O1.vectorP(); P2=O2.vectorP();
  Q1=O1.vectorQ(); Q2=O2.vectorQ();
  Pp = dot_product(P1,P2);
  Ps = dot_product(P1,Q2)*eta2;
  Sp = dot_product(Q1,P2)*eta1;
  Ss = dot_product(Q1,Q2)*eta1*eta2;

  P1w = dot_product(P1,w); P2w = dot_product(P2,w);
  Q1w = dot_product(Q1,w); Q2w = dot_product(Q2,w);

//  P1w = cos(i1)*sin(i2)*cos(w1)*cos(Om1-Om2)-sin(i1)*cos(i2)*cos(w1)+sin(i2)*sin(w1)*sin(Om2-Om1);
//  P2w =-cos(i2)*sin(i1)*cos(w2)*cos(Om1-Om2)+sin(i2)*cos(i1)*cos(w2)+sin(i1)*sin(w2)*sin(Om2-Om1);

//  Q1w =-cos(i1)*sin(i2)*sin(w1)*cos(Om1-Om2)+sin(i1)*cos(i2)*sin(w1)+sin(i2)*cos(w1)*sin(Om2-Om1);
//  Q2w = cos(i2)*sin(i1)*sin(w2)*cos(Om1-Om2)-sin(i2)*cos(i1)*sin(w2)+sin(i1)*cos(w2)*sin(Om2-Om1);
 }
};

//This is the implementation of the function g(u). It is not used anywhere.
template<typename realfp>
realfp func_g_ee(realfp u, const SAuxData<realfp>& data)
{
 const realfp x = cos(u);
 const realfp y = sin(u);
 const realfp A = data.Ps*y-data.Ss*x;
 const realfp B = data.Pp*y-data.Sp*x;
 const realfp C = data.e2*B-data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp M = data.Sp*y+data.Pp*(x-data.e1)+data.alpha2*data.e2;
 const realfp N =-data.Ss*y-data.Ps*(x-data.e1);
 const realfp C2 = C*C;
 const realfp A2 = A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*Ac*Bc+2*data.K*C*(N*A*Ac+M*B*Bc)-(A2+B2)*(N*N*Ac+M*M*Bc-2*N*M*A*B);
}

template<typename realfp>
realfp func_g_he(realfp u, const SAuxData<realfp>& data)
{
 const realfp eu= exp(u);
 const realfp x = (eu+1/eu)/2;
 const realfp y = (eu-1/eu)/2;
 const realfp A = data.Ps*y-data.Ss*x;
 const realfp B = data.Pp*y-data.Sp*x;
 const realfp C = data.e2*B-data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp M =-data.Sp*y+data.Pp*(x-data.e1)+data.alpha2*data.e2;
 const realfp N = data.Ss*y-data.Ps*(x-data.e1);
 const realfp C2 = C*C;
 const realfp A2 = A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*Ac*Bc+2*data.K*C*(N*A*Ac+M*B*Bc)-(A2+B2)*(N*N*Ac+M*M*Bc-2*N*M*A*B);
}

template<typename realfp>
realfp func_g_eh(realfp u, const SAuxData<realfp>& data)
{
 const realfp x = cos(u);
 const realfp y = sin(u);
 const realfp A = data.Ps*y-data.Ss*x;
 const realfp B = data.Pp*y-data.Sp*x;
 const realfp C = data.e2*B-data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp M = data.Sp*y+data.Pp*(x-data.e1)+data.alpha2*data.e2;
 const realfp N =-data.Ss*y-data.Ps*(x-data.e1);
 const realfp C2 = C*C;
 const realfp A2 =-A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*Ac*Bc+2*data.K*C*(-N*A*Ac+M*B*Bc)-(A2+B2)*(-N*N*Ac+M*M*Bc+2*N*M*A*B);
}

template<typename realfp>
realfp func_g_hh(realfp u, const SAuxData<realfp>& data)
{
 const realfp eu= exp(u);
 const realfp x = (eu+1/eu)/2;
 const realfp y = (eu-1/eu)/2;
 const realfp A = data.Ps*y-data.Ss*x;
 const realfp B = data.Pp*y-data.Sp*x;
 const realfp C = data.e2*B-data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp M =-data.Sp*y+data.Pp*(x-data.e1)+data.alpha2*data.e2;
 const realfp N = data.Ss*y-data.Ps*(x-data.e1);
 const realfp C2 = C*C;
 const realfp A2 =-A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*Ac*Bc+2*data.K*C*(-N*A*Ac+M*B*Bc)-(A2+B2)*(-N*N*Ac+M*M*Bc+2*N*M*A*B);
}

//This is a fast algorithm for evaluating g(u) at u = 2*M_PI*j/N. This algorithm does not
//use calls of trigonometric functions, because these data were pre-calculated.
template<typename realfp>
realfp func_g_ee_fast(int j, const SAuxData<realfp>& data)
{
 const realfp x = real(exps<realfp>()[j]); //cos(u)
 const realfp y = imag(exps<realfp>()[j]); //sin(u)
 const realfp A = data.Ps*y - data.Ss*x;
 const realfp B = data.Pp*y - data.Sp*x;
 const realfp C = data.e2*B - data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp x_= x - data.e1;
 const realfp M = data.Sp*y + data.Pp*x_+data.alpha2*data.e2;
 const realfp N =-data.Ss*y - data.Ps*x_;
 const realfp C2 = C*C;
 const realfp A2 = A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*(Ac*Bc) + 2*data.K*C*(N*A*Ac+M*B*Bc) - (A2+B2)*(N*N*Ac+M*M*Bc-2*N*M*A*B);
}

// computes complex-valued g for imaginary hyperbolic argument
template<typename realfp>
complex<realfp> func_g_he_fast(int j, const SAuxData<realfp>& data)
{
 const realfp x = real(exps<realfp>()[j]); //cos(u)
 const realfp y = imag(exps<realfp>()[j]); //sin(u)
 const complex<realfp> A(data.Ps*y, -data.Ss*x);
 const complex<realfp> B(data.Pp*y, -data.Sp*x);
 const complex<realfp> C = data.e2*B - data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp x_= x - data.e1;
 const complex<realfp> M( data.Pp*x_+data.alpha2*data.e2, data.Sp*y);
 const complex<realfp> N(-data.Ps*x_,                    -data.Ss*y);
 const complex<realfp> C2 = C*C;
 const complex<realfp> A2 = A*A;
 const complex<realfp> B2 = B*B;
 const complex<realfp> Ac = A2-C2;
 const complex<realfp> Bc = B2-C2;
 return data.K*data.K*(Ac*Bc) + 2*data.K*C*(N*A*Ac+M*B*Bc) - (A2+B2)*(N*N*Ac+M*M*Bc-2*N*M*A*B);
}

template<typename realfp>
realfp func_g_eh_fast(int j, const SAuxData<realfp>& data)
{
 const realfp x = real(exps<realfp>()[j]); //cos(u)
 const realfp y = imag(exps<realfp>()[j]); //sin(u)
 const realfp A = data.Ps*y - data.Ss*x;
 const realfp B = data.Pp*y - data.Sp*x;
 const realfp C = data.e2*B - data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp x_=         x - data.e1;
 const realfp M = data.Sp*y + data.Pp*x_+data.alpha2*data.e2;
 const realfp N =-data.Ss*y - data.Ps*x_;
 const realfp C2 = C*C;
 const realfp A2 =-A*A;
 const realfp B2 = B*B;
 const realfp Ac = A2-C2;
 const realfp Bc = B2-C2;
 return data.K*data.K*(Ac*Bc) + 2*data.K*C*(-N*A*Ac+M*B*Bc) - (A2+B2)*(-N*N*Ac+M*M*Bc+2*N*M*A*B);
}

// computes complex-valued g for imaginary hyperbolic argument
template<typename realfp>
complex<realfp> func_g_hh_fast(int j, const SAuxData<realfp>& data)
{
 const realfp x = real(exps<realfp>()[j]); //cos(u)
 const realfp y = imag(exps<realfp>()[j]); //sin(u)
 const complex<realfp> A(data.Ps*y, -data.Ss*x);
 const complex<realfp> B(data.Pp*y, -data.Sp*x);
 const complex<realfp> C = data.e2*B - data.alpha1*data.e1*y*(1-data.e1*x);
 const realfp x_= x - data.e1;
 const complex<realfp> M( data.Pp*x_+data.alpha2*data.e2, data.Sp*y);
 const complex<realfp> N(-data.Ps*x_,                    -data.Ss*y);
 const complex<realfp> C2 = C*C;
 const complex<realfp> A2 =-A*A;
 const complex<realfp> B2 = B*B;
 const complex<realfp> Ac = A2-C2;
 const complex<realfp> Bc = B2-C2;
 return data.K*data.K*(Ac*Bc) + 2*data.K*C*(-N*A*Ac+M*B*Bc) - (A2+B2)*(-N*N*Ac+M*M*Bc+2*N*M*A*B);
}

/*
  This function evaluates the complex coefficients of the algebraic polynomial
  using Fourier transform of the function g(u).

  Input:
    data - structure containing pre-calculated data.
  Output:
    c    - array of complex coefficients of the polynomial.
  Return value: uncertainty of the coefficients.
*/
template<typename realfp>
realfp create_polynomial(const SAuxData<realfp>& data, complex<realfp> c[DEG+1], realfp (*func_g_fast)(int, const SAuxData<realfp>&))
{
 if(!(data.e1<=1)) return -1;

 //direct evaluating the highest coefficients
 const complex<realfp> hc = sqr(data.alpha1*sqr(data.e1)/16)*( data.e2<=1 ?
                      complex<realfp>(data.Pp-data.Ss-data.e1*data.e2, data.Sp+data.Ps)*
                      complex<realfp>(data.Pp-data.Ss+data.e1*data.e2, data.Sp+data.Ps)*
                      complex<realfp>(data.Pp+data.Ss-data.e1*data.e2, data.Sp-data.Ps)*
                      complex<realfp>(data.Pp+data.Ss+data.e1*data.e2, data.Sp-data.Ps) :

                      complex<realfp>(data.Pp-data.Ps-data.e1*data.e2, data.Sp-data.Ss)*
                      complex<realfp>(data.Pp-data.Ps+data.e1*data.e2, data.Sp-data.Ss)*
                      complex<realfp>(data.Pp+data.Ps-data.e1*data.e2, data.Sp+data.Ss)*
                      complex<realfp>(data.Pp+data.Ps+data.e1*data.e2, data.Sp+data.Ss) );
 int i, j;
 realfp vals[DIM];
 for(j=0; j<DIM; j++)
  vals[j] = func_g_fast(j, data);
 for(j=1; j<=DIM/2; j++)
 {
  const realfp sum = vals[j]+vals[DIM-j];
  const realfp dif = vals[j]-vals[DIM-j];
  vals[j]=sum; vals[DIM-j]=dif;
 }

 realfp a=vals[0], b;
 for(j=1; j<=DIM/2; j++) a+=vals[j];
 c[DEG/2]=a/DIM;

 for(i=1; i<DEG/2; i++)
 {
  a=vals[0]; b=0;
  for(j=1; j<=DIM/2; j++)
  {
   a+=vals[j]*real(exps<realfp>()[(j*i)%DIM]);
   b+=vals[DIM-j]*imag(exps<realfp>()[(j*i)%DIM]);
  }
  a/=DIM; b/=DIM;
  c[DEG/2-i]=complex<realfp>(a, b);
  c[DEG/2+i]=complex<realfp>(a,-b);
 }

 c[0]=hc; c[DEG]=conj(hc);
 realfp cerr=0;

 a=vals[0]; b=0;
 for(j=1; j<=DIM/2; j++)
 {
  a+=vals[j]*real(exps<realfp>()[(j*i)%DIM]);
  b+=vals[DIM-j]*imag(exps<realfp>()[(j*i)%DIM]);
 }
 a/=DIM; b/=DIM;

 cerr+=norm(complex<realfp>(a, b)-hc);

 for(i++; i<=DIM/2; i++)
 {
  a=vals[0]; b=0;
  for(j=1; j<=DIM/2; j++)
  {
   a+=vals[j]*real(exps<realfp>()[(j*i)%DIM]);
   b+=vals[DIM-j]*imag(exps<realfp>()[(j*i)%DIM]);
  }
  a/=DIM; b/=DIM;
  cerr+=(a*a+b*b);
 }

 cerr/=((DIM-DEG+1)/2);

 return max(sqrt(cerr),ierr<realfp>()*abs(hc));
}

// this is hyperbolic version for real-valued coefficients
template<typename realfp>
realfp create_polynomial(const SAuxData<realfp>& data, complex<realfp> c[DEG+1], complex<realfp> (*func_g_fast)(int, const SAuxData<realfp>&))
{
 if(!(data.e1>1)) return -1;

 //direct evaluating the highest coefficients
 const realfp d=sqr(data.e2); // from 0 to 1
 const realfp d1p=sqr((data.Pp+data.Sp)/data.e1);
 const realfp d2p=sqr((data.Ps+data.Ss)/data.e1);
 const realfp dpp=d1p+d2p;
 const realfp dpm=d1p-d2p;
 const realfp h  = sqr(data.alpha1*sqr(sqr(data.e1))/16)*( data.e2<=1 ?
                      (dpp*dpp+d*d-2*d*dpm) :
                      (dpm*dpm+d*d-2*d*dpp) );
 const realfp d1m=sqr((data.Pp-data.Sp)/data.e1);
 const realfp d2m=sqr((data.Ps-data.Ss)/data.e1);
 const realfp dmp=d1m+d2m;
 const realfp dmm=d1m-d2m;
 const realfp hc = sqr(data.alpha1*sqr(sqr(data.e1))/16)*( data.e2<=1 ?
                      (dmp*dmp+d*d-2*d*dmm) :
                      (dmm*dmm+d*d-2*d*dmp) );

 int i, j;
 complex<realfp> vals[DIM/2+1];
 const realfp vals0=real(func_g_fast(0, data));  // the only real value is g(0)
// second half of the array is complex-conjugate of the first one taken in reverse order
 for(j=1; j<=DIM/2; j++) vals[j] = func_g_fast(j, data);

 realfp a=0,b;
 for(j=1; j<=DIM/2; j++) a+=real(vals[j]);
 c[DEG/2]=(2*a+vals0)/DIM;

 for(i=1; i<DEG/2; i++)
 {
  a=0; b=0;
  for(j=1; j<=DIM/2; j++)
  {
   a+=real(vals[j])*real(exps<realfp>()[(j*i)%DIM]);
   b+=imag(vals[j])*imag(exps<realfp>()[(j*i)%DIM]);
  }
  c[DEG/2-i]=(2*(a-b)+vals0)/DIM;
  c[DEG/2+i]=(2*(a+b)+vals0)/DIM;
 }

 c[0]=hc; c[DEG]=h;

 a=0; b=0;
 for(j=1; j<=DIM/2; j++)
 {
  a+=real(vals[j])*real(exps<realfp>()[(j*i)%DIM]);
  b+=imag(vals[j])*imag(exps<realfp>()[(j*i)%DIM]);
 }
 realfp cerr=sqr((2*(a-b)+vals0)/DIM-hc)+sqr((2*(a+b)+vals0)/DIM-h);

 for(i++; i<=DIM/2; i++)
 {
  a=0; b=0;
  for(j=1; j<=DIM/2; j++)
  {
   a+=real(vals[j])*real(exps<realfp>()[(j*i)%DIM]);
   b+=imag(vals[j])*imag(exps<realfp>()[(j*i)%DIM]);
  }
  cerr+=sqr((2*a+vals0)/DIM)+sqr(2*b/DIM);
 }

 cerr/=(DIM-DEG+1);

 return max(sqrt(cerr),ierr<realfp>()*max(fabs(hc),fabs(h)));
}

/*
  This function evaluates the distance between two fixed points of the two orbits.

  Input:
    data  - array of the pre-calculated data;
    u1,u2 - eccentric anomalies of the points.
  Return value: distance.
*/
/* template<typename realfp>
realfp distance_between(const SAuxData<realfp>& data, realfp u1, realfp u2)
{
 const realfp x1 = cos(u1);
 const realfp _x1 = x1-data.e1;
 const realfp y1 = sin(u1);
 const realfp x2 = cos(u2);
 const realfp _x2 = x2-data.e2;
 const realfp y2 = sin(u2);
 const realfp R1 = data.a1*(1-data.e1*x1);
 const realfp R2 = data.a2*(1-data.e2*x2);
 return safe_sqrt(R1*R1+R2*R2-2*data.a1*data.a2*(data.Pp*_x1*_x2+data.Sp*y1*_x2+data.Ps*y2*_x1+data.Ss*y1*y2));
}*/

/*
  Computes the radius-vector in an orbit

  Input:
    P,Q - orbital vectors;
    a,e - other Keplerian elements;
    u - eccentric anomaly.
  Output:
    r - radius-vector.
*/
template<typename realfp>
inline void radius_vector(bool br, const realfp P[3], const realfp Q[3], realfp a, realfp e, realfp u, realfp r[3])
{
 realfp p,q;
 if(e<=1) {p=cos(u);q=sin(u)*sqrt(1-e*e);}
 else {const realfp tmp=exp(fabs(u)); p=(tmp+1/tmp)/2;q=(1/tmp-tmp)/2*sign(u)*sqrt(e*e-1);}
 if(!br) {p=-p;q=-q;}
 p-=e;
 r[0] = a*(P[0]*p+Q[0]*q);
 r[1] = a*(P[1]*p+Q[1]*q);
 r[2] = a*(P[2]*p+Q[2]*q);
}

/*
  Computes the radius-vector r and the derivative dr/du

  Input:
    P,Q - orbital vectors;
    a,e - other Keplerian elements;
    u   - eccentric anomaly.
  Output:
    r  - radius-vector.
    rd - dr/du / a.
*/
template<typename realfp>
inline void radius_vector_valder(bool br, const realfp P[3], const realfp Q[3], realfp a, realfp e, realfp u, realfp r[3], realfp rd[3])
{
 realfp x,y,eta,p,q,h;
 if(e<=1) {x=cos(u);y=sin(u);eta=sqrt(1-e*e);h=x*eta;}
 else {const realfp tmp=exp(fabs(u)); x=(tmp+1/tmp)/2;y=(1/tmp-tmp)/2*sign(u);eta=sqrt(e*e-1);h=-x*eta;}
 if(!br) {x=-x;y=-y;h=-h;}
 p=x-e;q=y*eta;
 r[0] = a*(P[0]*p+Q[0]*q);
 r[1] = a*(P[1]*p+Q[1]*q);
 r[2] = a*(P[2]*p+Q[2]*q);
 rd[0] = -P[0]*y+Q[0]*h;
 rd[1] = -P[1]*y+Q[1]*h;
 rd[2] = -P[2]*y+Q[2]*h;
}

/*
  Computes the radius-vector r and the derivatives dr/du, d2r/du2

  Input:
    P,Q - orbital vectors;
    a,e - other Keplerian elements;
    u   - eccentric anomaly.
  Output:
    r   - radius-vector.
    rd  - dr/du / a.
    rdd - d2r/du2 / a.
*/
template<typename realfp>
inline void radius_vector_valder2(bool br, const realfp P[3], const realfp Q[3], realfp a, realfp e, realfp u, realfp r[3], realfp rd[3], realfp rdd[3])
{
 realfp x,y,eta,p,q,h;
 if(e<=1) {x=cos(u);y=sin(u);eta=sqrt(1-e*e);h=x*eta;}
 else {const realfp tmp=exp(fabs(u)); x=(tmp+1/tmp)/2;y=(1/tmp-tmp)/2*sign(u);eta=sqrt(e*e-1);h=-x*eta;}
 if(!br) {x=-x;y=-y;h=-h;}
 p=x-e;q=y*eta;
 r[0] = a*(P[0]*p+Q[0]*q);
 r[1] = a*(P[1]*p+Q[1]*q);
 r[2] = a*(P[2]*p+Q[2]*q);
 rd[0] = -P[0]*y+Q[0]*h;
 rd[1] = -P[1]*y+Q[1]*h;
 rd[2] = -P[2]*y+Q[2]*h;
 rdd[0] = P[0]*x+Q[0]*q;
 rdd[1] = P[1]*x+Q[1]*q;
 rdd[2] = P[2]*x+Q[2]*q;
 if(e<=1) {rdd[0]=-rdd[0];rdd[1]=-rdd[1];rdd[2]=-rdd[2];}
}

/*
  More numerically careful version of distance_between.

  Input:
    data  - array of the pre-calculated data;
    u1,u2 - eccentric anomalies of the points.
  Return value: distance.
*/
template<typename realfp>
realfp distance_between(const SAuxData<realfp>& data, bool br1, bool br2, realfp u1, realfp u2)
{
 realfp r1[3],r2[3];
 radius_vector(br1,data.P1,data.Q1,data.a1,data.e1,u1,r1);
 radius_vector(br2,data.P2,data.Q2,data.a2,data.e2,u2,r2);
 const realfp dr[3] = {r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
 return sqrt(norm(dr));
}

/*
  Computes rho(u1,u2) and its gradient in a numerically careful manner.

  Input:
    data  - array of the pre-calculated data;
    u1,u2 - eccentric anomalies of the points.
  Output:
    g - gradient rho'.
  Return value: rho.
*/
template<typename realfp>
realfp sqdist_vg(const SAuxData<realfp>& data, bool br1, bool br2, realfp u1, realfp u2, realfp g[2])
{
 realfp r1[3],rd1[3],r2[3],rd2[3];
 radius_vector_valder(br1,data.P1,data.Q1,data.a1,data.e1,u1,r1,rd1);
 radius_vector_valder(br2,data.P2,data.Q2,data.a2,data.e2,u2,r2,rd2);
 const realfp dr[3] = {r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
 g[0] = -dot_product(dr,rd1)/data.a2;
 g[1] =  dot_product(dr,rd2)/data.a1;
 return norm(dr)/(2*data.a1*data.a2);
}

/*
  Computes rho(u1,u2), its gradient and Hessian matrix in a numerically careful manner.

  Input:
    data  - array of the pre-calculated data;
    u1,u2 - eccentric anomalies of the points.
  Output:
    g - gradient rho'.
    H - Hessian rho'' (3 independent elements).
  Return value: rho.
*/
template<typename realfp>
realfp sqdist_vgH(const SAuxData<realfp>& data, bool br1, bool br2, realfp u1, realfp u2, realfp g[2], realfp H[3])
{
 realfp r1[3],rd1[3],rdd1[3],r2[3],rd2[3],rdd2[3];
 radius_vector_valder2(br1,data.P1,data.Q1,data.a1,data.e1,u1,r1,rd1,rdd1);
 radius_vector_valder2(br2,data.P2,data.Q2,data.a2,data.e2,u2,r2,rd2,rdd2);
 const realfp dr[3] = {r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
 g[0] = -dot_product(dr,rd1)/data.a2;
 g[1] =  dot_product(dr,rd2)/data.a1;
 H[0] = (-dot_product(dr,rdd1) + norm(rd1)*data.a1)/data.a2;
 H[1] = ( dot_product(dr,rdd2) + norm(rd2)*data.a2)/data.a1;
 H[2] = -dot_product(rd1,rd2);
 return norm(dr)/(2*data.a1*data.a2);
}

/*
  This function evaluates eccentric anomaly on the second orbit. If u1 is a critical point of the distance,
  one of the output values of u2, u2_ (most probably u2) corresponds to the critical point also.

  Input:
    data    - array of the pre-calculated data;
    u1      - eccentric anomaly on the first orbit;

  Output:
    u2, u2_ - two values of the second eccentric anomalies;

  Return value: true if ok, false if there was negative discriminant.
*/
template<typename realfp>
bool eliminated_anomaly(const SAuxData<realfp>& data, realfp u1, realfp& u2, realfp& u2_, realfp* perror=0)
{
 realfp x,y,x_,M,N,V,D;
 if(data.e1<=1) {x=cos(u1); y=sin(u1); x_= x-data.e1;
                 M = data.Sp*y+data.Pp*x_+data.alpha2*data.e2;
                 N =-data.Ss*y-data.Ps*x_;}
 else           {const realfp tmp1=exp(fabs(u1)); x=(tmp1+1/tmp1)/2;y=(tmp1-1/tmp1)/2*sign(u1);x_= x-data.e1;
                 M =-data.Sp*y+data.Pp*x_+data.alpha2*data.e2;
                 N = data.Ss*y-data.Ps*x_;}
 const realfp A = data.Ps*y-data.Ss*x;
 const realfp B = data.Pp*y-data.Sp*x;
 const realfp C = data.e2*B-data.alpha1*data.e1*y*(1-data.e1*x);
 if(data.e2<=1) {V=A*A+B*B;D=V-C*C;}
 else {V=B*B-A*A;D=C*C-V;}

 bool res;
 if(res=(D<0)) {D=0; u2_=u2= data.e2<=1 ? (C>0 ? atan2(A,B) : atan2(-A,-B)) : atan2h(A,B);}
 else
 {
  D = sqrt(D);
  const realfp  AC = A*C;
  const realfp  BD = B*D;
  realfp  y1_= (AC-BD)/V;
  realfp  y2_= (AC+BD)/V;
  if(!(data.e2<=1)) swap(y1_,y2_);
  const realfp  BC = B*C;
  const realfp  AD = A*D;
  realfp  x1_= (BC+AD)/V;
  realfp  x2_= (BC-AD)/V;
  if(data.e2<=1) {u2 = atan2(y1_, x1_);
                  u2_= atan2(y2_, x2_);}
  else           {if(!(x1_>0) && x2_>0) {swap(x1_,x2_);swap(y1_,y2_);}
                  u2 = atan2h(y1_,x1_);
                  u2_= atan2h(y2_,x2_);}

  const realfp err = M*y1_+N*x1_-data.K*y1_*x1_;
  const realfp err_= M*y2_+N*x2_-data.K*y2_*x2_;
  if(!(fabs(err)<fabs(err_)) && (data.e2<=1 || x2_>0)) swap(u2,u2_);
 }
 if(perror!=0) {
    const realfp sigma = *perror*sqrt(1+sqr(data.e2)+sqr(data.alpha1*data.e1));
    *perror = sigma/sqrt(D*D+abs(C)*sigma/2);}

 return res;
}

/*
  This function evaluates candidate minimum distance for a given value
  of the eccentric anomaly on the first orbit.

  Input:
    data     - pre-calculated data;
    u1       - eccentric anomaly on the first orbit;

  Output:
    u2       - candidate value of the second eccentric anomaly.

  Return value: candidate distance (or a negative value if there is no such distance for the requested u1).
*/
template<typename realfp>
realfp min_distance_for(const SAuxData<realfp>& data, realfp u1, realfp& u2, realfp* perror=0)
{
 realfp _u2, _u2_;
 if(eliminated_anomaly(data, u1, _u2, _u2_,perror))
  return -distance_between(data,true,true,u1,_u2);
 const realfp d = distance_between(data,true,true,u1,_u2);
 const realfp d_= distance_between(data,true,true,u1,_u2_);
 if(d<d_)
 {
  u2 = _u2;
  return d;
 }
 else
 {
  u2 = _u2_;
  return d_;
 }
}

//This function calculates the second linking coefficient,
//which is defined by (3) formula in the article (Kholshevnikov, Vassiliev, 1999; CMDA 75, 67-74).
template<typename realfp>
inline realfp l2(const SAuxData<realfp>& data)
{
 const realfp fct1 = data.p2*(data.abs_w+data.e1*data.P1w)-data.p1*(data.abs_w+data.e2*data.P2w);
 const realfp fct2 = data.p2*(data.abs_w-data.e1*data.P1w)-data.p1*(data.abs_w-data.e2*data.P2w);
 return fct1*fct2;
}

//This function calculates the first linking coefficient,
//which is defined by (2) formula in the article (Kholshevnikov, Vassiliev, 1999; CMDA 75, 67-74).
template<typename realfp>
inline pair<realfp,realfp> l1(const SAuxData<realfp>& data)
{
 const realfp costheta1 = data.P1w/data.abs_w;
 const realfp costheta2 = data.P2w/data.abs_w;
 realfp rr = data.p1/(1+data.e1*costheta1)-data.p2/(1+data.e2*costheta2); //calculating r-r'
 realfp RR = data.p1/(1-data.e1*costheta1)-data.p2/(1-data.e2*costheta2); //calculating R-R'
 const realfp mlt = rr*RR;
 const realfp sq = fabs(RR)<fabs(rr) ? sqr(RR) : sqr(rr);
 return make_pair(mlt,sq*sign(mlt));
}

//This function calculates the third linking coefficient,
//which is defined by (5) formula in the article (Kholshevnikov, Vassiliev, 1999; CMDA 75, 67-74).
//This function is supposed to be used when the mutual inclination of a pair of orbits is sufficiently small.
template<typename realfp>
inline realfp l3(const SAuxData<realfp>& data)
{
 return data.a1*data.p1+data.a2*data.p2-2*data.a1*data.a2*(1-data.e1*data.e2*cos(data.w1-data.w2));
}

template<typename realfp>
SLCResult<realfp>::SLCResult(): I(-1.), l(0.0), l2(0.0) {}

template<typename realfp>
SLCResult<realfp> LC(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2, realfp min_mut_incl)
{
 SLCResult<realfp> result;
 SAuxData<realfp> data(O1, O2);
 result.I = data.I;
 result.l2= l2(data);
 if (data.I < min_mut_incl) result.l = l3(data);
 else
 {
  const pair<realfp,realfp> l = l1(data);
  result.l = l.first;
  result.lmod = l.second;
 }
 return result;
}

template<typename realfp>
inline realfp sqdist_2Diter(const SAuxData<realfp>& data, bool br1, bool br2, realfp& u1, realfp& u2, realfp& rho, realfp g[2], realfp H[3], realfp& detH)
{
 rho = sqdist_vgH(data,br1,br2,u1,u2,g,H);
 detH = H[0]*H[1]-H[2]*H[2];
 if(detH==0) return 0;
 const realfp du1 = -(H[1]*g[0]-H[2]*g[1])/detH;
 const realfp du2 = -(H[0]*g[1]-H[2]*g[0])/detH;
 u1 += du1;
 u2 += du2;
 return sqr(du1)+sqr(du2);
}

/*
  This function refines the distance by a single 2D Newtonian iteration and estimates residual uncertainty

  Input:
    data  - array of the pre-calculated data;
    u1,u2 - eccentric anomalies.
    maxcnt - maximum number of iterations

  Output:
    u1, u2  - refined eccentric anomalies;
    u1_err, u2_err  - their numeric uncertainties;
    rho     - value of function rho
    rho_err - its numeric uncertainty
    eps     - desired precision in u1,u2
    iter_cnt- number of iterations
    g       - gradient of rho;
    H       - Hessian of rho.

  Return value: +2(+1) if H is positive-(semi)definite, -2(-1) if negative-(semi)definite, 0 otherwise.
*/
template<typename realfp>
short int Newton_sqdist(const SAuxData<realfp>& data,
                        realfp& u1, realfp& u2,
                        realfp& u1_err, realfp& u2_err,
                        realfp& rho, realfp& rho_err,
                        realfp eps, unsigned long& iter_cnt, int maxcount,
                        realfp g[2], realfp H[3])
{
 static const realfp circ=2*pi<realfp>();
 realfp detH;
 realfp derr2;
 realfp err2=numeric_limits<realfp>::infinity();

 iter_cnt=0;
 do
 {
  const realfp err2_ = sqdist_2Diter(data,true,true,u1,u2,rho,g,H,detH);
  if(data.e1<=1 && fabs(u1)>pi<realfp>()) u1=angle_wrap(u1);
  if(data.e2<=1 && fabs(u2)>pi<realfp>()) u2=angle_wrap(u2);
  if(detH==0) break;
  derr2 = err2-err2_;
  err2 = err2_;
  iter_cnt++;
 }
 while(derr2 > ierr<realfp>()*err2 && err2 > sqr(eps) && iter_cnt <= maxcount);

 rho = sqdist_vg(data,true,true,u1,u2,g);
 u1_err = fabs((H[1]*g[0]-H[2]*g[1])/detH);
 u2_err = fabs((H[0]*g[1]-H[2]*g[0])/detH);
 rho_err = fabs((H[1]*sqr(g[0])+H[0]*sqr(g[1])-2*H[2]*g[0]*g[1])/(2*detH));

 if(detH==0)
 {
  u1_err=0;
  u2_err=0;
  rho_err=0;
  if(H[0]>0 || H[1]>0) return +1;
  if(H[0]<0 || H[1]<0) return -1;
 }

 if(detH>0)
 {
  if(H[0]>0) return +2;
  if(H[0]<0) return -2;
 }
 return 0;
}

template<typename realfp>
SMOIDResult<realfp>::SMOIDResult(): good(true), distance(-1), time(0),
                            distance_error(0), u1_error(0), min_delta(-1),
                            root_count(0), iter_count(0), iter_count_2D(0) {}

template<typename realfp>
SMOIDResult<realfp> MOID_fast(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2,
                              realfp maxrooterr, realfp minrooterr,
                              realfp nu)
{
 const clock_t time = clock();
 SMOIDResult<realfp> result;
 SAuxData<realfp> data(O1, O2);
 if(!(data.e1!=1 && data.e2!=1)) {
   result.good = false;
   result.time = static_cast<realfp>(clock()-time)/CLOCKS_PER_SEC;
   return result;}

 const static bool meticulous = false;  // tests revealed meticulous mode is actually bad

 complex<realfp> c[DEG+1]; //an array of complex Fourier coefficients
 complex<realfp> c_[DEG+1];
// initializing the coefficients and calculating their precision
 realfp cerr = data.e1<=1 ? create_polynomial(data, c, data.e2<=1 ? func_g_ee_fast<realfp> : func_g_eh_fast<realfp>) :
                            create_polynomial(data, c, data.e2<=1 ? func_g_he_fast<realfp> : func_g_hh_fast<realfp>);
 const realfp znrm = data.e1<=1 ? 1 : pow(norm(c[0])/norm(c[DEG]),static_cast<realfp>(1)/(2*DEG));  // geometric mean of all |z_k|
 {
  realfp tmp=1;
  for(int i=1;i<=DEG/2;i++) {tmp*=znrm;c[DEG/2+i]*=tmp;c[DEG/2-i]/=tmp;}
  cerr*=max(tmp,1/tmp);
 }
 copy(c,c+DEG+1,c_);

 complex<realfp> roots[DEG]; //an array of the complex roots
 fill(roots,roots+DEG,0);

 {
// first four root are seeked automatically
  int idx= 4;
  const realfp maxabsroot= max_root_bound(DEG,c,true)+ierr<realfp>();
  const realfp minabsroot= data.e1<=1 ? 1/maxabsroot : 1/max_root_bound(DEG,c,false)+ierr<realfp>();
  if(idx>0)
  {
   if(data.e1<=1)
    roots[0] = polar<realfp>(minabsroot,random<realfp>()*2*pi<realfp>());
   else
   {
    roots[0] = polar<realfp>(1/znrm, 2*pi<realfp>()/5);
    roots[2] = polar<realfp>(1/znrm, data.e2<=1 ? 3*pi<realfp>()/5 : pi<realfp>()/2);
   }
  }

// first roots are seeked near the orbital nodes and +-90 degrees from them
  const realfp eta1=sqrt(fabs(1-sqr(data.e1)));
  const realfp eta2=sqrt(fabs(1-sqr(data.e2)));
  const realfp wPQ1=hypot(data.P1w,data.Q1w);
  const realfp wPQ2=hypot(data.P2w,data.Q2w);
  realfp u1__, u2__, rho, g[2], H[3], detH;
  const realfp u1min=data.e1<=0?0:-log(maxabsroot*znrm);
  const realfp u1max=data.e1<=0?0:-log(minabsroot*znrm);

 // ascending node
  const realfp u1 =atan2_smart( data.Q1w*eta1, data.e1*wPQ1+data.P1w, data.e1<=1); // always exists
  const bool br1 = wPQ1+data.e1*data.P1w>=0;
  const realfp u2 =atan2_smart( data.Q2w*eta2, data.e2*wPQ2+data.P2w, data.e2<=1); // always exists
  const bool br2 = wPQ2+data.e2*data.P2w>=0;

 // descending node
  const realfp u1_ =atan2_smart(-data.Q1w*eta1, data.e1*wPQ1-data.P1w, data.e1<=1); // always exists
  const bool br1_ = wPQ1-data.e1*data.P1w>=0;
  const realfp u2_ =atan2_smart(-data.Q2w*eta2, data.e2*wPQ2-data.P2w, data.e2<=1); // always exists
  const bool br2_ = wPQ2-data.e2*data.P2w>=0;

  {
   u1__=u1;u2__=u2;
   sqdist_2Diter(data,br1,br2,u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1 ?1:-1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }

  if((data.e1<=1 || br1_!=br1) && (data.e2<=1 || br2_!=br2))
  {
   u1__=u1_;u2__=u2_;
   sqdist_2Diter(data,br1_,br2_,u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1_?1:-1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }
  else
  {
   u1__= data.e1<=1 || br1_!=br1 ? u1_ : (u1+u1_)/2;
   u2__= data.e2<=1 || br2_!=br2 ? u2_ : (u2+u2_)/2;
   sqdist_2Diter(data, data.e1<=1 || !br1, data.e2<=1 || !br2, u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1 ?-1:1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }

// in the EE-case set up 2 additional roots
  if(data.e1<=1 && data.e2<=1)
  {
 // +pi/2 from asc.node (in true anomalies)
   u1__ =atan2(-data.P1w*eta1, data.e1*wPQ1+data.Q1w);
   u2__ =atan2(-data.P2w*eta2, data.e2*wPQ2+data.Q2w);
   realfp u2_=atan2( data.P2w*eta2, data.e2*wPQ2-data.Q2w);
   if(2*data.I>pi<realfp>()) swap(u2__,u2_);
   sqdist_2Diter(data,true,true,u1__,u2__,rho,g,H,detH);
   roots[idx++] = polar<realfp>(1, u1__);

 // -pi/2 from asc.node (in true anomalies)
   u1__ =atan2( data.P1w*eta1, data.e1*wPQ1-data.Q1w);
   sqdist_2Diter(data,true,true,u1__,u2_,rho,g,H,detH);
   roots[idx++] = polar<realfp>(1,u1__);
  }
  else
  {

  if(data.e1<=1 || br1_!=br1)
  {
   u1__=u1_;u2__=u2;
   sqdist_2Diter(data,br1_,br2,u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1_?1:-1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }
  else
  {
   u1__=(u1+u1_)/2;u2__=u2;
   sqdist_2Diter(data,!br1,br2,u1__,u2__,rho,g,H,detH);
   u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = exp(-u1__)*(br1 ?-1:1)/znrm;
  }

  if(data.e2<=1 || br2_!=br2)
  {
   u1__=u1;u2__=u2_;
   sqdist_2Diter(data,br1,br2_,u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1 ?1:-1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }
  else
  {
   u1__=u1;u2__=(u2+u2_)/2;
   sqdist_2Diter(data,br1,!br2,u1__,u2__,rho,g,H,detH);
   if(!(data.e1<=1)) u1__=min(max(u1__,u1min),u1max);
   roots[idx++] = polar<realfp>(data.e1<=1 ? 1 : exp(-u1__)*(br1 ?1:-1)/znrm, (data.e1<=1 ? u1__ : random<realfp>()*1e-2+1e-3));
  }

  }

  if(data.e1<=1)
   roots[idx] = polar<realfp>(minabsroot,random<realfp>()*2*pi<realfp>());  // remaining roots
  else
  {
   roots[idx]   = polar<realfp>(1/znrm, 2*pi<realfp>()/5);
   roots[idx+2] = polar<realfp>(1/znrm, data.e2<=1 ? 3*pi<realfp>()/5 : pi<realfp>()/2);
  }
 }

 result.iter_count = polynomial_roots(DEG, c_, roots, minrooterr, maxrooterr, data.e1<=1);

 int root_index = 0;
 realfp u1, u2;

 for(int i=0; i<DEG; i++)
 {
  const realfp rooterr = root_error(DEG, c, cerr, roots[i])*nu; // relative uncertainty for the root
  const realfp dlt = data.e1<=1 ? fabs(log(norm(roots[i])))/(2*rooterr) :    //deviation of the root from the unit circle
                                                                             //(relatively to the expected error)
                                  (roots[i]!=complex<realfp>(0) ? fabs(arg(roots[i])) : 0)/rooterr;    //deviation of the root from real positive ray
  if(meticulous || dlt<10)
  {
   result.good &= (rooterr<maxrooterr);
   if(meticulous || dlt<3)
   {
    result.root_count++;
    u1 = data.e1<=1 ? arg(roots[i]) : (real(roots[i])>0 ? -log(real(roots[i])*znrm) : numeric_limits<realfp>::infinity());

    const realfp dst = fabs(min_distance_for(data,u1,u2));

    if(dst>=0 && ((result.distance)<0 || dst<(result.distance)))
    {
     result.distance = dst;
     result.u1 = u1;
     result.u2 = u2;
     root_index = i;
    }
   }
  }
  if(!(meticulous || dlt<3) && ((result.min_delta)<0 || dlt<(result.min_delta)))
    result.min_delta = dlt;
 }

 realfp g[2];
 realfp H[3];
 realfp du1=result.u1, du2=result.u2;
 const short int Hsign = Newton_sqdist(data,
                                       result.u1,result.u2,result.u1_error,result.u2_error,
                                       result.distance,result.distance_error,
                                       minrooterr,result.iter_count_2D,30,g,H)*(data.a1*data.a2>=0?1:-1);
 du1-=result.u1; if(data.e1<=1 && fabs(du1)>pi<realfp>()) du1=angle_wrap(du1);
 du2-=result.u2; if(data.e2<=1 && fabs(du2)>pi<realfp>()) du2=angle_wrap(du2);
 const realfp du = hypot(du1,du2);

 const realfp angerr=ierr<realfp>()*pi<realfp>()*nu;
 result.u1_error += angerr;
 result.u2_error += angerr;

 const realfp tmp = 2*data.a1*data.a2;
 const realfp x1 = data.e1<=1 ? cos(result.u1) : cosh(result.u1);
 const realfp x2 = data.e2<=1 ? cos(result.u2) : cosh(result.u2);
 const realfp tmp1 = data.e1*x1;
 const realfp tmp2 = data.e2*x2;
 const realfp r1 = data.a1*(1-tmp1);
 const realfp r2 = data.a2*(1-tmp2);
 const realfp rd1 = fabs(data.a1)*sqrt(fabs(1-sqr(tmp1)));
 const realfp rd2 = fabs(data.a2)*sqrt(fabs(1-sqr(tmp2)));
 const realfp raderr = hypot(r1,r2)*ierr<realfp>()*nu;
 const realfp graderr = raderr*hypot(rd1,rd2)/fabs(tmp);

 const realfp eigenval_max = fabs(H[0]+H[1])/2+hypot((H[0]-H[1])/2,H[2]);
 const realfp detH = fabs(H[0]*H[1]-H[2]*H[2]);
 const realfp uerr3 = eigenval_max/detH*graderr;
 result.u1_error += uerr3;
 result.u2_error += uerr3;
 result.distance_error += eigenval_max/2*(sqr(angerr)+sqr(graderr)/detH);

 result.distance *= tmp;
 result.distance_error *= fabs(tmp);

 result.distance = sqrt(result.distance);
 result.distance_error += 2*result.distance*raderr+sqr(raderr);
 result.distance_error /= sqrt(sqr(result.distance)+result.distance_error/2);

 //some diagnostic of the quality:
 result.good = result.good &&
               result.root_count>=(data.e1<=1?4:2) && result.root_count%2==0 &&
               (!meticulous && result.min_delta>10) &&
               Hsign==+2 &&
               fabs(result.u1_error)<=maxrooterr;

 result.time = static_cast<realfp>(clock()-time)/CLOCKS_PER_SEC;

 return result;
}

//This function finds the minimum value of the distance in the segment u1\in (a,b) on a grid of N points.
template<typename realfp>
int search_in_segment(const SAuxData<realfp>& data, realfp a, realfp b, unsigned int N, realfp& min_val)
{
 min_val = -1;
 int min_i;
 for(int i=0; i<=N; i++) {
   const realfp u1 = (a*(N-i)+b*i)/N;
   realfp u2;
   const realfp d = min_distance_for(data,u1,u2);
   if(d<0) continue;
   if(d<min_val || min_val<0) {
     min_val=d;
     min_i=i;}}

 return min_i;
}

template<typename realfp>
int restrict_search_range(const SAuxData<realfp>& data, pair<realfp,realfp>& u1, pair<realfp,realfp>& u1_)
{
 const realfp theta1 = atan2(data.Q1w, data.P1w);
 const realfp theta2 = atan2(data.Q2w, data.P2w);
 const realfp costheta1 = cos(theta1);
 const realfp costheta2 = cos(theta2);
 const realfp rp = data.p1/(1+data.e1*costheta1);
 const realfp rm = data.p1/(1-data.e1*costheta1);
 const realfp Rp = data.p2/(1+data.e2*costheta2);
 const realfp Rm = data.p2/(1-data.e2*costheta2);
 const realfp d[4] = {fabs(rp-Rp), fabs(rp+Rm), fabs(rm+Rp), fabs(rm-Rm)};
 const bool flag[4]= {rp>0 && Rp>0, rp>0 && Rm>0, rm>0 && Rp>0, rm>0 && Rm>0};
 realfp dO=-1;
 for(int i=0;i<4;i++)
  if(flag[i] && (dO>d[i] || dO<0))
   dO=d[i];
 const realfp A = sqrt(fabs(1-sqr(data.e1*costheta1)));
 const realfp k = dO/(A*data.abs_w*fabs(data.a1));
 const realfp sintheta1 = sin(theta1);
 const realfp eta1 = sqrt(fabs(1-sqr(data.e1)));
 realfp phi = atan2_smart(sintheta1, costheta1*eta1, data.e1<=1);
 realfp esth= data.e1*sintheta1/A;  // e sin(phi) for the E-case
 if(data.e1<=1)
  if(min<realfp>(fabs(esth),1)<=fabs(1-k))
   if(k<1)
   {
    const realfp tmpm = safe_asin(esth-k);
    const realfp tmpp = safe_asin(esth+k);
    u1 = make_pair(phi-tmpp,phi-tmpm);
    phi += pi<realfp>();
    u1_= make_pair(phi+tmpm,phi+tmpp);
    return 2;
   }
   else u1=make_pair(-pi<realfp>(),pi<realfp>());
  else if(esth>0)
       {
        const realfp tmpm = safe_asin(esth-k);
        u1 = make_pair(phi+tmpm,phi+pi<realfp>()-tmpm);
        if(u1.first>u1.second) u1.first=u1.second=phi+pi<realfp>()/2;
       }
       else
       {
        const realfp tmpp = safe_asin(esth+k);
        u1 = make_pair(phi-tmpp,phi+pi<realfp>()+tmpp);
        if(u1.first>u1.second) u1.first=u1.second=phi-pi<realfp>()/2;
       }
 else if(fabs(costheta1)*eta1<=fabs(sintheta1))
       {
        if(sintheta1<0) esth=-esth;
        const realfp tmpm = safe_acosh(esth-k);
        const realfp tmpp = safe_acosh(esth+k);
        if(tmpm>0)
        {
         u1 = make_pair(-tmpp-phi,-tmpm-phi);
         u1_= make_pair( tmpm-phi, tmpp-phi);
         return 2;
        }
        else u1 = make_pair(-tmpp-phi,tmpp-phi);
       }
      else
      {
       if(costheta1<0) esth=-esth;
       u1 = make_pair(asinh(esth-k)-phi,asinh(esth+k)-phi);
      }
 return 1;
}

template<typename realfp>
int restrict_search_range(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2, pair<realfp,realfp>& u1, pair<realfp,realfp>& u1_)
{
 return restrict_search_range(SAuxData<realfp>(O1,O2),u1,u1_);
}

template<typename realfp>
SMOIDResult<realfp> MOID_direct_search(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2,
                                       const unsigned int densities[],
                                       realfp max_dist_error, realfp max_anom_error)
{
 const clock_t time = clock();
 SMOIDResult<realfp> result[2];

 SAuxData<realfp> data(O1, O2);
 pair<realfp,realfp> u1[2];
 int num;
 bool swapped=false;

// select the order (O1,O2) or (O2,O1) that implies the smallest length of the scan range
 {
  pair<realfp,realfp> u1_[2];
  SAuxData<realfp> data_(O2, O1);
  num = restrict_search_range(data,u1[0],u1[1]);
  const int num_ = restrict_search_range(data_,u1_[0],u1_[1]);
  realfp length=0,length_=0;
  for(int idx=0; idx<num; idx++)
   length += max<realfp>(u1[idx].second-u1[idx].first,0);
  for(int idx=0; idx<num_; idx++)
   length_ += max<realfp>(u1_[idx].second-u1_[idx].first,0);
  if(length_<length)
  {
   data=data_;
   u1[0]=u1_[0];
   u1[1]=u1_[1];
   num=num_;
   swapped=true;
  }
 }

 for(int idx=0; idx<num; idx++)
 {
  realfp h = (u1[idx].second-u1[idx].first)/2;
  realfp a = u1[idx].first;
  realfp b = u1[idx].second;

  if(h>0)
  {
   realfp min_val;
   realfp delta;
   bool is_dens = (densities != NULL);
   unsigned int _dens = 3;
   for(int i=0; i<1 || h>max_anom_error || fabs(delta)>max_dist_error; i++)
   {
     if(is_dens) if(densities[i]!=0) if(i>0) _dens = densities[i];
                                     else _dens = ceil(h*densities[i]/pi<realfp>());
                 else is_dens = false;
     if(_dens<3) _dens = 3;
     h = 2*h/_dens;
     realfp min_val_;
     a += h*search_in_segment(data, a, b, _dens, min_val_);
     b = a+h; a -= h;
     if(i>0) delta = min_val_-min_val;
     min_val = min_val_;
   }
   result[idx].u1_error = h;
  }
  else result[idx].u1_error = 0;

  result[idx].u1 = (a+b)/2;
  result[idx].distance = fabs(min_distance_for(data, result[idx].u1, result[idx].u2));
 }

 SMOIDResult<realfp>& res = result[num<2?0:((result[1].distance)<(result[0].distance)?1:0)];

 const realfp eps= sqrt(ierr<realfp>());
 realfp u2;
 const realfp d = fabs(min_distance_for(data, res.u1+eps, u2));
 res.distance_error = fabs(d-res.distance)*sqr(res.u1_error/eps);  // sqr assumes that first derivative vanished due to minimum
 res.u2_error = fabs(u2-res.u2)*res.u1_error/eps;                  // here first derivative is nonzero

 if(swapped)
 {
  swap(res.u1,res.u2);
  swap(res.u1_error,res.u2_error);
 }

  //diagnostic of the quality:
 res.good = res.good && res.distance_error<max_dist_error;

 res.time = static_cast<realfp>(clock()-time)/CLOCKS_PER_SEC;
 return res;
}

#define instantiate_all(realfp) \
template class COrbitData<realfp>; \
template COrbitData<realfp>::COrbitData(const COrbitData<float>&); \
template COrbitData<realfp>::COrbitData(const COrbitData<double>&); \
template COrbitData<realfp>::COrbitData(const COrbitData<long double>&); \
template class SMOIDResult<realfp>; \
template class SLCResult<realfp>; \
template void detect_suitable_options(realfp&, realfp&, realfp&); \
template bool test_peri_apo(const COrbitData<realfp>&, const COrbitData<realfp>&, realfp); \
template SLCResult<realfp> LC(const COrbitData<realfp>&, const COrbitData<realfp>&, realfp); \
template SMOIDResult<realfp> MOID_fast(const COrbitData<realfp>&, const COrbitData<realfp>&, realfp, realfp, realfp); \
template SMOIDResult<realfp> MOID_direct_search(const COrbitData<realfp>&, const COrbitData<realfp>&, const unsigned int[], realfp, realfp); \
template int restrict_search_range(const COrbitData<realfp>& O1, const COrbitData<realfp>& O2, pair<realfp,realfp>& u1, pair<realfp,realfp>& u1_);

template class CInitializer<DIM>;

instantiate_all(float)
instantiate_all(double)
instantiate_all(long double)
