#include "MatrixBLAS.h"
#include "spinblock.h"
#include "initblocks.h"
#include "input.h"
#include "timer.h"
#include <ctime>
#include "rotationmat.h"
#include "wavefunction.h"
#include "global.h"
#include "sweep.h"
#include <unordered_map>
#include <unordered_set>
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <random>
#include "IntegralMatrix.h"
#include <chrono>
#include "stochasticpt_new.h"
#include "sampling.h"
#include "heatbath.h"
#include "nonspinmps.h"
#include <algorithm>
#include <array>
#include <set>
#include <stack>
#include <tuple>
#include <vector>
#include "couplingCoeffs.h"
#include "CIcoeff_from_mps_TRIE_new.h"
#include <math.h>
#include <numeric> 
using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;
int bitstring::n_orb;
void ReadInput(char* conf);

struct MPS2CI;

struct MPS2CI{
    std::vector<double> vals;
    std::vector<std::pair< std::vector<int>, std::vector<int> > > idx_ex;
    int n_sites;
    bool enable_look_up;

    MPS2CI(int n_sites, bool enable_look_up = false)
        : n_sites(n_sites), enable_look_up(enable_look_up) {
    }

//lsh need to change
    // clear trie
    void clear() {
        vals.clear(), idx_ex.first.clear(), idx_ex.second.clear(), idx_ex.clear();
    }
    // deep copy
    std::shared_ptr<MPS2CI> copy() {
        std::shared_ptr<MPS2CI> dett =
            std::make_shared<MPS2CI>(n_sites, enable_look_up);
        dett->type = std::vector<int>(type.begin(), type.end());
        dett->ex   = std::vector<std::vector<int>>(ex.begin(), ex.end());
        dett->vals = std::vector<double>(vals.begin(), vals.end());
        return dett;
    }
    // number of determinants
    size_t size() const noexcept { return vals.size(); }
//    // add a determinant to trie
//    // det[i] = 0 (empty) 1 (beta) 2 (alpha) 3 (alpha beta)
//    void push_back(const std::vector<uint8_t> &det) {
//        assert((int)det.size() == n_sites);
//        size_t cur = 0;
//        for (int i = 0; i < n_sites; i++) {
//            uint8_t j = det[i];
//            if (data[cur][j] == 0) {
//                data[cur][j] = (size_t)data.size();
//                data.push_back(std::array<size_t, 4>{0, 0, 0, 0});
//            }
////            //lsh test
////            std::cout << "i, j, cur = " << i << " "
////             << unsigned(j) << " " << cur << " data[cur][j] = " << data[cur][j] << std::endl;
//
//            cur = data[cur][j];
//
//        }
//        // cannot push_back repeated determinants
//        assert(dets.size() == 0 || cur > dets.back());
//        dets.push_back(cur);
//        if (enable_look_up) {
//            invs.resize(data.size());
//            cur = 0;
//            for (int i = 0; i < n_sites; i++) {
//                uint8_t j = det[i];
//                invs[data[cur][j]] = cur;
//                cur = data[cur][j];
//            }
//        }
//
//
//    }
//    // find the index of a determinant
//    size_t find(const std::vector<uint8_t> &det) {
//        assert((int)det.size() == n_sites);
//        size_t cur = 0;
//        for (int i = 0; i < n_sites; i++) {
//            uint8_t j = det[i];
//            if (data[cur][j] == 0)
//                return -1;
//            cur = data[cur][j];
//        }
//        return (size_t)(lower_bound(dets.begin(), dets.end(), cur) - dets.begin());
//    }
//    // get a determinant in trie
//    std::vector<uint8_t> operator[](size_t idx) const {
//        assert(enable_look_up && idx < dets.size());
//        std::vector<uint8_t> r(n_sites, 0);
//        for (size_t cur = dets[idx], i = n_sites - 1, ir; i >= 0; i--, cur = ir) {
//            ir = invs[cur];
//            for (uint8_t j = 0; j < (uint8_t)data[ir].size(); j++)
//                if (data[ir][j] == cur) {
//                    r[i] = j;
//                    break;
//                }
//        }
//        return r;
//    }
    // set the value for each determinant to the overlap between mps
    void evaluate(const std::shared_ptr<simplemps> &mps, int nocc, int nvir, double cutoff = 1e-10, int max_ex = 4) {
        std::vector<partial_ex> l_idxs, r_idxs;
        std::vector<partial_wv> lwaves, rwaves;
        // initialize lwaves
        partial_wv lw;
        RowVector oldwave(1);
        oldwave(1) = 1.0;
        lw[std::make_tuple(0, 0, 0, 0)] = oldwave;
        int Ms = 0;//This is 2*Sz of left quanta num.
        lwaves.push_back(lw);

        for (int i_site = 0; i_site < n_sites; i_site++)
        {
            int rw_len = 0;
            for (auto& [lwave, l_idx] : zip(lwaves, l_idxs))
            { 
                int na_ex = l_idx.first.size();
                int nb_ex = l_idx.second.size();
                // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                for (int dot = 0; dot < 4; dot++)
                {
                    // update indices of configurations
                    if (i_site < nocc)   // idx for hole
                    {
                        if (na_ex+nb_ex == max_ex   && dot < 3)  continue;
                        if (na_ex+nb_ex == max_ex-1 && dot == 0) continue;
                        rw_len += 1;
                        partial_ex r_idx(l_idx);  // copy constructor
                        if (dot == 0){
                            r_idx.first.resize(rw_len);
                            r_idx.first.push_back(i_site);
                            r_idx.second.resize(rw_len);
                            r_idx.second.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                        else if (dot == 1){
                            r_idx.first.resize(rw_len);
                            r_idx.first.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                        else if (dot == 2){
                            r_idx.second.resize(rw_len);
                            r_idx.second.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                    } 
                    else                // idx for particle 
                    {
                        if (na_ex+nb_ex == max_ex   && dot > 0)  continue;
                        if (na_ex+nb_ex == max_ex-1 && dot == 3) continue;
                        rw_len += 1;
                        partial_ex r_idx(l_idx);  // copy constructor
                        if (dot == 1){
                            r_idx.second.resize(rw_len);
                            r_idx.second.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                        else if (dot == 2){
                            r_idx.first.resize(rw_len);
                            r_idx.first.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                        else if (dot == 3){
                            r_idx.first.resize(rw_len);
                            r_idx.first.push_back(i_site);
                            r_idx.second.resize(rw_len);
                            r_idx.second.push_back(i_site);
                            r_idxs.push_back(r_idx);
                        } 
                    } 

            //   RowVector(lwave) SiteMatrices                     RowVector(rwave)
            //                        dot
            //                         |                        =  
            // leftq --*-- leftq .   --*-- rightq = LefttoRight      --*--
            //      i_site          i_site+1                        i_site+1
            //input:  i_site, dot, leftq, Ms, dotMs,
            //output: rwave<SQ,RowVector>
                    partial_wv rwave;
                    mps->lw_dot_mat_new(i_site, dot, lwave, rwave);

                    if (cutoff != 0) 
                    {
                        double sqsum = 0;
                        for (auto &m : rwave) 
                        { 
                            sqsum += dotproduct(m.second.first, m.second.first); 
                        }
                        if (sqrt(sqsum) < cutoff)
                            continue;
                    }
                    rwaves.resize(rw_len);
                    rwaves.push_back(rwave);
                }
            }

            if (i_sites < n_sites-1){ 
                for (auto& [lwave, l_idx] : zip(lwaves, l_idxs))
                { 
                    lwave.clear();
                    l_idx.clear();
                } 
                lwaves = rwaves;
                l_idxs = r_idxs;

                for (auto& [rwave, r_idx] : zip(rwaves, r_idxs))
                { 
                    rwave.clear();
                    r_idx.clear();
                } 
            }
            else{ 
                vals.resize(rw_len);

                int tmp = 0;
                for (auto &m : rwave) { 
                    assert(m.second.Ncols() == 1);
                    vals.push_back(m.second.element(0));
                    tmp += 1; 
                }
                assert(tmp == 1);
                idx_ex = r_idxs
            }
        }

    }
};

int alpha_count(std::vector<uint8_t> det, const uint8_t qi)
{
    int ncount = 0;
    for (int i=0; i<qi; i++)
        if (det[i] == 3 || det[i] == 1)
            ncount += 1;
    return ncount;
}

// det[i] = 0 (empty) 1 (beta) 2 (alpha) 3 (alpha beta)
double parity_ab_str(std::vector<uint8_t> det)
{
    // | 1_alpha 1_beta ... noc_alpha noc_beta > 
    // = (-1)**n * | 1_beta ...  noc_beta > | 1_alpha ...  noc_alpha > 
    // = (-1)**n * | noc_beta ...  1_beta > | noc_alpha ...  1_alpha > 
    int n=0;
    for (int i=0; i<det.size(); i++)
        if (det[i]==3 || det[i]==2)
           n += alpha_count(det, i);
    return pow(double(-1), n);
}

double parity_ci_to_cc(int sum_ijkl, int n_excite, int nocc)
{
    // For example of singly-excited configuration, parity is
    // | a noc noc-1 ... i+1 i-1 ... 2 1 > = parity_ci_to_cc * | 1 2 ... i-1 a i+1 ... noc-1 noc >
    return pow(double(-1), n_excite * nocc - sum_ijkl - n_excite*(n_excite+1)/2);
}

double parity_reorder(std::vector<uint8_t> det, std::vector<int> reorder, int nocc)
{
    // parity from reorder | 1_beta ...  noc_beta > and | 1_alpha ...  noc_alpha >

    int num_swaps_alpha = 0, num_swaps_beta = 0;
    int tmp;
    std::vector<int> alpha, beta;

    // |0>, |beta>, |alpha>, |alpha beta>  (for j = 0, 1, 2, 3)
    for (int i=0; i<det.size(); i++) {
        if (det[i]==3 || det[i]==2) alpha.push_back(reorder[i]);
        if (det[i]==3 || det[i]==1) beta.push_back(reorder[i]);
    }     

//    std::cout << "alpha" << std::endl;
//    for (auto idet : alpha)
//       std::cout << " " << unsigned(idet);
//    std::cout << std::endl;
//    std::cout << "beta " << std::endl;
//    for (auto idet : beta)
//       std::cout << " " << unsigned(idet);
//    std::cout << std::endl;

    for (int i=0;   i<nocc-1; i++){
    for (int j=i+1; j<nocc; j++){
        if (alpha[i] > alpha[j]){
           tmp = alpha[i];
           alpha[i] = alpha[j];
           alpha[j] = tmp; 
           num_swaps_alpha += 1;
        }
        if (beta[i] > beta[j]){
           tmp = beta[i];
           beta[i] = beta[j];
           beta[j] = tmp; 
           num_swaps_beta += 1;
        }
    }
    }

//    std::cout << "alpha swap: " << num_swaps_alpha << std::endl;
//    std::cout << "beta  swap: " << num_swaps_beta  << std::endl;

    return pow(double(-1), num_swaps_alpha + num_swaps_beta);
}


//int factorial(int n) 
//{ 
//    return (n==1 || n==0) ? 1: n * factorial(n - 1);  
//} 

void MPS2CI_run(double cutoff, int max_ex)
{
  std::shared_ptr<simplemps> zeromps;
  if(dmrginp.spinAdapted())
  {
    std::shared_ptr<NonAbelianmps> mps = std::make_shared<NonAbelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    mps->build(0);
    zeromps = mps;
  }
  else{
    assert(False);
  }
  int nsite = dmrginp.last_site(); // n orb
  int nelec = dmrginp.total_particle_number();
  int nocc  = dmrginp.total_particle_number()/2;
  int nvir  = nsite - nocc; 

//  int idxS = nocc * nvir;
//  int idxD = nocc * (nocc-1) * nvir * (nvir-1) / 4; // pow(factorial(2), 2)
//  size_t idxT = (size_t) nocc * (nocc-1) * (nocc-2) * nvir * (nvir-1) * (nvir-2) / 36; // pow(factorial(3), 2)
 // int idxQ = nocc * (nocc-1) * (nocc-2) * (nocc-3) * nvir * (nvir-1) * (nvir-2) * (nvir-3) / 576; // pow(factorial(4), 2)

//  std::vector<double> Ref(1, 0.0), S_a(idxS, 0.0), D_aa(idxD, 0.0), T_aaa(idxT, 0.0);
//  std::vector<std::vector<double>> D_ab(idxS, std::vector<double>(idxS, 0.0));
//  std::vector<std::vector<double>> T_aab(idxD, std::vector<double>(idxS, 0.0));
//  std::vector<std::vector<double>> Q_aaab(idxT, std::vector<double>(idxS, 0.0)); 
//  std::vector<std::vector<double>> Q_aabb(idxD, std::vector<double>(idxD, 0.0));

  std::vector<double>  parity, reorder_parity;

////============================== test ========================================  
//  //lsh test: ref
//  MPS2CI dtrie(nsite, true); 
//  std::vector<int> reorder = dmrginp.reorder_vector(), rev_reorder(nocc+nvir, 0);
//  // reorder [ lattice_index ] = mo_index, rev_reorder [ mo_index ] = lattice_index
//  // Ref append
//  std::vector<uint8_t> det(nocc+nvir, 0), tmpdet, tmpdet2;
//  for (int i = 0; i < nocc+nvir; i++){
//      std::cout << i << " " << reorder[i] << std::endl; 
//      rev_reorder[reorder[i]] = i;
//      if (reorder[i]<nocc)
//          det[i] = 3;  
//  }
//
////  std::cout << "nocc = " << nocc << " nvir = " << nvir << std::endl;
////  for (auto idet : det)
////     std::cout << " " << unsigned(idet);
////  std::cout << std::endl;
////  dtrie.push_back(det);  
////  reorder_parity.push_back(1.0);
////
////  //lsh test: 1st D
////  tmpdet.assign(det.begin(), det.end());
////  tmpdet[rev_reorder[nocc-1]] = 0;
////  tmpdet[rev_reorder[nocc]] = 3;
////  for (auto idet : tmpdet)
////     std::cout << " " << unsigned(idet);
////  std::cout << std::endl;
////
////  dtrie.push_back(tmpdet);  
////  double pD1 = parity_ab_str(tmpdet);
////  double pD2 = parity_reorder(tmpdet, reorder, nocc);
////  double pD3 = parity_ci_to_cc(nocc-1, 1, nocc);
////  double pD4 = parity_ci_to_cc(nocc-1, 1, nocc);
//
//  //lsh test: 1st S 
//  tmpdet.assign(det.begin(), det.end());
//
//  int i = 2;
//  int a = 0; 
//
//  std::cout << i << " " << a << ": " << rev_reorder[i] << " " << rev_reorder[a+nocc] << std::endl;
//  tmpdet[rev_reorder[i]] = 1;
//  tmpdet[rev_reorder[a+nocc]] = 2;
//  for (auto idet : tmpdet)
//     std::cout << " " << unsigned(idet);
//  std::cout << std::endl;
//
//  dtrie.push_back(tmpdet);  
//  double pS1 = parity_ab_str(tmpdet);
//  double pS2 = parity_reorder(tmpdet, reorder, nocc);
//  double pS3 = parity_ci_to_cc(i, 1, nocc);
//
//
//  tmpdet.assign(det.begin(), det.end());
//  i = 1;
//  a = 0; 
//
//  std::cout << i << " " << a << ": " << rev_reorder[i] << " " << rev_reorder[a+nocc] << std::endl;
//  tmpdet[rev_reorder[i]] = 1;
//  tmpdet[rev_reorder[a+nocc]] = 2;
//  for (auto idet : tmpdet)
//     std::cout << " " << unsigned(idet);
//  std::cout << std::endl;
//
//  dtrie.push_back(tmpdet);  
//  double pD1 = parity_ab_str(tmpdet);
//  double pD2 = parity_reorder(tmpdet, reorder, nocc);
//  double pD3 = parity_ci_to_cc(i, 1, nocc);
//
//
//  tmpdet.assign(det.begin(), det.end());
//  i = 3;
//  a = 5; 
//
//  std::cout << i << " " << a << ": " << rev_reorder[i] << " " << rev_reorder[a+nocc] << std::endl;
//  tmpdet[rev_reorder[i]] = 1;
//  tmpdet[rev_reorder[a+nocc]] = 2;
//  for (auto idet : tmpdet)
//     std::cout << " " << unsigned(idet);
//  std::cout << std::endl;
//
//  dtrie.push_back(tmpdet);  
//  double pT1 = parity_ab_str(tmpdet);
//  double pT2 = parity_reorder(tmpdet, reorder, nocc);
//  double pT3 = parity_ci_to_cc(i, 1, nocc);
//
//
//  tmpdet.assign(det.begin(), det.end());
//  i = 0;
//  a = 5; 
//
//  std::cout << i << " " << a << ": " << rev_reorder[i] << " " << rev_reorder[a+nocc] << std::endl;
//  tmpdet[rev_reorder[i]] = 1;
//  tmpdet[rev_reorder[a+nocc]] = 2;
//  for (auto idet : tmpdet)
//     std::cout << " " << unsigned(idet);
//  std::cout << std::endl;
//
//  dtrie.push_back(tmpdet);  
//  double pQ1 = parity_ab_str(tmpdet);
//  double pQ2 = parity_reorder(tmpdet, reorder, nocc);
//  double pQ3 = parity_ci_to_cc(i, 1, nocc);
//
//
//
////  //lsh test: major Q 
////  tmpdet.assign(det.begin(), det.end());
////  tmpdet[rev_reorder[nocc-2]] = 0;
////  tmpdet[rev_reorder[nocc-1]] = 0;
////  tmpdet[rev_reorder[nocc]] = 3;
////  tmpdet[rev_reorder[nocc+1]] = 3;
////  for (auto idet : tmpdet)
////     std::cout << " " << unsigned(idet);
////  std::cout << std::endl;
////  dtrie.push_back(tmpdet);  
////  double pQ1 = parity_ab_str(tmpdet);
////  double pQ2 = parity_reorder(tmpdet, reorder, nocc);
////  double pQ3 = parity_ci_to_cc(nocc-2+nocc-1, 2, nocc);
////  double pQ4 = parity_ci_to_cc(nocc-2+nocc-1, 2, nocc);
//
//  dtrie.evaluate(zeromps);  
////  std::cout << "ref = " << dtrie.vals[0] << std::endl; 
////
////  std::cout << "D   = " << dtrie.vals[1] * pD1 * pD2 * pD3 * pD4 << std::endl; 
//
//  std::cout << "S2,0= " << dtrie.vals[0] * pS1 * pS2 * pS3 << std::endl; 
//  std::cout << pS1 << " " << pS2 << " " << pS3 << std::endl;
//
//  std::cout << "S1,0= " << dtrie.vals[1] * pD1 * pD2 * pD3 << std::endl; 
//  std::cout << pD1 << " " << pD2 << " " << pD3 << std::endl;
//
//  std::cout << "S3,5= " << dtrie.vals[2] * pT1 * pT2 * pT3 << std::endl; 
//  std::cout << pT1 << " " << pT2 << " " << pT3 << std::endl;
//
//  std::cout << "S0,5= " << dtrie.vals[3] * pQ1 * pQ2 * pQ3 << std::endl; 
//  std::cout << pQ1 << " " << pQ2 << " " << pQ3 << std::endl;
//
////  std::cout << "Q  = " << dtrie.vals[3] * pQ1 * pQ2 * pQ3 * pQ4 << std::endl; 
////============================== test ========================================  


  // ================= Determinant TRIE ====================
  MPS2CI dtrie(nsite, true); 

//  // ================= append ====================
//  // reorder [ lattice_index ] = mo_index, rev_reorder [ mo_index ] = lattice_index
//  std::vector<int> reorder = dmrginp.reorder_vector(), rev_reorder(nocc+nvir, 0);
//
//  // Ref append
//  std::vector<uint8_t> Refdet(nocc+nvir, 0), tmpdet, tmpdet2;
//  for (int i = 0; i < nocc+nvir; i++){
//      rev_reorder[reorder[i]] = i;
//      if (reorder[i]<nocc)
//          Refdet[i] = 3;  
//  }
//
//  dtrie.push_back(Refdet);  
//  parity.push_back(parity_ab_str(Refdet));
//  reorder_parity.push_back(1.0);
//  // Sa append
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      dtrie.push_back(tmpdet);  
//      parity.push_back(parity_ab_str(tmpdet));
//      reorder_parity.push_back(parity_reorder(tmpdet, reorder, nocc));
//  }
//  }
//  // Daa append
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[j]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      tmpdet[rev_reorder[b+nocc]] = 2;
//      dtrie.push_back(tmpdet);  
//      parity.push_back(parity_ab_str(tmpdet));
//      reorder_parity.push_back(parity_reorder(tmpdet, reorder, nocc));
//  }
//  }
//  }
//  }
//  // Dab append
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      for (int b = 0; b < nvir; b++){
//      for (int j = nocc-1; j > -1; j--){
//          tmpdet2.assign(tmpdet.begin(), tmpdet.end());
//          if (i != j) tmpdet2[rev_reorder[j]] = 2;
//          else  tmpdet2[rev_reorder[j]] = 0;
//          if (a != b) tmpdet2[rev_reorder[b+nocc]] = 1;
//          else  tmpdet2[rev_reorder[b+nocc]] = 3;
//          dtrie.push_back(tmpdet2);  
//          parity.push_back(parity_ab_str(tmpdet2));
//          reorder_parity.push_back(parity_reorder(tmpdet2, reorder, nocc));
//      }
//      }
//  }
//  }
//  // Taaa append
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[j]] = 1;
//      tmpdet[rev_reorder[k]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      tmpdet[rev_reorder[b+nocc]] = 2;
//      tmpdet[rev_reorder[c+nocc]] = 2;
//      dtrie.push_back(tmpdet);  
//      parity.push_back(parity_ab_str(tmpdet));
//      reorder_parity.push_back(parity_reorder(tmpdet, reorder, nocc));
//  }
//  }
//  }
//  }
//  }
//  }
//  // Taab append
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[j]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      tmpdet[rev_reorder[b+nocc]] = 2;
//      for (int c = 0; c < nvir; c++){
//      for (int k = nocc-1; k > -1; k--){
//          tmpdet2.assign(tmpdet.begin(), tmpdet.end());
//          if (i != k && j != k) tmpdet2[rev_reorder[k]] = 2;
//          else  tmpdet2[rev_reorder[k]] = 0;
//          if (a != c && b != c) tmpdet2[rev_reorder[c+nocc]] = 1;
//          else  tmpdet2[rev_reorder[c+nocc]] = 3;
//          dtrie.push_back(tmpdet2);  
//          parity.push_back(parity_ab_str(tmpdet2));
//          reorder_parity.push_back(parity_reorder(tmpdet2, reorder, nocc));
//      }
//      }
//  }
//  }
//  }
//  }
//  // Qaaab append
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[j]] = 1;
//      tmpdet[rev_reorder[k]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      tmpdet[rev_reorder[b+nocc]] = 2;
//      tmpdet[rev_reorder[c+nocc]] = 2;
//      for (int d = 0; d < nvir; d++){
//      for (int l = nocc-1; l > -1; l--){
//          tmpdet2.assign(tmpdet.begin(), tmpdet.end());
//          if (i != l && j != l && k != l) tmpdet2[rev_reorder[l]] = 2;
//          else  tmpdet2[rev_reorder[l]] = 0;
//          if (a != d && b != d && c != d) tmpdet2[rev_reorder[d+nocc]] = 1;
//          else  tmpdet2[rev_reorder[d+nocc]] = 3;
//          dtrie.push_back(tmpdet2);  
//          parity.push_back(parity_ab_str(tmpdet2));
//          reorder_parity.push_back(parity_reorder(tmpdet2, reorder, nocc));
//      }
//      }
//  }
//  }
//  }
//  }
//  }
//  }
//  // Qaabb append
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      tmpdet.assign(Refdet.begin(), Refdet.end());
//      tmpdet[rev_reorder[i]] = 1;
//      tmpdet[rev_reorder[j]] = 1;
//      tmpdet[rev_reorder[a+nocc]] = 2;
//      tmpdet[rev_reorder[b+nocc]] = 2;
//      for (int d = 1; d < nvir; d++){
//      for (int c = 0; c < d; c++){
//      for (int l = nocc-1; l > 0; l--){
//      for (int k = l-1; k > -1; k--){
//          tmpdet2.assign(tmpdet.begin(), tmpdet.end());
//          if (i != k && j != k) tmpdet2[rev_reorder[k]] = 2;
//          else  tmpdet2[rev_reorder[k]] = 0;
//          if (i != l && j != l) tmpdet2[rev_reorder[l]] = 2;
//          else  tmpdet2[rev_reorder[l]] = 0;
//          if (a != c && b != c) tmpdet2[rev_reorder[c+nocc]] = 1;
//          else  tmpdet2[rev_reorder[c+nocc]] = 3;
//          if (a != d && b != d) tmpdet2[rev_reorder[d+nocc]] = 1;
//          else  tmpdet2[rev_reorder[d+nocc]] = 3;
//          dtrie.push_back(tmpdet2);  
//          parity.push_back(parity_ab_str(tmpdet2));
//          reorder_parity.push_back(parity_reorder(tmpdet2, reorder, nocc));
//      }
//      }
//      }
//      }
//  }
//  }
//  }
//  }
//
  // ================ calculate ===================
  dtrie.evaluate(zeromps, nocc, nvir, cutoff, max_ex);

//  // ================= assign =====================
//  int ia, jb, kc, ld, ijab, klcd, ijkabc;
//  int idet = 0;
//  // Ref assign 
//  Ref[0] = dtrie.vals[idet] * reorder_parity[idet] * parity[idet];
//  // Sa assign 
//  ia = -1;
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      double parity_ia = parity_ci_to_cc(i, 1, nocc);
//      idet += 1;
//      ia   += 1;
//      S_a[ia] = parity_ia * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//  }
//  }
//  // Daa assign 
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijab = parity_ci_to_cc(i+j, 2, nocc);
//      idet += 1;
//      ijab += 1;
//      D_aa[ijab] = parity_ijab * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//  }
//  }
//  }
//  }
//  // Dab assign 
//  ia = -1;
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      double parity_ia = parity_ci_to_cc(i, 1, nocc);
//      ia   += 1;
//      jb    =-1;
//      for (int b = 0; b < nvir; b++){
//      for (int j = nocc-1; j > -1; j--){
//          double parity_jb = parity_ci_to_cc(j, 1, nocc);
//          idet += 1;
//          jb   += 1;
//          D_ab[ia][jb] = parity_ia * parity_jb * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//      }
//      }
//  }
//  }
//  // Taaa assign
//  ijkabc = -1;
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijkabc = parity_ci_to_cc(i+j+k, 3, nocc);
//      idet   += 1;
//      ijkabc += 1;
//      T_aaa[ijkabc] = parity_ijkabc * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//  }
//  }
//  }
//  }
//  }
//  }
//  // Taab assign
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijab = parity_ci_to_cc(i+j, 2, nocc);
//      ijab += 1;
//      kc = -1;
//      for (int c = 0; c < nvir; c++){
//      for (int k = nocc-1; k > -1; k--){
//          double parity_kc = parity_ci_to_cc(k, 1, nocc);
//          idet += 1; 
//          kc   += 1; 
//          T_aab[ijab][kc] = parity_ijab * parity_kc * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//      }
//      }
//  }
//  }
//  }
//  }
//  // Qaaab assign
//  ijkabc = -1;
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijkabc = parity_ci_to_cc(i+j+k, 3, nocc);
//      ijkabc += 1;
//      ld = -1;
//      for (int d = 0; d < nvir; d++){
//      for (int l = nocc-1; l > -1; l--){
//          double parity_ld = parity_ci_to_cc(l, 1, nocc);
//          idet += 1; 
//          ld   += 1; 
//          Q_aaab[ijkabc][ld] = parity_ijkabc * parity_ld * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//      }
//      }
//  }
//  }
//  }
//  }
//  }
//  }
//  // Qaabb assign
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijab = parity_ci_to_cc(i+j, 2, nocc);
//      ijab += 1;
//      klcd = -1;
//      for (int d = 1; d < nvir; d++){
//      for (int c = 0; c < d; c++){
//      for (int l = nocc-1; l > 0; l--){
//      for (int k = l-1; k > -1; k--){
//          double parity_klcd = parity_ci_to_cc(k+l, 2, nocc);
//          idet += 1; 
//          klcd += 1; 
//          Q_aabb[ijab][klcd] = parity_ijab * parity_klcd * parity[idet] * reorder_parity[idet] * dtrie.vals[idet]; 
//      }
//      }
//      }
//      }
//  }
//  }
//  }
//  }
//
  // ================= print =====================
//  double norm = 0.0;
//
//  double normRef = std::inner_product(Ref.begin(), Ref.end(), Ref.begin(), 0.0);
//  norm += normRef;
//  std::cout << "Ref norm = " << norm << std::endl;
//
//  double normS = std::inner_product(S_a.begin(), S_a.end(), S_a.begin(), 0.0);
//  norm += normS;
//  std::cout << "S   norm = " << norm << " ( " << normS << " )"  << std::endl;
//
//  double normD = std::inner_product(D_aa.begin(), D_aa.end(), D_aa.begin(), 0.0);
//  normD *= 2.0;
//  for (int i=0; i<D_ab.size(); i++)
//    normD = std::inner_product(D_ab[i].begin(), D_ab[i].end(), D_ab[i].begin(), normD);
//  norm += normD;
//  std::cout << "D   norm = " << norm  << " ( " << normD << " )"  << std::endl;
//
//  double normT = std::inner_product(T_aaa.begin(), T_aaa.end(), T_aaa.begin(), 0.0);
//  for (int i=0; i<T_aab.size(); i++)
//    normT = std::inner_product(T_aab[i].begin(), T_aab[i].end(), T_aab[i].begin(), normT);
//  normT *= 2.0;
//  norm += normT;
//  std::cout << "T   norm = " << norm  << " ( " << normT << " )"  << std::endl;
//
//  double normQ = 0.0;
//  for (int i=0; i<Q_aaab.size(); i++)
//    normQ = std::inner_product(Q_aaab[i].begin(), Q_aaab[i].end(), Q_aaab[i].begin(), normQ);
//  normQ *= 2.0;
//  for (int i=0; i<Q_aabb.size(); i++)
//    normQ = std::inner_product(Q_aabb[i].begin(), Q_aabb[i].end(), Q_aabb[i].begin(), normQ);
//  norm += normQ;
//  std::cout << "Q   norm = " << norm  << " ( " << normQ << " )"  << std::endl;
//
//  //TODO: interface by wrapping
  // ================= print =====================
//
//  //double cutoff = 1e-7;
//  std::cout << std::endl;
  std::cout << "extracted CI coeff from MPS start" << std::endl;
  std::cout << "typ,1,2,3,4,5,6,7,8,9" << std::endl;
//  std::cout << "rf,     " << Ref[0] << std::endl;
//  ia = -1;
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      ia   += 1;
//      if (std::fabs(S_a[ia]) > thresh) std::cout << "a," << i << "," << a<< ",      " << S_a[ia] << std::endl;
//  }
//  }
//  // Daa assign 
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      ijab += 1;
//      if (std::fabs(D_aa[ijab]) > thresh) std::cout << "aa," << i << ","  << j << ","
//                << a<< ","  << b<< ",      " << D_aa[ijab] << std::endl;
//  }
//  }
//  }
//  }
//  // Dab assign 
//  ia = -1;
//  for (int a = 0; a < nvir; a++){
//  for (int i = nocc-1; i > -1; i--){
//      double parity_ia = parity_ci_to_cc(i, 1, nocc);
//      ia   += 1;
//      jb    =-1;
//      for (int b = 0; b < nvir; b++){
//      for (int j = nocc-1; j > -1; j--){
//          jb   += 1;
//          if (std::fabs(D_ab[ia][jb]) > thresh)
//              std::cout << "ab," << i << ","  << a<< "," << j
//                        << ","  << b<< ",      " << D_ab[ia][jb] << std::endl;
//      }
//      }
//  }
//  }
//  // Taaa assign
//  ijkabc = -1;
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      ijkabc += 1;
//      if (std::fabs(T_aaa[ijkabc]) > thresh)
//      std::cout << "aaa," << i << ","  << j << "," << k 
//                << ","  << a<< ","  << b<< ","  << c
//                << ",      " << T_aaa[ijkabc] << std::endl;
//  }
//  }
//  }
//  }
//  }
//  }
//  // Taab assign
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijab = parity_ci_to_cc(i+j, 2, nocc);
//      ijab += 1;
//      kc = -1;
//      for (int c = 0; c < nvir; c++){
//      for (int k = nocc-1; k > -1; k--){
//          kc   += 1; 
//          if (std::fabs(T_aab[ijab][kc]) > thresh)
//          std::cout << "aab," << i << ","  << j << "," << a
//                    << ","  << b<< ","  << k << ","  << c
//                    << ",      " << T_aab[ijab][kc] << std::endl;
//      }
//      }
//  }
//  }
//  }
//  }
//  // Qaaab assign
//  ijkabc = -1;
//  for (int c = 2; c < nvir; c++){
//  for (int b = 1; b < c; b++){
//  for (int a = 0; a < b; a++){
//  for (int k = nocc-1; k > 1; k--){
//  for (int j = k-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijkabc = parity_ci_to_cc(i+j+k, 3, nocc);
//      ijkabc += 1;
//      ld = -1;
//      for (int d = 0; d < nvir; d++){
//      for (int l = nocc-1; l > -1; l--){
//          ld   += 1; 
//          if (std::fabs(Q_aaab[ijkabc][ld]) > thresh)
//          std::cout << "aaab," << i << ","  << j << "," << k 
//                    << ","  << a<< ","  << b<< ","  << c
//                    << ","  << l << ","  << d
//                    << ",      " << Q_aaab[ijkabc][ld] << std::endl;
//      }
//      }
//  }
//  }
//  }
//  }
//  }
//  }
//  // Qaabb assign
//  ijab = -1;
//  for (int b = 1; b < nvir; b++){
//  for (int a = 0; a < b; a++){
//  for (int j = nocc-1; j > 0; j--){
//  for (int i = j-1; i > -1; i--){
//      double parity_ijab = parity_ci_to_cc(i+j, 2, nocc);
//      ijab += 1;
//      klcd = -1;
//      for (int d = 1; d < nvir; d++){
//      for (int c = 0; c < d; c++){
//      for (int l = nocc-1; l > 0; l--){
//      for (int k = l-1; k > -1; k--){
//          double parity_klcd = parity_ci_to_cc(k+l, 2, nocc);
//          klcd += 1; 
//          if (std::fabs(Q_aabb[ijab][klcd]) > thresh)
//          std::cout << "aabb," << i << ","  << j << "," << a
//                    << ","  << b<< ","  << k << ","  << l 
//                    << ","  << c<< ","  << d
//                    << ",      " << Q_aabb[ijab][klcd] << std::endl;
//      }
//      }
//      }
//      }
//  }
//  }
//  }
//  }
  std::cout << "extracted CI coeff from MPS end" << std::endl;
  std::cout << std::endl;

  return;
}


int main(int argc, char* argv[])
{
  //test(argc,argv);
//  for(auto i: argv)
//    cout <<string(i)<<endl;
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  //MPI_Comm Calc;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

//  for(auto i: argv)
//    cout <<string(i)<<endl;

  ReadInput(argv[1]);
  //dmrginp.matmultFlops.resize(1, 0.);
  dmrginp.initCumulTimer();
  int max_ex = 4;
  double cutoff = dmrginp.stochasticpt_tol();

  MPS2CI_run(cutoff, max_ex);

  return 0;
}

