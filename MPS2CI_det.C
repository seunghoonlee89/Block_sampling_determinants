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
#include "MPS2CI.h"
#include <math.h>
#include <numeric> 
#include <boost/iterator/zip_iterator.hpp>

// 210625: revise for S!=0 determinant sampling

using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;
int bitstring::n_orb;
void ReadInput(char* conf);

// consider only for rf, a, aa, ab, aaa, aab, aaab, aabb
// dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
bool duplicate_check_ec(partial_ex l_idx, int dot, int i_site, int nocc, int n_sites, int max_ex){
     int na_h_ex = l_idx.first.first.size(); 
     int na_p_ex = l_idx.first.second.size(); 
     int nb_h_ex = l_idx.second.first.size(); 
     int nb_p_ex = l_idx.second.second.size(); 

//     //lsh test: temporarily block for S!=0 sampling
//     // check 
//     if (i_site < nocc){
//         if(dot == 0 || dot == 1) na_h_ex += 1; 
//         if(dot == 0 || dot == 2) nb_h_ex += 1; 
//         // # of holes can't be larger than maximum excitation
//         if(na_h_ex + nb_h_ex > max_ex) return true;
//     } 
//     else{ 
//         if(dot == 3 || dot == 2) na_p_ex += 1; 
//         if(dot == 3 || dot == 1) nb_p_ex += 1; 
//         // # of particles can't be larger than # of holes
//         if(na_p_ex > na_h_ex) return true;
//         if(nb_p_ex > nb_h_ex) return true;
//     }
// 
//     // a, aa, ab, aaa, aab, aaab, aabb
//     if (na_h_ex > 3) return true;
//     if (nb_h_ex > 2) return true;
//
//     if (i_site == n_sites-1){
//        if (!( (na_h_ex==0 && nb_h_ex==0) || (na_h_ex==1 && nb_h_ex==0) || (na_h_ex==2 && nb_h_ex==0) ||
//               (na_h_ex==1 && nb_h_ex==1) || (na_h_ex==3 && nb_h_ex==0) || (na_h_ex==2 && nb_h_ex==1) ||  
//               (na_h_ex==3 && nb_h_ex==1) || (na_h_ex==2 && nb_h_ex==2) ) ) return true;
//     } 
//     if (i_site == n_sites-2){
//        if ( (na_h_ex==0 && nb_h_ex==2) ) return true;
//        
//        // <==> if (!( (na_h_ex==0 && nb_h_ex==0) || (na_h_ex==1 && nb_h_ex==0) || (na_h_ex==2 && nb_h_ex==0) ||
//        //             (na_h_ex==1 && nb_h_ex==1) || (na_h_ex==0 && nb_h_ex==1) || (na_h_ex==3 && nb_h_ex==0) ||
//        //             (na_h_ex==2 && nb_h_ex==1) || (na_h_ex==3 && nb_h_ex==1) || (na_h_ex==2 && nb_h_ex==2) ||
//        //             (na_h_ex==1 && nb_h_ex==2) ) ) return true;
//     } 

     return false;
}

bool skip_check (partial_occ l_occ, int dot, int i_site, int n_sites, int neleca, int nelecb){
    // check particle qn
    int neleca_l = 0, nelecb_l = 0;
    for (char& occ: l_occ){
        if ( occ == '2'){
            neleca_l += 1;
            nelecb_l += 1;
        }
        else if (occ == 'a'){
            neleca_l += 1;
        }
        else if (occ == 'b'){
            nelecb_l += 1;
        }
    }
    // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
    if (dot == 3){
        neleca_l += 1;
        nelecb_l += 1;
    }
    else if (dot == 2){
        neleca_l += 1;
    }
    else if (dot == 1){
        nelecb_l += 1;
    }

    //std::cout << neleca_l << " / " << neleca << " " << nelecb_l << " / " << nelecb << std::endl;
    if (neleca_l > neleca) return true; 
    if (nelecb_l > nelecb) return true; 

    // check Ms 



    return false;

}

struct MPS2CI;

struct MPS2CI{
    std::vector<double> vals;
    std::vector<std::pair < std::pair< std::vector<int>, std::vector<int> >, std::pair< std::vector<int>, std::vector<int> >>> idx_ex;
    std::vector<std::string> occ_pattern;
    int n_sites;
    bool enable_look_up;

    MPS2CI(int n_sites, bool enable_look_up = false)
        : n_sites(n_sites), enable_look_up(enable_look_up) {
    }

//lsh need to change
    // clear trie
    void clear() {
        vals.clear();
        for (auto &idx : idx_ex){
            idx.first.first.clear();
            idx.first.second.clear();
            idx.second.first.clear();
            idx.second.second.clear();
        }
        idx_ex.clear();
        for (auto &occ : occ_pattern){
            occ.clear();
        }
        occ_pattern.clear();
    }
    // deep copy
    std::shared_ptr<MPS2CI> copy() {
        std::shared_ptr<MPS2CI> dett =
            std::make_shared<MPS2CI>(n_sites, enable_look_up);
        dett->vals  = std::vector<double>(vals.begin(), vals.end());
        dett->idx_ex= std::vector<std::pair < std::pair< std::vector<int>, std::vector<int> >, std::pair< std::vector<int>, std::vector<int> >>> (idx_ex.begin(), idx_ex.end());
        dett->occ_pattern= std::vector<std::string> (occ_pattern.begin(), occ_pattern.end());
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
    void evaluate(const std::shared_ptr<simplemps> &mps, const std::vector<int> reorder, int nsite, int neleca, int nelecb, bool noreorder, double cutoff = 1e-10, int max_ex = 4) {
        double timer_i = 0.0;
        double timer_f = 0.0;

        //mpi
        boost::mpi::communicator world;
        int nprocs = world.size();
        int iproc  = world.rank();
        int partition;
        std::vector<partial_wv> lwaves, rwaves;
        std::vector<partial_occ> l_occs, r_occs;

        if(iproc == 0) {
            // initialize lwaves
            partial_wv lw;
            RowVector oldwave(1);
    
            //210625
            //In singlet embedding, there is an ensemble of 2*s+1 state with different Sz. The norm of each state is 1/sqrt(2*s+1);
            oldwave(1) = sqrt(dmrginp.molecule_quantum().get_s().getirrep()+1);
            //oldwave(1) = 1.0;
    
            //210625
            // tuple (num_site d, leftq_for_d, rightq_for_d+1, Ms_for_d + dotMs_for_d+1) 
            int Ms2 = -dmrginp.molecule_quantum().get_s().getirrep();  //This is 2*Sz of left quanta num : target 2 * Ms.
            lw[std::make_tuple(0, 0, 0, Ms2)] = oldwave;
            //lw[std::make_tuple(0, 0, 0, 0)] = oldwave;
    
            lwaves.push_back(lw);
            // initialize l_partial_occ
            l_occs.push_back( "" );

        }
//        //lsh test
//        std::cout << "nsite: " << n_sites << std::endl;

        for (int i_site = 0; i_site < n_sites; i_site++)
        {
            int rw_len = 0;
            // broadcast lwaves
            if(iproc == 0) 
            {
              std::cout << " Sampling at " << i_site << " / " << n_sites << ". # of partials: " << lwaves.size() << ". Time " << (timer_f - timer_i) / (double) CLOCKS_PER_SEC << " (s). " << std::endl;
              timer_i = std::clock();

              //std::cout << " test " << iproc << std::endl;
              partition = (int) (lwaves.size() / nprocs);
              for (int i=1; i<nprocs; i++){
                  world.send(i, 4, partition);
              }
              int partition_master;

              if (partition > 0){
                  for (int i=1; i<nprocs-1; i++){
                      //std::cout << " test " << iproc << " to " << i << std::endl;
                      std::vector<partial_wv> lwaves_tmp(lwaves.cbegin()+i*partition, lwaves.cbegin()+(i+1)*partition);
                      world.send(i, 0, lwaves_tmp);
                      std::vector<partial_occ> loccs_tmp(l_occs.cbegin()+i*partition, l_occs.cbegin()+(i+1)*partition);
                      world.send(i, 1, loccs_tmp);
                  }
                  {
                      int i = nprocs-1;
                      //std::cout << " test " << iproc << " to " << i << std::endl;
                      std::vector<partial_wv> lwaves_tmp(lwaves.cbegin()+i*partition, lwaves.cend());
                      world.send(i, 0, lwaves_tmp);
                      std::vector<partial_occ> loccs_tmp(l_occs.cbegin()+i*partition, l_occs.cend());
                      world.send(i, 1, loccs_tmp);
                  }
                  partition_master = partition;
              }
              else{
                  partition_master = lwaves.size();
              }
              //std::cout << " send done " << iproc << std::endl;

              auto lwi = lwaves.begin();
              auto lxi = l_occs.begin();
              assert (lwaves.size() == l_occs.size());
              while (lwi != lwaves.begin()+partition_master)
              { 
                  auto lwave = *lwi++;
                  auto l_occ = *lxi++;
                  // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                  for (int dot = 0; dot < 4; dot++)
                  {
                      // check nptcl, Ms
                      // TODO: add time-reversal symm 
//                      bool skip_symm;
//                      skip_symm = skip_check ( l_occ, dot, i_site, n_sites, neleca, nelecb );
//  
//                      if (skip_symm) continue;
  
              //   RowVector(lwave) SiteMatrices                     RowVector(rwave)
              //                        dot
              //                         |                        =  
              // leftq --*-- leftq .   --*-- rightq = LefttoRight      --*--
              //      i_site          i_site+1                        i_site+1
              //input:  i_site, dot, leftq, Ms, dotMs,
              //output: rwave<SQ,RowVector>
                      partial_wv rwave;
                      mps->lw_dot_mat_new(i_site, dot, lwave, rwave);
  
                      double sqsum = 0;
                      if (cutoff != 0) 
                      {
                          for (auto &m : rwave) 
                          { 
                              sqsum += dotproduct(m.second, m.second); 
                          }
  
                          if (sqrt(sqsum) < cutoff){
                              continue;
                          }
                      }
                      rw_len += 1;
                      rwaves.push_back(rwave);
  
                      // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                      if (noreorder) {
                          partial_occ r_occ(l_occ);  // copy constructor
                          if (dot == 0) r_occ.append("0");
                          if (dot == 1) r_occ.append("b");
                          if (dot == 2) r_occ.append("a");
                          if (dot == 3) r_occ.append("2");
                          r_occs.push_back(r_occ);
                      }
                  }
              }

            }
            else {
              world.recv(0, 4, partition);
              if (partition == 0) {continue;}

              world.recv(0, 0, lwaves);
              world.recv(0, 1, l_occs);
              //std::cout << " test " << iproc << std::endl;

              //lsh test
              //std::cout << "Received numbers of lwaves for rank " << iproc << ": " << lwaves.size() << std::endl;

              auto lwi = lwaves.begin();
              auto lxi = l_occs.begin();
              assert (lwaves.size() == l_occs.size());
  
              //while (lwi != lwaves.end())
              while (lwi != lwaves.end() && lxi != l_occs.end())
              { 
                  auto lwave = *lwi++;
                  auto l_occ = *lxi++;
                  // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                  for (int dot = 0; dot < 4; dot++)
                  {
                      // check nptcl, Ms
                      // TODO: add time-reversal symm 
//                      bool skip_symm;
//                      skip_symm = skip_check ( l_occ, dot, i_site, n_sites, neleca, nelecb );
//  
//                      if (skip_symm) continue;
  
              //   RowVector(lwave) SiteMatrices                     RowVector(rwave)
              //                        dot
              //                         |                        =  
              // leftq --*-- leftq .   --*-- rightq = LefttoRight      --*--
              //      i_site          i_site+1                        i_site+1
              //input:  i_site, dot, leftq, Ms, dotMs,
              //output: rwave<SQ,RowVector>
                      partial_wv rwave;
                      mps->lw_dot_mat_new(i_site, dot, lwave, rwave);
  
                      double sqsum = 0;
                      if (cutoff != 0) 
                      {
                          for (auto &m : rwave) 
                          { 
                              sqsum += dotproduct(m.second, m.second); 
                          }
  
                          if (sqrt(sqsum) < cutoff){
                              continue;
                          }
                      }
                      rw_len += 1;
                      rwaves.push_back(rwave);
  
                      // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                      if (noreorder) {
                          partial_occ r_occ(l_occ);  // copy constructor
                          if (dot == 0) r_occ.append("0");
                          if (dot == 1) r_occ.append("b");
                          if (dot == 2) r_occ.append("a");
                          if (dot == 3) r_occ.append("2");
                          r_occs.push_back(r_occ);
                      }
                  }
  
              }

            }

            if(iproc != 0) 
            {
                assert(partition != 0);
                if (i_site < n_sites-1){ 
                    auto lwi = lwaves.begin();
                    auto lxi = l_occs.begin();
                    //while (lwi != lwaves.end())
                    while (lwi != lwaves.end() && lxi != l_occs.end())
                    { 
                        auto lwave = *lwi++;
                        auto l_occ = *lxi++;
                        lwave.clear();
                        l_occ.clear();
                    } 
                    lwaves.clear();
                    l_occs.clear();

                    world.send(0, 2, rwaves);
                    world.send(0, 3, r_occs);
                    auto rwi = rwaves.begin();
                    auto rxi = r_occs.begin();
                    //while (rwi != rwaves.end())
                    while (rwi != rwaves.end() && rxi != r_occs.end())
                    { 
    										auto rwave = *rwi++;
                        auto r_occ = *rxi++;
    
                        rwave.clear();
                        r_occ.clear();
                    } 
                    rwaves.clear();
                    r_occs.clear();
                }
                else{ 
                    std::vector<double> vals_slav;
                    vals_slav.reserve(rw_len);
                    int tmp = 0;
                    for (auto &rwave : rwaves) { 
                        for (auto &m : rwave) { 
                            assert(m.second.Ncols() == 1);
                            vals_slav.push_back(m.second.element(0));
                            tmp += 1; 
                        }
                    }
                    assert(tmp == 1); //?

                    world.send(0, 2, vals_slav);
                    world.send(0, 3, r_occs);
                }
            } 
            else 
            {
                if (i_site < n_sites-1){ 
                    auto lwi = lwaves.begin();
                    auto lxi = l_occs.begin();
                    //while (lwi != lwaves.end())
                    while (lwi != lwaves.end() && lxi != l_occs.end())
                    { 
                        auto lwave = *lwi++;
                        auto l_occ = *lxi++;
                        lwave.clear();
                        l_occ.clear();
                    } 
                    lwaves.clear();
                    l_occs.clear();
                    lwaves = rwaves;
                    l_occs = r_occs;
    
                    auto rwi = rwaves.begin();
                    auto rxi = r_occs.begin();
                    //while (rwi != rwaves.end())
                    while (rwi != rwaves.end() && rxi != r_occs.end())
                    { 
    										auto rwave = *rwi++;
                        auto r_occ = *rxi++;
    
                        rwave.clear();
                        r_occ.clear();
                    } 
                    rwaves.clear();
                    r_occs.clear();

                    if (partition > 0){
                       for (int i=1; i<nprocs; i++){
                           std::vector<partial_wv> lwaves_tmp;
                           world.recv(i, 2, lwaves_tmp);
                           lwaves.insert( lwaves.end(), std::make_move_iterator(lwaves_tmp.begin()),
                                          std::make_move_iterator(lwaves_tmp.end()) );
                       }
                       for (int i=1; i<nprocs; i++){
                           std::vector<partial_occ> l_occs_tmp;
                           world.recv(i, 3, l_occs_tmp);
                           l_occs.insert( l_occs.end(), std::make_move_iterator(l_occs_tmp.begin()),
                                          std::make_move_iterator(l_occs_tmp.end()) );
                       }
                    } 

                }

                else{ 
                    vals.reserve(rw_len);
    
                    int tmp = 0;
                    for (auto &rwave : rwaves) { 
                        for (auto &m : rwave) { 
    
                            assert(m.second.Ncols() == 1);
                            vals.push_back(m.second.element(0));
                            tmp += 1; 
                        }
                    }
                    assert(tmp == 1);

                    if (partition > 0){
                       for (int i=1; i<nprocs; i++){
                           std::vector<double> vals_slav;
                           world.recv(i, 2, vals_slav);
                           vals.insert( vals.end(), std::make_move_iterator(vals_slav.begin()),
                                        std::make_move_iterator(vals_slav.end()) );
                       }
                       for (int i=1; i<nprocs; i++){
                           std::vector<partial_occ> r_occs_tmp;
                           world.recv(i, 3, r_occs_tmp);
                           r_occs.insert( r_occs.end(), std::make_move_iterator(r_occs_tmp.begin()),
                                          std::make_move_iterator(r_occs_tmp.end()) );
                       }
                    } 

                    occ_pattern = r_occs;
                }

              timer_f = std::clock();

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
  boost::mpi::communicator world;
  std::shared_ptr<simplemps> zeromps;
  if(dmrginp.spinAdapted())
  {
    std::shared_ptr<NonAbelianmps> mps = std::make_shared<NonAbelianmps>();
    //Use their own shared_ptr for serialization in boost broadcast;
    if (dmrginp.targetState() == -1){
        mps->build(0);
    }
    else{
        mps->build(dmrginp.targetState());
    }
    zeromps = mps;
  }
  else{
    assert(False);
  }
  int nsite = dmrginp.last_site(); // n orb
  //int nocc  = dmrginp.total_particle_number()/2;
  //int nvir  = nsite - nocc; 
  int Ms2   = dmrginp.molecule_quantum().get_s().getirrep();  //This is 2*Sz (n_alpha - n_beta) 
  int neleca= dmrginp.alpha_particle_number();
  int nelecb= dmrginp.beta_particle_number();
  int nelec = neleca+nelecb;
  //int nelec = dmrginp.total_particle_number();
  //int neleca= (nelec + Ms2)/2;
  //int nelecb= (nelec - Ms2)/2;

  if (world.rank() == 0) {
  std::cout << "nelec tot, alpha, beta =" << nelec << ", " << neleca << ", " << nelecb << std::endl;
  }

  // reorder [ lattice_index ] = mo_index, rev_reorder [ mo_index ] = lattice_index
  double parity;
  std::vector<int> reorder = dmrginp.reorder_vector(), rev_reorder(nsite, 0);

  // det[i] = 0 (empty) 1 (beta) 2 (alpha) 3 (alpha beta)

  // first no consider reorder
  // Ref append
  bool noreorder = true;
  std::vector<uint8_t> Refdet(nsite, 0), tmpdet;
  assert(neleca > nelecb);
  for (int i = 0; i < nsite; i++){
      if (i<nelecb){
          Refdet[i] = 3;  
      }
      else if (i<neleca && i>=nelecb){
          Refdet[i] = 2;  
      }
      else if (i>=neleca){
          Refdet[i] = 0;  
      }
      else{
          assert(.False.);
      }
  }

//  bool noreorder = true;
//  std::vector<uint8_t> Refdet(nsite, 0), tmpdet;
//  for (int i = 0; i < nsite; i++){
//      rev_reorder[reorder[i]] = i;
//      if (reorder[i]<neleca && reorder[i]<nelecb){
//          Refdet[i] = 3;  
//          if (i>=nocc) noreorder = false; 
//
//      }
//  }

  MPS2CI dtrie(nsite, true);
   
  dtrie.evaluate(zeromps, reorder, nsite, neleca, nelecb, noreorder, cutoff, max_ex);
  if (world.rank() == 0) 
  {
      std::cout << "noreorder: " << noreorder << std::endl;
      // ================= parity and print =====================
      std::cout << "extracted CI coeff from MPS" << std::endl;
      std::cout << "total " << dtrie.vals.size() << " CI coeffs are sampled" << std::endl;
      //std::cout << "typ,1,2,3,4,5,6,7,8,9" << std::endl;
    
      for (auto it: sort_indices(dtrie.vals)){
    
          std::cout << dtrie.occ_pattern[it]; 
          std::cout << "       " << dtrie.vals[it] << std::endl; 
      }
    
      std::cout << "extracted CI coeff from MPS end" << std::endl;
      std::cout << std::endl;

  }
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
  int max_ex = 10;
  double cutoff = dmrginp.stochasticpt_tol();

  //lsh test
  MPS2CI_run(cutoff, max_ex);

//  //lsh test : mpi
//  if(world.rank() == 0) {
//    //std::vector<partial_wv> lwaves;
//    partial_wv lw;
//    RowVector test(1);
//    test(1) = 7;
//    SQ a (0, 0, 0, 7);
//    //a=std::make_tuple(); 
//    lw[a] = test;
//    //Matrix test(1,1);
//    std::cout << "rank " << world.rank() << " throw" << lw[a].element(0) << std::endl;
//
//    {
//    //std::vector<int> test_tmp(test.cbegin(), test.cbegin()+2);
//    world.send(1, 0, lw);
//    }
////    {
////    std::vector<int> test_tmp(test.cbegin()+2, test.cbegin()+4);
////    world.send(2, 0, test_tmp);
////    }
//
////    for (int i=1; i<3; i++){
////        world.send(i, i, test);
////    }
//  }
//  else {
//    //std::vector<int> test;
//    partial_wv lw;
//    //RowVector test;
//    SQ a (0, 0, 0, 7);
//    //std::tuple a;
//    //a=std::make_tuple(0, 0, 0, 7); 
//    //Matrix test;
//    world.recv(0, 0, lw);
//    std::cout << "Received GPS positions:" << std::endl;
//    std::cout << "rank " << world.rank() << " received " << lw[a].element(0) << std::endl;
//
////    for(int i=0;i<test.size(); i++) {
////      std::cout << "rank " << world.rank() << " received " << test[i] << std::endl;
////    }
//    }


  return 0;
}

