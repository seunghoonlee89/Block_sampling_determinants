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

bool duplicate_check_ec (partial_ex l_idx, int dot, int i_site, int nocc, int n_sites, int max_ex, const std::vector<int> reorder){
     int na_h_ex = l_idx.first.first.size(); 
     int na_p_ex = l_idx.first.second.size(); 
     int nb_h_ex = l_idx.second.first.size(); 
     int nb_p_ex = l_idx.second.second.size(); 

//     //lsh test: temporarily block for S!=0 sampling
//     // check 
//     if (reorder[i_site] < nocc){
//         if(dot == 0 || dot == 1) na_h_ex += 1; 
//         if(dot == 0 || dot == 2) nb_h_ex += 1; 
//         // # of holes can't be larger than maximum excitation
//         if(na_h_ex + nb_h_ex > max_ex) return true;
//     } 
//     else{ 
//         if(dot == 3 || dot == 2) na_p_ex += 1; 
//         if(dot == 3 || dot == 1) nb_p_ex += 1; 
//         // # of particles can't be larger than # of holes
//         //if(na_p_ex > na_h_ex) return true;
//         //if(nb_p_ex > nb_h_ex) return true;
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

struct MPS2CI;

struct MPS2CI{
    std::vector<double> vals;
    std::vector<std::pair < std::pair< std::vector<int>, std::vector<int> >, std::pair< std::vector<int>, std::vector<int> >>> idx_ex;
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
    }
    // deep copy
    std::shared_ptr<MPS2CI> copy() {
        std::shared_ptr<MPS2CI> dett =
            std::make_shared<MPS2CI>(n_sites, enable_look_up);
        dett->vals  = std::vector<double>(vals.begin(), vals.end());
        dett->idx_ex= std::vector<std::pair < std::pair< std::vector<int>, std::vector<int> >, std::pair< std::vector<int>, std::vector<int> >>> (idx_ex.begin(), idx_ex.end());
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
        std::vector<partial_wv> lwaves, rwaves;
        // initialize lwaves
        partial_wv lw;
        RowVector oldwave(1);

        //210625
        oldwave(1) = sqrt(dmrginp.molecule_quantum().get_s().getirrep()+1);
        //oldwave(1) = 1.0;

        //210625
        // tuple (num_site d, leftq_for_d, rightq_for_d+1, Ms_for_d + dotMs_for_d+1) 
        int Ms = -dmrginp.molecule_quantum().get_s().getirrep();  //This is 2*Sz of left quanta num : target 2 * Ms.
        lw[std::make_tuple(0, 0, 0, Ms)] = oldwave;
        //lw[std::make_tuple(0, 0, 0, 0)] = oldwave;

        //int Ms = 0;//This is 2*Sz of left quanta num.

        lwaves.push_back(lw);
        // initialize l_idxs 
        //std::vector<partial_ex> l_idxs, r_idxs;
        std::vector<int> vacidx ( 0, 0);
        //l_idxs.push_back( std::make_pair( std::make_pair ( vacidx, vacidx ), std::make_pair ( vacidx, vacidx )));

        int na_ex, nb_ex;

//        //lsh test
//        std::cout << "nsite: " << n_sites << std::endl;

        for (int i_site = 0; i_site < n_sites; i_site++)
        {
            int rw_len = 0;
            auto lwi = lwaves.begin();
            //auto lxi = l_idxs.begin();

            //while (lwi != lwaves.end())

//            //lsh test
//            int itest = 0; 
            //while (lwi != lwaves.end() && lxi != l_idxs.end())
            while (lwi != lwaves.end())
            { 
                auto lwave = *lwi++;
                //auto l_idx = *lxi++;
//                itest += 1; 

//                if (i_site < nocc){
//                    // l_idx: pair< pair<alpha_hole_idx, alpha_ptcl_idx>, pair<beta_hole_idx, beta_ptcl_idx> >
//                    na_ex = l_idx.first.first.size(); 
//                    nb_ex = l_idx.second.first.size(); 
//                }
//                else{
//                    na_ex = l_idx.first.second.size(); 
//                    nb_ex = l_idx.second.second.size(); 
//                }
//                //lsh test
//                std::cout << "na_ex: " << na_ex << ", nb_ex: " << nb_ex << std::endl;

                // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
                for (int dot = 0; dot < 4; dot++)
                {
                    // consider only for a, aa, ab, aaa, aab, aaab, aabb for ecCCSD
//                    bool duplicate;
//                    if (noreorder) {
//                        duplicate = duplicate_check_ec ( l_idx, dot, i_site, nocc, n_sites, max_ex );
//                    }
//                    else{
//                        duplicate = duplicate_check_ec ( l_idx, dot, i_site, nocc, n_sites, max_ex, reorder );
//                    }
//
//                    if (duplicate) continue;


//                    // update indices of configurations
//                    if (i_site < nocc)   // idx for hole
//                    {
//                        // skip for excitation
//                        if (na_ex+nb_ex == max_ex   && dot < 3)  continue;
//                        if (na_ex+nb_ex == max_ex-1 && dot == 0) continue;
//
//
////                        //lsh test
////                        if (i_site == 2){
////                            if (!(dot == 1 || dot == 3)) continue;
////                            //if (!(dot == 1)) continue;
////                        }
////                        else if (i_site != 2){
//////                            if (i_site == 0){
//////
//////                            }
//////                            if (i_site != 0){
////                                if (!(dot == 3)) continue;
//////                            }
////                        }
//
////                        //lsh test
////                        std::cout << "isite, dot: " << i_site << " " << dot << std::endl;
//                    } 
//                    else                // idx for particle 
//                    {
//                        if (na_ex+nb_ex == max_ex   && dot > 0)  continue;
//                        if (na_ex+nb_ex == max_ex-1 && dot == 3) continue;
//
////                        //lsh test
////                        if (i_site == 4){
////                            if (!(dot == 2 || dot == 0)) continue;
////                            //if (!(dot == 2)) continue;
////                        }
////                        else if (i_site != 4){
////                            if (!(dot == 0)) continue;
////                        }
//
////                        //lsh test
////                        std::cout << "isite, dot: " << i_site << " " << dot << std::endl;
//                    } 

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

//lsh test: for debugging

                    if (cutoff != 0) 
                    {
                        //l_idx: particle - hole index
                        //int ref = l_idx.first.first.size() + l_idx.first.second.size() + l_idx.second.first.size() + l_idx.second.second.size();
//                        bool protect = false;
//                        // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>
//                        //if ( ref == 0 ){
//                            if ( i_site < nelecb ){
//                               if ( dot == 3 ) protect = true;
//                            }
//                            else if ( i_site < neleca && i_site >= nelecb ){
//                               if ( dot == 2 ) protect = true;
//                            }
//                            else
//                            { 
//                               if ( dot == 0 ) protect = true;
//                            }
//                        //}

//                        if ( not protect ) {
                            for (auto &m : rwave) 
                            { 
                                sqsum += dotproduct(m.second, m.second); 
                            }
    
                            if (sqrt(sqsum) < cutoff){
                                continue;
                            }
                              //continue;
//                        }
//                        else 
//                        {
//                            //std::cout << i_site << " " << dot << " " << itest << std::endl; 
//                        }

                    }

// abandon CI coeff less than thresh
//                    if (cutoff != 0) 
//                    {
//                        int ref = l_idx.first.first.size() + l_idx.first.second.size() + l_idx.second.first.size() + l_idx.second.second.size();
//                        bool protect = false;
//                        if ( ref == 0 ){
//                            if ( reorder[i_site] < nocc ){
//                               if ( dot == 3 ) protect = true;
//                            }
//                            else
//                            { 
//                               if ( dot == 0 ) protect = true;
//                            }
//                        }
//
//                        if ( not protect ) {
//                            for (auto &m : rwave) 
//                            { 
//                                sqsum += dotproduct(m.second, m.second); 
//                            }
//    
//                            if (sqrt(sqsum) < cutoff){
//                                continue;
//                            }
//                        }
//                        else 
//                        {
//                            //std::cout << i_site << " " << dot << " " << itest << std::endl; 
//                        }
//
//                    }
                    rw_len += 1;

                    rwaves.push_back(rwave);

                    // l_idx: pair< pair<alpha_hole_idx, alpha_ptcl_idx>, pair<beta_hole_idx, beta_ptcl_idx> >
                    // dot 0: |0>, 1: |b>, 2: |a>. 3:|ab>

// index append 
//                    if (noreorder) {
//                        partial_ex r_idx(l_idx);  // copy constructor
//                        if (i_site < nocc)   // idx for hole
//                        {
//                            if (dot == 0){
//                                r_idx.first.first.push_back(i_site);
//                                r_idx.second.first.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 1){
//                                r_idx.first.first.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 2){
//                                r_idx.second.first.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            }
//                            else if (dot == 3){
//                                r_idxs.push_back(r_idx);
//                            }   
//                        } 
//                        else                // idx for particle 
//                        {
//                            if (dot == 0){
//                                r_idxs.push_back(r_idx);
//                            }   
//                            else if (dot == 1){
//                                r_idx.second.second.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 2){
//                                r_idx.first.second.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 3){
//                                r_idx.first.second.push_back(i_site);
//                                r_idx.second.second.push_back(i_site);
//                                r_idxs.push_back(r_idx);
//                            } 
//                        } 
//                    }
//                    else {
//                        partial_ex r_idx(l_idx);  // copy constructor
//                        int i_mo = reorder[i_site];
//                        if (i_mo < nocc)   // idx for hole
//                        {
//                            if (dot == 0){
//                                r_idx.first.first.push_back(i_mo);
//                                r_idx.second.first.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 1){
//                                r_idx.first.first.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 2){
//                                r_idx.second.first.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            }
//                            else if (dot == 3){
//                                r_idxs.push_back(r_idx);
//                            }   
//                        } 
//                        else                // idx for particle 
//                        {
//                            if (dot == 0){
//                                r_idxs.push_back(r_idx);
//                            }   
//                            else if (dot == 1){
//                                r_idx.second.second.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 2){
//                                r_idx.first.second.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            } 
//                            else if (dot == 3){
//                                r_idx.first.second.push_back(i_mo);
//                                r_idx.second.second.push_back(i_mo);
//                                r_idxs.push_back(r_idx);
//                            } 
//                        } 
//                    }



//                    //lsh test
//                    std::cout << "isite: " << i_site << ", dot: " << dot << ", rw_len: " << rw_len << ", sqsum: " << sqsum << std::endl; 



                }

            }

            if (i_site < n_sites-1){ 
//                //lsh test
//                //std::cout << "rwaves size: " << rwaves.size() << ", r_idxs size: " << r_idxs.size() << std::endl;
//                std::cout << "rwaves size: " << rwaves.size() << std::endl;
                auto lwi = lwaves.begin();
                //auto lxi = l_idxs.begin();
                //while (lwi != lwaves.end())
                //while (lwi != lwaves.end() && lxi != l_idxs.end())
                while (lwi != lwaves.end())
                { 
                    auto lwave = *lwi++;
                    //auto l_idx = *lxi++;
                    lwave.clear();
                    //l_idx.first.first.clear();
                    //l_idx.first.second.clear();
                    //l_idx.second.first.clear();
                    //l_idx.second.second.clear();
                } 
                lwaves.clear();
                //l_idxs.clear();
                lwaves = rwaves;
                //l_idxs = r_idxs;


////            //lsh test
//            std::cout << "ith site: " << i_site << " rwaves size: " << rwaves.size() << std::endl;
//            //int test = 0;

                auto rwi = rwaves.begin();
                //auto rxi = r_idxs.begin();
                //while (rwi != rwaves.end())
                //while (rwi != rwaves.end() && rxi != r_idxs.end())
                while (rwi != rwaves.end())
                { 
										auto rwave = *rwi++;
                    //auto r_idx = *rxi++;


//            int extest = r_idx.first.first.size() + r_idx.first.second.size() + r_idx.second.first.size() + r_idx.second.second.size();
//            test += 1;
//            std::cout << test << " " << extest << std::endl; 

                    rwave.clear();
//                    r_idx.first.first.clear();
//                    r_idx.first.second.clear();
//                    r_idx.second.first.clear();
//                    r_idx.second.second.clear();
                } 
                rwaves.clear();
//                r_idxs.clear();
            }
            else{ 
                vals.reserve(rw_len);

                int tmp = 0;
                for (auto &rwave : rwaves) { 
//                    std::cout << "test1" << std::endl;
                    for (auto &m : rwave) { 
//                        std::cout << "test2" << std::endl;

                        assert(m.second.Ncols() == 1);
                        vals.push_back(m.second.element(0));
                        tmp += 1; 
                    }
                }
                assert(tmp == 1);
//                idx_ex = r_idxs;

//            //lsh test
//            std::cout << "ith site: " << i_site << " rwaves size: " << rwaves.size() << " vals size:" << vals.size() << std::endl;
//            int test = 0;
//
//                auto rwi = rwaves.begin();
//                auto rxi = r_idxs.begin();
//                //while (rwi != rwaves.end())
//                while (rwi != rwaves.end() && rxi != r_idxs.end())
//                { 
//										auto rwave = *rwi++;
//                    auto r_idx = *rxi++;
//
//
//            int extest = r_idx.first.first.size() + r_idx.first.second.size() + r_idx.second.first.size() + r_idx.second.second.size();
//            test += 1;
//            std::cout << test << " " << extest << std::endl; 
//
//                    rwave.clear();
//                    r_idx.first.first.clear();
//                    r_idx.first.second.clear();
//                    r_idx.second.first.clear();
//                    r_idx.second.second.clear();
//                } 

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
  //int nocc  = dmrginp.total_particle_number()/2;
  //int nvir  = nsite - nocc; 
  int Ms2   = dmrginp.molecule_quantum().get_s().getirrep();  //This is 2*Sz (n_alpha - n_beta) 
  int neleca= dmrginp.alpha_particle_number();
  int nelecb= dmrginp.beta_particle_number();
  int nelec = neleca+nelecb;
  //int nelec = dmrginp.total_particle_number();
  //int neleca= (nelec + Ms2)/2;
  //int nelecb= (nelec - Ms2)/2;
  std::cout << "nelec tot, alpha, beta =" << nelec << ", " << neleca << ", " << nelecb << std::endl;

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
  std::cout << "noreorder: " << noreorder << std::endl;
  // ================= parity and print =====================
  std::cout << "extracted CI coeff from MPS" << std::endl;
  std::cout << "total " << dtrie.vals.size() << " CI coeffs are sampled" << std::endl;
  //std::cout << "typ,1,2,3,4,5,6,7,8,9" << std::endl;

  for (int it = 0; it < dtrie.vals.size(); it++){ 
//      int na_h_ex = dtrie.idx_ex[it].first.first.size(); 
//      int na_p_ex = dtrie.idx_ex[it].first.second.size(); 
//      int nb_h_ex = dtrie.idx_ex[it].second.first.size(); 
//      int nb_p_ex = dtrie.idx_ex[it].second.second.size(); 
//      if (na_h_ex == 0 && nb_h_ex == 0){

//         std::cout << "rf,"; 
//      }
      std::cout << "       " << dtrie.vals[it] << std::endl; 
  }

//  for (int it = 0; it < dtrie.vals.size(); it++){ 
//      int na_h_ex = dtrie.idx_ex[it].first.first.size(); 
//      int na_p_ex = dtrie.idx_ex[it].first.second.size(); 
//      int nb_h_ex = dtrie.idx_ex[it].second.first.size(); 
//      int nb_p_ex = dtrie.idx_ex[it].second.second.size(); 
//
//      assert(na_h_ex == na_p_ex);
//      assert(nb_h_ex == nb_p_ex);
//
//      if (na_h_ex == 0 && nb_h_ex == 0){
//
//         parity = parity_ab_str(Refdet);
//         std::cout << "rf,"; 
//      }
//      else if (na_h_ex == 1 && nb_h_ex == 0){
//
//         int i = dtrie.idx_ex[it].first.first[0];
//         int a = dtrie.idx_ex[it].first.second[0];
//          
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//         std::cout << "a," << i << "," << a << ","; 
//      }
//      else if (na_h_ex == 2 && nb_h_ex == 0){
//
//         int i = dtrie.idx_ex[it].first.first[0];
//         int j = dtrie.idx_ex[it].first.first[1];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int b = dtrie.idx_ex[it].first.second[1];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[j]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         tmpdet[rev_reorder[b]] = 2;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::sort(dtrie.idx_ex[it].first.first.begin(), dtrie.idx_ex[it].first.first.end());
//         std::sort(dtrie.idx_ex[it].first.second.begin(),dtrie.idx_ex[it].first.second.end());
//
//         std::cout << "aa,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//      }
//      else if (na_h_ex == 1 && nb_h_ex == 1){
//
//         int i = dtrie.idx_ex[it].first.first[0];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int j = dtrie.idx_ex[it].second.first[0];
//         int b = dtrie.idx_ex[it].second.second[0];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         if (i != j) tmpdet[rev_reorder[j]] = 2;
//         else  tmpdet[rev_reorder[j]] = 0;
//         if (a != b) tmpdet[rev_reorder[b]] = 1;
//         else  tmpdet[rev_reorder[b]] = 3;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::cout << "ab,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.second)
//             std::cout << idx << ",";
//      }
//      else if (na_h_ex == 3 && nb_h_ex == 0){
//
//         int i = dtrie.idx_ex[it].first.first[0];
//         int j = dtrie.idx_ex[it].first.first[1];
//         int k = dtrie.idx_ex[it].first.first[2];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int b = dtrie.idx_ex[it].first.second[1];
//         int c = dtrie.idx_ex[it].first.second[2];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[j]] = 1;
//         tmpdet[rev_reorder[k]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         tmpdet[rev_reorder[b]] = 2;
//         tmpdet[rev_reorder[c]] = 2;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::sort(dtrie.idx_ex[it].first.first.begin(), dtrie.idx_ex[it].first.first.end());
//         std::sort(dtrie.idx_ex[it].first.second.begin(),dtrie.idx_ex[it].first.second.end());
//
//         std::cout << "aaa,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//      }
//      else if (na_h_ex == 2 && nb_h_ex == 1){
//         int i = dtrie.idx_ex[it].first.first[0];
//         int j = dtrie.idx_ex[it].first.first[1];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int b = dtrie.idx_ex[it].first.second[1];
//         int k = dtrie.idx_ex[it].second.first[0];
//         int c = dtrie.idx_ex[it].second.second[0];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[j]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         tmpdet[rev_reorder[b]] = 2;
//         if (i != k && j != k) tmpdet[rev_reorder[k]] = 2;
//         else  tmpdet[rev_reorder[k]] = 0;
//         if (a != c && b != c) tmpdet[rev_reorder[c]] = 1;
//         else  tmpdet[rev_reorder[c]] = 3;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::sort(dtrie.idx_ex[it].first.first.begin(), dtrie.idx_ex[it].first.first.end());
//         std::sort(dtrie.idx_ex[it].first.second.begin(),dtrie.idx_ex[it].first.second.end());
//
//         std::cout << "aab,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.second)
//             std::cout << idx << ",";
//      }
//      else if (na_h_ex == 3 && nb_h_ex == 1){
//         int i = dtrie.idx_ex[it].first.first[0];
//         int j = dtrie.idx_ex[it].first.first[1];
//         int k = dtrie.idx_ex[it].first.first[2];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int b = dtrie.idx_ex[it].first.second[1];
//         int c = dtrie.idx_ex[it].first.second[2];
//         int l = dtrie.idx_ex[it].second.first[0];
//         int d = dtrie.idx_ex[it].second.second[0];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[j]] = 1;
//         tmpdet[rev_reorder[k]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         tmpdet[rev_reorder[b]] = 2;
//         tmpdet[rev_reorder[c]] = 2;
//         if (i != l && j != l && k != l) tmpdet[rev_reorder[l]] = 2;
//         else  tmpdet[rev_reorder[l]] = 0;
//         if (a != d && b != d && c != d) tmpdet[rev_reorder[d]] = 1;
//         else  tmpdet[rev_reorder[d]] = 3;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::sort(dtrie.idx_ex[it].first.first.begin(), dtrie.idx_ex[it].first.first.end());
//         std::sort(dtrie.idx_ex[it].first.second.begin(),dtrie.idx_ex[it].first.second.end());
//
//         std::cout << "aaab,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.second)
//             std::cout << idx << ",";
//      }
//      else if (na_h_ex == 2 && nb_h_ex == 2){
//         int i = dtrie.idx_ex[it].first.first[0];
//         int j = dtrie.idx_ex[it].first.first[1];
//         int a = dtrie.idx_ex[it].first.second[0];
//         int b = dtrie.idx_ex[it].first.second[1];
//         int k = dtrie.idx_ex[it].second.first[0];
//         int l = dtrie.idx_ex[it].second.first[1];
//         int c = dtrie.idx_ex[it].second.second[0];
//         int d = dtrie.idx_ex[it].second.second[1];
//
//         tmpdet.assign(Refdet.begin(), Refdet.end());
//         tmpdet[rev_reorder[i]] = 1;
//         tmpdet[rev_reorder[j]] = 1;
//         tmpdet[rev_reorder[a]] = 2;
//         tmpdet[rev_reorder[b]] = 2;
//         if (i != k && j != k) tmpdet[rev_reorder[k]] = 2;
//         else  tmpdet[rev_reorder[k]] = 0;
//         if (i != l && j != l) tmpdet[rev_reorder[l]] = 2;
//         else  tmpdet[rev_reorder[l]] = 0;
//         if (a != c && b != c) tmpdet[rev_reorder[c]] = 1;
//         else  tmpdet[rev_reorder[c]] = 3;
//         if (a != d && b != d) tmpdet[rev_reorder[d]] = 1;
//         else  tmpdet[rev_reorder[d]] = 3;
//         parity  = parity_ab_str(tmpdet);
//         parity *= parity_reorder(tmpdet, reorder, nocc);
//
//         std::sort(dtrie.idx_ex[it].first.first.begin(), dtrie.idx_ex[it].first.first.end());
//         std::sort(dtrie.idx_ex[it].first.second.begin(),dtrie.idx_ex[it].first.second.end());
//         std::sort(dtrie.idx_ex[it].second.first.begin(), dtrie.idx_ex[it].second.first.end());
//         std::sort(dtrie.idx_ex[it].second.second.begin(),dtrie.idx_ex[it].second.second.end());
//
//         std::cout << "aabb,";
//         for (auto idx : dtrie.idx_ex[it].first.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].first.second)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.first)
//             std::cout << idx << ",";
//         for (auto idx : dtrie.idx_ex[it].second.second)
//             std::cout << idx << ",";
//      }
//      else{
//         std::cout << "what is it?" << std::endl;
//         continue;
//      }
//      std::cout << "       " << parity * dtrie.vals[it] << std::endl; 
//  }




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
  int max_ex = 10;
  double cutoff = dmrginp.stochasticpt_tol();

  MPS2CI_run(cutoff, max_ex);

  return 0;
}

