#include "MatrixBLAS.h"
#include "rotationmat.h"
#include "wavefunction.h"
#include "global.h"
#include "sweep.h"
#include "nonspinmps.h"
#include <random>
#include <ctime>
#include <chrono>
#include "couplingCoeffs.h"
#include <unordered_map>

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;

void Abelianmps::convertblockfiles()
{
    physicaldim =dmrginp.spinAdapted()?3:4; 
    num_spatial_orbs = dmrginp.spinAdapted()?dmrginp.last_site():dmrginp.last_site()/2;
    SiteMatrices.resize(num_spatial_orbs);
    Leftindex.resize(num_spatial_orbs);
    Rightindex.resize(num_spatial_orbs);
    LefttoRight.resize(num_spatial_orbs);
    RighttoLeft.resize(num_spatial_orbs);
    for(int i=0;i<num_spatial_orbs;i++)
    {
      SiteMatrices[i].resize(physicaldim);
      Leftindex[i].resize(physicaldim);
      Rightindex[i].resize(physicaldim);
      LefttoRight[i].resize(physicaldim);
      RighttoLeft[i].resize(physicaldim);
    }
    for(int j=0;j<physicaldim;j++)
    {
      Matrix dummy(1,1);
      dummy(1,1) = 1.0;
      SiteMatrices[0][j] = std::vector<Matrix>(1, dummy);
      Leftindex[0][j] = std::vector<int>(1,0);
      Rightindex[0][j] = std::vector<int>(1,0);
      LefttoRight[0][j] = std::vector<int>(1,0);
      RighttoLeft[0][j] = std::vector<int>(1,0);
    }


    
    std::vector<int> sites, spinsites;

    int new_site, wave_site;
    new_site = 0;
    sites.push_back(new_site);
    if (dmrginp.spinAdapted())
      spinsites.push_back(new_site);
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    StateInfo stateinfo, combinedstateinfo; 
    std::vector<Matrix> RotationMatrix;

    makeStateInfo(stateinfo, new_site);
    //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
    //LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

    for (int i=1; i<num_spatial_orbs-1; i++) {
      new_site++;
      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }

      pout <<"Sites: ";
      for (int k : sites)
        pout << k<<" ";
      pout <<endl;


      StateInfo siteState;
      StateInfo leftStateInfo=stateinfo;
      //TensorProduct provides a pointer to left and right StateInfo.
      //Load Stateinfo will overwrite the left StateInfo pointed.

      makeStateInfo(siteState, new_site);
      //TODO
      //TensorProduct of StateInfo will not removed previous data in StateInfo.
      combinedstateinfo = StateInfo(); 
      TensorProduct(leftStateInfo, siteState, combinedstateinfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

      combinedstateinfo.CollectQuanta();
      StateInfo uncollectedstateinfo = *(combinedstateinfo.unCollectedStateInfo);
      pout << uncollectedstateinfo <<endl;
      pout << combinedstateinfo<<endl;

      LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

      SpinAdapted::StateInfo::transform_state(RotationMatrix, combinedstateinfo, stateinfo);
      //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
      pout << stateinfo <<endl;



      for(int j=0;j<physicaldim;j++)
      {
        SiteMatrices[i][j].resize(leftStateInfo.quanta.size());
        Leftindex[i][j] = std::vector<int>(leftStateInfo.quanta.size(), -1);
        LefttoRight[i][j] = std::vector<int>(leftStateInfo.quanta.size(), -1);
        Rightindex[i][j] = std::vector<int>(stateinfo.quanta.size(), -1);
        RighttoLeft[i][j] = std::vector<int>(stateinfo.quanta.size(), -1);
      }



          pout <<combinedstateinfo<<endl;
          pout <<stateinfo<<endl;
      for (int totalquanta=0; totalquanta<combinedstateinfo.quanta.size();totalquanta++)
      {
        //Some quanta num discarded after renormalization.
        auto iter = std::find(stateinfo.quanta.begin(), stateinfo.quanta.end(), combinedstateinfo.quanta[totalquanta]);
        if (iter !=stateinfo.quanta.end() )
        {
          int rightq = std::distance( stateinfo.quanta.begin(), iter );
          int quantastate =0;

          for(int k=0; k < combinedstateinfo.oldToNewState[totalquanta].size(); k++)
          {
            int oldquanta = combinedstateinfo.oldToNewState[totalquanta][k];
            int leftq = uncollectedstateinfo.leftUnMapQuanta[oldquanta];
            int dotq = uncollectedstateinfo.rightUnMapQuanta[oldquanta];

            Leftindex[i][dotq][leftq] = leftq;
            Rightindex[i][dotq][rightq] = leftq;
            LefttoRight[i][dotq][leftq] = rightq;
            RighttoLeft[i][dotq][rightq] = leftq;


            //TODO
            //Discarded quanta is still in Rotation matrix.
            int leftquantastates = uncollectedstateinfo.quantaStates[oldquanta];
            int rightquantastates = stateinfo.quantaStates[rightq];


            SiteMatrices[i][dotq][leftq].ReSize(leftquantastates, rightquantastates);
            for(int l=0;l<leftquantastates;l++)
            for(int r=0;r<rightquantastates;r++)
            {
             SiteMatrices[i][dotq][leftq](l+1,r+1) = RotationMatrix[totalquanta](l+1+quantastate,r+1);
            }
            quantastate += leftquantastates;


          }
        }
      }
    }
    new_site++;

    //TODO
    //The wave function.
    std::vector<int> complementarySites;
    if (dmrginp.spinAdapted())
    {
      complementarySites.assign(1,dmrginp.last_site()-1);
    }
    else
    {
      complementarySites.push_back(dmrginp.last_site()/2-1);

    }
    //wave function and rotation matrix files are named by spin orbitals.
    //StateInfo files are named by spatial orbitals.
    getComplementarySites(spinsites,complementarySites);

    Wavefunction lastsitewave;
    StateInfo waveinfo;
    //TODO
    lastsitewave.LoadWavefunctionInfo (waveinfo, spinsites, currentstate);
    norm = DotProduct(lastsitewave,lastsitewave);

    pout <<"Last site state info " <<endl;
    pout <<stateinfo<<endl;
    pout <<waveinfo<<endl;
    pout <<*(waveinfo.leftStateInfo)<<endl;
    pout <<*(waveinfo.rightStateInfo)<<endl;
    pout <<lastsitewave<<endl;

    //TODO
    //Assuming quanta num waveinfo is not collected.
    //There is one entangle bond between the last site and other sites in each quanta num.

    StateInfo leftstateinfo = combinedstateinfo;
    combinedstateinfo = StateInfo(); 
    StateInfo siteState;
    makeStateInfo(siteState, new_site);
    TensorProduct(stateinfo, siteState, combinedstateinfo, PARTICLE_SPIN_NUMBER_CONSTRAINT);


    for(int j=0;j<physicaldim;j++)
    {
      SiteMatrices.back()[j].resize(siteState.quanta.size());
      Leftindex.back()[j] = std::vector<int>(stateinfo.quanta.size(), -1);
      LefttoRight.back()[j] = std::vector<int>(stateinfo.quanta.size(), -1);
      Rightindex.back()[j] = std::vector<int>(1, -1);
      RighttoLeft.back()[j] = std::vector<int>(1, -1);
      //Rightindex.back()[j] = std::vector<int>(combinedstateinfo.quanta.size(), -1);
      //RighttoLeft.back()[j] = std::vector<int>(combinedstateinfo.quanta.size(), -1);
    }
    //cout <<"lastsite"<< lastsitewave<<endl;
    pout <<"combinedstateinfo" <<combinedstateinfo<<endl;
    //cout <<<"wavestateinfo" <<waveinfo<<endl;


    //TODO
    //combinedstate and stateinfo need renormed quant num;
    //Wave function on the last site are three index wave function, built on unrenormed basis;
    std::vector<int> renormed_to_unrenormed;

    for(int j=0;j<RotationMatrix.size();j++)
    {
      if(RotationMatrix[j].Ncols()!=0)
        renormed_to_unrenormed.push_back(j);
    }
    for (int rightq=0; rightq<combinedstateinfo.quanta.size();rightq++)
    {
      int quantastate =0;

      int leftq = combinedstateinfo.leftUnMapQuanta[rightq];
      int dotq = combinedstateinfo.rightUnMapQuanta[rightq];

      Leftindex.back()[dotq][leftq] = leftq;
      Rightindex.back()[dotq][0] = leftq;
      LefttoRight.back()[dotq][leftq] = 0;
      RighttoLeft.back()[dotq][0] = leftq;

      int leftquantastates =  combinedstateinfo.leftStateInfo->quantaStates[leftq];
      int rightquantastates = combinedstateinfo.rightStateInfo->quantaStates[rightq];
      SiteMatrices.back()[dotq][leftq].ReSize(leftquantastates, rightquantastates);



      //int uncollectl = stateinfo.unCollectedStateInfo->leftUnMapQuanta[leftq];
      //int uncollectr = stateinfo.unCollectedStateInfo->rightUnMapQuanta[leftq];
      int unrenormedleftq =renormed_to_unrenormed[leftq]; 
	    Matrix& lM = lastsitewave.operator_element(unrenormedleftq, dotq);
	    Matrix tM = RotationMatrix[unrenormedleftq];
      if (tM.Ncols()==0) continue;
	    Matrix nM(1,1);
      nM(1,1) = 0.0;
	    MatrixMultiply(tM, 't', lM, 'n', nM, 1.0);
      SiteMatrices.back()[dotq][leftq] = nM;
      
    }

}

double Abelianmps::getcoeff(const bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters = num_spatial_orbs-1;

    // Index for spinAdapted  |0> : 0, |beta> : 1, |alpha> : 2, |alpha beta> : 3
    // Index for !spinAdapted |0> : 0, |1> : 1, 
    // leftq : quantum number of the first site with spinAdapted index
    int leftq = dmrginp.spinAdapted()? ci[0]: 2*ci[0]+ci[1];
    //TODO
    //call DGEMV rather than DGEMM
    Matrix lbasis(1,1);
    lbasis(1,1) = 1.0;
    for(int i=1; i<num_spatial_orbs;i++)
    {
      // dotq: quantum number of the ith site with spinAdapted index
      int dotq = dmrginp.spinAdapted()? ci[i]: 2*ci[2*i]+ci[2*i+1];
      int rightq = LefttoRight[i][dotq][leftq];
      if(rightq==-1 ) return 0.0;

      Matrix& rotationmatrix = SiteMatrices[i][dotq][leftq];
      
      Matrix newlbasis;
      newlbasis.ReSize(lbasis.Nrows(), rotationmatrix.Ncols());
      newlbasis= 0.0;
      MatrixMultiply(lbasis,'n', rotationmatrix,'n',newlbasis,1.0);

      leftq = rightq;
      lbasis = newlbasis;
      //pout <<lbasis<<endl;

    }
    assert(lbasis.Ncols()==1);
    assert(lbasis.Nrows()==1);
    //pout <<lbasis<<endl;
    mpscoeffT += (std::clock() - startTime)/ (double) CLOCKS_PER_SEC;
    return lbasis(1,1);
}

double Abelianmps::sampling(bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    // tmp for sampling determinant 
    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);



    ColumnVector oldwave(1);
    oldwave(1) = 1.0;
    int rightquanta = 0;
    // may be generating conditional probability for each site
    for(int i=num_spatial_orbs-1; i>0 ;i--)
    {
      std::vector<double> coeff2(physicaldim,0.0);
      std::vector<ColumnVector> newwave(physicaldim);
      for(int localquanta=0;localquanta<physicaldim;localquanta++)
      {

        int leftq = RighttoLeft[i][localquanta][rightquanta];
        if(leftq==-1) {coeff2[localquanta] = 0.0; continue;}
        Matrix& rotM = SiteMatrices[i][localquanta][leftq];
        newwave[localquanta].ReSize(rotM.Nrows());
        newwave[localquanta] = 0.0;
        MatrixMultiply(rotM,'n', oldwave,'n',newwave[localquanta], 1.0, 0.0);
        coeff2[localquanta] = dotproduct(newwave[localquanta],newwave[localquanta]);
      }

      std::discrete_distribution<int> distribution(coeff2.begin(),coeff2.end());
      int localquanta = distribution(generator);
      oldwave = newwave[localquanta];
      reverseci.push_back(localquanta%2);
      reverseci.push_back(localquanta/2);
      rightquanta = RighttoLeft[i][localquanta][rightquanta];

      if (i==1)
      {
        reverseci.push_back(rightquanta%2);
        reverseci.push_back(rightquanta/2);
      }
    }

    ci.reset();
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(oldwave.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    return oldwave(1);
}

void Abelianmps::initstate()
{
    bool direction=false;
    SweepParams sweepParams;
    dmrginp.set_algorithm_method() = ONEDOT;
    sweepParams.current_root() = currentstate;
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
}

void NonAbelianmps::initstate()
{
    bool direction=false;
    SweepParams sweepParams;
    dmrginp.set_algorithm_method() = ONEDOT;
    sweepParams.current_root() = currentstate;
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, direction, currentstate);
    Sweep::CanonicalizeWavefunction(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, !direction, currentstate);
    Sweep::InitializeStateInfo(sweepParams, direction, currentstate);
}

void Abelianmps::build(int state)
{
  currentstate = state;
  if (mpigetrank() == 0) {
    std::string file;
    file = str(boost::format("%s%s%d%s") % dmrginp.load_prefix() % "/simplemps-" % state % ".tmp" );
    if (!dmrginp.stochasticpt_restart())
    {
      initstate();
      SweepParams sweepParams;
      sweepParams.set_sweep_iter() = 0;
      sweepParams.current_root() = state;
      convertblockfiles();
      std::ofstream ofs(file.c_str(), std::ios::binary);
      boost::archive::binary_oarchive save_state(ofs);
      save_state << *this;
    }
    else{
      std::ifstream ifs(file.c_str(), std::ios::binary);
      boost::archive::binary_iarchive save_state(ifs);
      save_state >> *this;

    }
  }

#ifndef SERIAL
      boost::mpi::communicator world;
      broadcast(world, *this, 0);
#endif
}

void NonAbelianmps::build(int state)
{
  currentstate = state;
  if (mpigetrank() == 0) {
    std::string file;
    file = str(boost::format("%s%s%d%s") % dmrginp.load_prefix() % "/simplemps-" % state % ".tmp" );
    if (!dmrginp.stochasticpt_restart())
    {
      initstate();
      SweepParams sweepParams;
      sweepParams.set_sweep_iter() = 0;
      sweepParams.current_root() = state;
      convertblockfiles();
      std::ofstream ofs(file.c_str(), std::ios::binary);
      boost::archive::binary_oarchive save_state(ofs);
      save_state << *this;
    }
    else{
      std::ifstream ifs(file.c_str(), std::ios::binary);
      boost::archive::binary_iarchive save_state(ifs);
      save_state >> *this;
    }
  }

  //TODO
  //boost serialization does not work for map (or unordered_map), even though it is claimed working.
  //Change map to vector of keys and vector of values.
  std::vector<std::vector<std::vector<std::pair<int,int> > > > keys;
  std::vector<std::vector<std::vector<Matrix> > > values;
#ifndef SERIAL
  boost::mpi::communicator world;

  broadcast(world, num_spatial_orbs, 0);
  broadcast(world, spin_vec, 0);
  broadcast(world, norm, 0);
  broadcast(world, physicaldim, 0);

  world.barrier();


  if(mpigetrank()!=0)
  {
    Leftindex.resize(num_spatial_orbs);
    Rightindex.resize(num_spatial_orbs);
    LefttoRight.resize(num_spatial_orbs);
    RighttoLeft.resize(num_spatial_orbs);
  }
  world.barrier();

  for(int i=0;i<num_spatial_orbs;i++)
  {
      broadcast(world, LefttoRight[i], 0);
      broadcast(world, RighttoLeft[i], 0);
      broadcast(world, Leftindex[i], 0);
      broadcast(world, Rightindex[i], 0);
  }
  world.barrier();


// broadcast in 1st loop 
  keys.resize(num_spatial_orbs);
  values.resize(num_spatial_orbs);
  if(mpigetrank()==0)
  {
     for(int i=0;i<SiteMatrices.size();i++)
     {
        keys[i].resize(SiteMatrices[i].size());
        values[i].resize(SiteMatrices[i].size());
        for(int j=0;j<SiteMatrices[i].size();j++)
        {
          //keys[i][j].resize(SiteMatrices[i][j].size());
          //values[i][j].resize(SiteMatrices[i][j].size());
            for(auto iter=SiteMatrices[i][j].begin();iter!=SiteMatrices[i][j].end();iter++)
            {
              keys[i][j].push_back(iter->first);
              values[i][j].push_back(iter->second);
            }
        }
     }
  }
  for(int i=0;i<num_spatial_orbs;i++)
  {
     broadcast(world, keys[i], 0);
     broadcast(world, values[i], 0);
  }

// broadcast in 2nd loop 
//   keys.resize(num_spatial_orbs);
//   values.resize(num_spatial_orbs);
// 
//   for(int i=0;i<num_spatial_orbs;i++)
//   {
//      keys[i].resize(physicaldim);
//      values[i].resize(physicaldim);
//      for(int j=0;j<SiteMatrices[i].size();j++)
//      {
//          if(mpigetrank()==0)
//          {
//             //keys[i][j].resize(SiteMatrices[i][j].size());
//             //values[i][j].resize(SiteMatrices[i][j].size());
//               for(auto iter=SiteMatrices[i][j].begin();iter!=SiteMatrices[i][j].end();iter++)
//               {
//                 keys[i][j].push_back(iter->first);
//                 values[i][j].push_back(iter->second);
//               }
//          }
//          broadcast(world, keys[i][j], 0);
//          broadcast(world, values[i][j], 0);
//      }
//   }

  if(mpigetrank()!=0)
  {
    SiteMatrices.resize(keys.size());
    for(int i=0;i<keys.size();i++)
    {
      SiteMatrices[i].resize(keys[i].size());
      for(int j=0;j<keys[i].size();j++)
      {
        //keys[i][j].resize(SiteMatrices[i][j].size());
        //values[i][j].resize(SiteMatrices[i][j].size());
          for(int k=0;k<keys[i][j].size();k++)
          {
            SiteMatrices[i][j][keys[i][j][k]] = values[i][j][k];
          }
      }
    }
  }
  world.barrier();

#endif
}

void NonAbelianmps::convertblockfiles()
{
    physicaldim =dmrginp.spinAdapted()?3:4; 
    num_spatial_orbs = dmrginp.spinAdapted()?dmrginp.last_site():dmrginp.last_site()/2;
    SiteMatrices.resize(num_spatial_orbs);
    Leftindex.resize(num_spatial_orbs);
    Rightindex.resize(num_spatial_orbs);
    LefttoRight.resize(num_spatial_orbs);
    RighttoLeft.resize(num_spatial_orbs);
    spin_vec.resize(num_spatial_orbs+1);
    spin_vec.back() = {0};//Singlet embedding is used.
    for(int i=0;i<num_spatial_orbs;i++)
    {
      SiteMatrices[i].resize(physicaldim);
      Leftindex[i].resize(physicaldim);
      Rightindex[i].resize(physicaldim);
      LefttoRight[i].resize(physicaldim);
      RighttoLeft[i].resize(physicaldim);
    }

    
    std::vector<int> sites, spinsites;

    int new_site, wave_site;
    new_site = 0;
    sites.push_back(new_site);
    if (dmrginp.spinAdapted())
      spinsites.push_back(new_site);
    else {
      spinsites.push_back(2*new_site);
      spinsites.push_back(2*new_site+1);
    }
    StateInfo stateinfo, combinedstateinfo; 
    std::vector<Matrix> RotationMatrix;

    makeStateInfo(stateinfo, new_site);
    //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
    //LoadRotationMatrix (spinsites, RotationMatrix, currentstate);
    
    Matrix one(1,1);
    spin_vec[0] = {dmrginp.molecule_quantum().get_s().getirrep() };
    spin_vec[1].resize(stateinfo.quanta.size());
    for(int j=0;j<stateinfo.quanta.size();j++)
    {
      spin_vec[1][j] = stateinfo.quanta[j].get_s().getirrep();
    }
    if (dmrginp.add_noninteracting_orbs() && dmrginp.molecule_quantum().get_s().getirrep() != 0 && dmrginp.spinAdapted())
    {
        for(int j=0;j<physicaldim;j++)
        {
          Leftindex[0][j] = std::vector<std::vector<int> >(1);
          Rightindex[0][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
          LefttoRight[0][j] = std::vector<std::vector<int> >(1);
          RighttoLeft[0][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
        }

        for (int totalquanta=0; totalquanta<stateinfo.quanta.size();totalquanta++)
        {
            int dotq = stateinfo.leftUnMapQuanta[totalquanta];
            int leftq = stateinfo.rightUnMapQuanta[totalquanta];//It is the non-interacting orbs for singlet embeding.
            Matrix dummy(1,1);
            dummy(1,1) = 1.0;
            SiteMatrices[0][dotq][make_pair(leftq,totalquanta)] = dummy;
            Leftindex[0][dotq][leftq].push_back(totalquanta);
            LefttoRight[0][dotq][leftq].push_back(totalquanta);
            RighttoLeft[0][dotq][totalquanta].push_back(leftq);
            Rightindex[0][dotq][totalquanta].push_back(leftq);

        }
    }
    else{
        Matrix dummy(1,1);
        dummy(1,1) = 1.0;
        for(int j=0;j<physicaldim;j++)
        {
          SiteMatrices[0][j][make_pair(0,j)] = dummy;
          Leftindex[0][j] = std::vector<std::vector<int> >(1);
          Rightindex[0][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
          LefttoRight[0][j] = std::vector<std::vector<int> >(1);
          RighttoLeft[0][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
          LefttoRight[0][j][0] = {j};
          RighttoLeft[0][j][j] = {0};
        }

    }


    for (int i=1; i<num_spatial_orbs-1; i++) {
      new_site++;
      sites.push_back(new_site);
      if (dmrginp.spinAdapted())
        spinsites.push_back(new_site);
      else {
        spinsites.push_back(2*new_site);
        spinsites.push_back(2*new_site+1);
      }

      pout <<"Sites: ";
      for (int k : sites)
        pout << k<<" ";
      pout <<endl;


      StateInfo siteState;
      StateInfo leftStateInfo=stateinfo;
      //TensorProduct provides a pointer to left and right StateInfo.
      //Load Stateinfo will overwrite the left StateInfo pointed.

      makeStateInfo(siteState, new_site);
      //TODO
      //TensorProduct of StateInfo will not removed previous data in StateInfo.
      combinedstateinfo = StateInfo(); 
      TensorProduct(leftStateInfo, siteState, combinedstateinfo, NO_PARTICLE_SPIN_NUMBER_CONSTRAINT);

      combinedstateinfo.CollectQuanta();
      StateInfo uncollectedstateinfo = *(combinedstateinfo.unCollectedStateInfo);
      pout << uncollectedstateinfo <<endl;
      pout << combinedstateinfo<<endl;

      LoadRotationMatrix (spinsites, RotationMatrix, currentstate);

      SpinAdapted::StateInfo::transform_state(RotationMatrix, combinedstateinfo, stateinfo);
      //StateInfo::restore(forward, sites, stateinfo, sweepParams.current_root());
      pout << stateinfo <<endl;



      for(int j=0;j<physicaldim;j++)
      {
        //SiteMatrices[i][j].resize(leftStateInfo.quanta.size());
        Leftindex[i][j] = std::vector<std::vector<int> >(leftStateInfo.quanta.size());
        LefttoRight[i][j] = std::vector<std::vector<int> >(leftStateInfo.quanta.size());
        Rightindex[i][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
        RighttoLeft[i][j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
      }
      spin_vec[i+1].resize(stateinfo.quanta.size());
      for(int j=0;j<stateinfo.quanta.size();j++)
      {
        spin_vec[i+1][j] = stateinfo.quanta[j].get_s().getirrep();
      }



          pout <<combinedstateinfo<<endl;
          pout <<stateinfo<<endl;
      for (int totalquanta=0; totalquanta<combinedstateinfo.quanta.size();totalquanta++)
      {
        //Some quanta num discarded after renormalization.
        auto iter = std::find(stateinfo.quanta.begin(), stateinfo.quanta.end(), combinedstateinfo.quanta[totalquanta]);
        if (iter !=stateinfo.quanta.end() )
        {
          int rightq = std::distance( stateinfo.quanta.begin(), iter );
          int quantastate =0;

          for(int k=0; k < combinedstateinfo.oldToNewState[totalquanta].size(); k++)
          {
            int oldquanta = combinedstateinfo.oldToNewState[totalquanta][k];
            int leftq = uncollectedstateinfo.leftUnMapQuanta[oldquanta];
            int dotq = uncollectedstateinfo.rightUnMapQuanta[oldquanta];

            //TODO
            //Leftindex[i][dotq][leftq].push_back(oldquanta);
            Rightindex[i][dotq][rightq].push_back(leftq);
            LefttoRight[i][dotq][leftq].push_back(rightq);
            RighttoLeft[i][dotq][rightq].push_back(leftq);


            //TODO
            //Discarded quanta is still in Rotation matrix.
            int leftquantastates = uncollectedstateinfo.quantaStates[oldquanta];
            int rightquantastates = stateinfo.quantaStates[rightq];


            Matrix RotM;
            RotM.ReSize(leftquantastates, rightquantastates);
            for(int l=0;l<leftquantastates;l++)
            for(int r=0;r<rightquantastates;r++)
            {
             RotM(l+1,r+1) = RotationMatrix[totalquanta](l+1+quantastate,r+1);
            }
            SiteMatrices[i][dotq][make_pair(leftq, rightq)] = RotM;
            quantastate += leftquantastates;


          }
        }
      }
    }
    new_site++;

    //TODO
    //The wave function.
    std::vector<int> complementarySites;
    if (dmrginp.spinAdapted())
    {
      complementarySites.assign(1,dmrginp.last_site()-1);
    }
    else
    {
      complementarySites.push_back(dmrginp.last_site()/2-1);

    }
    //wave function and rotation matrix files are named by spin orbitals.
    //StateInfo files are named by spatial orbitals.
    getComplementarySites(spinsites,complementarySites);

    Wavefunction lastsitewave;
    StateInfo waveinfo;
    //TODO
    lastsitewave.LoadWavefunctionInfo (waveinfo, spinsites, currentstate);
    norm = DotProduct(lastsitewave,lastsitewave);

    pout <<"Last site state info " <<endl;
    pout <<stateinfo<<endl;
    pout <<waveinfo<<endl;
    pout <<*(waveinfo.leftStateInfo)<<endl;
    pout <<*(waveinfo.rightStateInfo)<<endl;
    pout <<lastsitewave<<endl;

    //TODO
    //Assuming quanta num waveinfo is not collected.
    //There is one entangle bond between the last site and other sites in each quanta num.

    StateInfo leftstateinfo = combinedstateinfo;
    combinedstateinfo = StateInfo(); 
    StateInfo siteState;
    makeStateInfo(siteState, new_site);
    TensorProduct(stateinfo, siteState, combinedstateinfo, PARTICLE_SPIN_NUMBER_CONSTRAINT);


    for(int j=0;j<physicaldim;j++)
    {
      //SiteMatrices.back()[j].resize(siteState.quanta.size());
      //Leftindex.back()[j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
      LefttoRight.back()[j] = std::vector<std::vector<int> >(stateinfo.quanta.size());
      //Rightindex.back()[j] = std::vector<std::vector<int> >(1);
      RighttoLeft.back()[j] = std::vector<std::vector<int> >(1);
      //Rightindex.back()[j] = std::vector<int>(combinedstateinfo.quanta.size(), -1);
      //RighttoLeft.back()[j] = std::vector<int>(combinedstateinfo.quanta.size(), -1);
    }
    //cout <<"lastsite"<< lastsitewave<<endl;
    pout <<"combinedstateinfo" <<combinedstateinfo<<endl;
    //cout <<<"wavestateinfo" <<waveinfo<<endl;


    //TODO
    //combinedstate and stateinfo need renormed quant num;
    //Wave function on the last site are three index wave function, built on unrenormed basis;
    std::vector<int> renormed_to_unrenormed;

    for(int j=0;j<RotationMatrix.size();j++)
    {
      if(RotationMatrix[j].Ncols()!=0)
        renormed_to_unrenormed.push_back(j);
    }
    for (int rightq=0; rightq<combinedstateinfo.quanta.size();rightq++)
    {
      int quantastate =0;

      int leftq = combinedstateinfo.leftUnMapQuanta[rightq];
      int dotq = combinedstateinfo.rightUnMapQuanta[rightq];

      //Leftindex.back()[dotq][leftq] = leftq;
      //Rightindex.back()[dotq][0] = leftq;
      LefttoRight.back()[dotq][leftq].push_back(0);
      RighttoLeft.back()[dotq][0].push_back(leftq);

      int leftquantastates =  combinedstateinfo.leftStateInfo->quantaStates[leftq];
      int rightquantastates = combinedstateinfo.rightStateInfo->quantaStates[rightq];
      //SiteMatrices.back()[dotq][leftq].ReSize(leftquantastates, rightquantastates);



      //int uncollectl = stateinfo.unCollectedStateInfo->leftUnMapQuanta[leftq];
      //int uncollectr = stateinfo.unCollectedStateInfo->rightUnMapQuanta[leftq];
      int unrenormedleftq =renormed_to_unrenormed[leftq]; 
	    Matrix& lM = lastsitewave.operator_element(unrenormedleftq, dotq);
	    Matrix tM = RotationMatrix[unrenormedleftq];
      if (tM.Ncols()==0) continue;
	    Matrix nM(1,1);
      nM(1,1) = 0.0;
	    MatrixMultiply(tM, 't', lM, 'n', nM, 1.0);
      SiteMatrices.back()[dotq][make_pair(leftq,0)] = nM;
      
    }

}


double NonAbelianmps::getcoeff(const bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    int niters = num_spatial_orbs-1;

    std::vector<int> leftq_vec = {0};
    //TODO
    //call DGEMV rather than DGEMM
    Matrix one(1,1);
    one(1,1) = 1.0;
    std::unordered_map<int, Matrix> lbasis;
    lbasis[0] = one;
    //int Ms = -ci.Sz();
    //TODO
    //In singlet embedding, there is an ensemble of 2*s+1 state with different Sz. The norm of each state is 1/sqrt(2*s+1);
    //It is possible to only use the state with Sz=S, however, the it is hard to only sample with certein Sz, without discarding determinant at last step.
    //If some theories required a pure state with Sz=S. Add the following line.
    lbasis[0]*=sqrt(dmrginp.molecule_quantum().get_s().getirrep()+1);
    //int Ms = -dmrginp.molecule_quantum().get_s().getirrep();
    int Ms = -ci.Sz();

    for(int i=0; i<num_spatial_orbs;i++)
    {
      std::unordered_map<int, Matrix>  new_lbasis;
      int dotq = ci[2*i]+ci[2*i+1];
      int dotspin = dotq%2;
      int dotMs = ci[2*i]-ci[2*i+1];
      for(auto iter=lbasis.begin();iter!=lbasis.end();iter++ )
      {
        int leftq = iter->first;
        for(int rightq: LefttoRight[i][dotq][leftq])
        {
          Matrix& rotationmatrix = SiteMatrices[i][dotq][make_pair(leftq,rightq)];
          double cgfactor = cg(spin_vec[i][leftq], dotspin, spin_vec[i+1][rightq], Ms, dotMs, Ms+dotMs);
          //TODO
          //In block, the noninteracting orbital is on the right of first physical spatial orbital.
          if(i==0)
          {
            cgfactor = cg(dotspin, spin_vec[i][leftq], spin_vec[i+1][rightq], dotMs, Ms, Ms+dotMs);
            if(dotspin%2 && spin_vec[i][leftq]%2)
              cgfactor *=-1;
          }
          //if(i==0) cgfactor =1.0;
          //cout <<"cgfactor"<< cgfactor<<endl;
          //cout <<"S: " << spin_vec[i][leftq]<<dotspin<<spin_vec[i+1][rightq]<<endl;
          //cout <<"Sz: " << Ms<<dotMs<<Ms+dotMs<<endl;
          if(new_lbasis.find(rightq)==new_lbasis.end())
          {
            new_lbasis[rightq] = Matrix();
            new_lbasis[rightq].ReSize(1, rotationmatrix.Ncols());
            new_lbasis[rightq] = 0.0;
          }
          MatrixMultiply(lbasis[leftq],'n', rotationmatrix,'n',new_lbasis[rightq], cgfactor);
        }
      }
      lbasis = new_lbasis;
      //cout <<"step "<<i<<endl;
      //cout <<"lbasis"<<endl;
      //for(auto iter=lbasis.begin();iter!=lbasis.end();iter++)
      //  cout <<iter->first<<" : "<<endl<<iter->second;
      if(lbasis.size()==0) return 0.0;
      Ms =Ms+ci[2*i]-ci[2*i+1];

    }
    assert(lbasis.size()==1);
    assert(lbasis[0].Nrows()==1);
    assert(lbasis[0].Ncols()==1);
    //pout <<lbasis<<endl;
    mpscoeffT += (std::clock() - startTime)/ (double) CLOCKS_PER_SEC;
    return lbasis[0](1,1);
}

double NonAbelianmps::sampling(bitstring& ci)
{
    std::clock_t startTime = std::clock();
    int statenum=0;

    std::vector<int> reverseci;

    std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> f(0.0,1.0);

    const int physicaldim =4;
    ColumnVector oldwave(1);
    oldwave(1) = 1.0;
    //TODO
    //In singlet embedding, there is an ensemble of 2*s+1 state with different Sz. The norm of each state is 1/sqrt(2*s+1);
    //It is possible to only use the state with Sz=S, however, the it is hard to only sample with certein Sz, without discarding determinant at last step.
    //If some theories required a pure state with Sz=S. Add the following line.
    oldwave(1) = sqrt(dmrginp.molecule_quantum().get_s().getirrep()+1);

    // rwave: 
    //
    //   |
    // --*-- 0 (rwave[0]: target quanta)
    int rightquanta = 0;
    std::unordered_map<int,ColumnVector> rwave;
    rwave[0] = oldwave;
    int Ms = 0;//This is 2*Sz of right quanta num.
    double prob = 1.0;
    //prob*=dmrginp.molecule_quantum().get_s().getirrep()+1;

    for(int i=num_spatial_orbs-1; i>=0 ;i--)
    {
      std::vector<double> coeff2(physicaldim,0.0);
      std::vector<std::unordered_map<int,ColumnVector>> newwave(physicaldim);
      for(int localquanta=0;localquanta<physicaldim;localquanta++)
      {
        // dotq   : # of electron of dot
        // dotspin: <S2> of dot
        // dotMs  : 2*Sz of dot
        //   |
        // --*-- 0 
        int dotq = (localquanta+1)/2;
        int dotspin = dotq%2;
        int dotMs = localquanta/2-localquanta%2;

        for(auto iter=rwave.begin();iter!=rwave.end();iter++)
        {
          int rightq = iter->first;
          //              dot
          //               |
          // RighttoLeft --*-- right
          //               i
          for(int leftq :RighttoLeft[i][dotq][rightq])
          {
            // SiteMatrice for each quantum #: ==*--
            Matrix& rotM = SiteMatrices[i][dotq][make_pair(leftq,rightq)];
            double cgfactor = cg(spin_vec[i][leftq], dotspin, spin_vec[i+1][rightq], -Ms-dotMs, dotMs, -Ms);

           
            if(newwave[localquanta].find(leftq)==newwave[localquanta].end())
            {
              newwave[localquanta][leftq] = ColumnVector();
              newwave[localquanta][leftq].ReSize(rotM.Nrows());
              newwave[localquanta][leftq] = 0.0;
            }
            MatrixMultiply(rotM,'n', rwave[rightq],'n',newwave[localquanta][leftq], cgfactor);

          }
        }
        for(auto iter=newwave[localquanta].begin();iter!=newwave[localquanta].end();iter++)
        {
          coeff2[localquanta] += dotproduct(iter->second,iter->second);
        }
      }

      std::discrete_distribution<int> distribution(coeff2.begin(),coeff2.end());
      //cout <<"Sum: "<< std::accumulate(coeff2.begin(),coeff2.end(),0.0)<<endl;
      //for(auto x: coeff2)
      //  cout <<x<<" , ";
      //cout <<endl;
      int localquanta = distribution(generator);
      prob *=distribution.probabilities()[localquanta];
      rwave = newwave[localquanta];
      if(localquanta==1) Ms-=1;
      else if(localquanta==2) Ms+=1;
      reverseci.push_back(localquanta%2);
      reverseci.push_back(localquanta/2);

      //if (i==1)
      //{
      //  assert(rwave.size()==1);
      //  int rightq = rwave.begin()->first;
      //  if(rightq ==0)
      //  {
      //    reverseci.push_back(0);
      //    reverseci.push_back(0);
      //  }
      //  else if(rightq==2)
      //  {
      //    reverseci.push_back(1);
      //    reverseci.push_back(1);
      //  }
      //  else{
      //    if(dmrginp.molecule_quantum().get_s().getirrep()>Ms)
      //    {
      //      reverseci.push_back(0);
      //      reverseci.push_back(1);
      //    }
      //    else{
      //      reverseci.push_back(1);
      //      reverseci.push_back(0);
      //    }

      //  }
      //}
    }
      //cout <<"Ms: "<<Ms<<endl;
      //for(auto i: reverseci)
      //  cout <<i<<" ";
      //cout <<endl;

    ci.reset();
    for(int i=0;i<reverseci.size();i++)
      ci[i] = reverseci[reverseci.size()-1-i];
    assert(rwave.size()==1);
    assert(rwave.begin()->second.Nrows()==1);
    mpssampleT += ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    //cout <<"Sample probability: " << prob<<endl;
    //if(ci.Sz()!=dmrginp.molecule_quantum().get_s().getirrep()) return 0.0;
    return rwave.begin()->second(1);
}

void NonAbelianmps::lw_dot_mat(const int d, const int j, SQ_RV &lw, SQ_RV &nw)
{
    //   RowVector(lw)   SiteMatrices                     RowVector(nw)
    //                        dotq
    //                         |                        =  
    // leftq --*-- leftq .   --*-- rightq = LefttoRight      --*--
    //         d              d+1                             d+1

    int dotq = (j+1)/2; 
    int dotspin = dotq%2;
    int dotMs = j/2-j%2;
    // dotq   : # of electron of dot
    // dotq    = 0,  1, 1, 2 (for j = 0, 1, 2, 3)  
    // dotspin: <S2> of dot
    // dotspin = 0,  1, 1, 0 (for j = 0, 1, 2, 3) 
    // dotMs  : 2*Sz of dot
    // dotMs   = 0, -1, 1, 0 (for j = 0, 1, 2, 3)
    // |0>, |beta>, |alpha>, |alpha beta>  (for j = 0, 1, 2, 3)


    for(auto &m : lw){
        int Ms    = std::get<3>(m.first);    // Ms for left 
        int leftq = std::get<2>(m.first);    // rightq in lw is leftq for nw

//        //lsh test ======================================================== 
//        std::cout << "d, leftq, dotq, rightq size = " << d << " "
//                  << leftq << " " << dotq << " " << LefttoRight[d][dotq][leftq].size()
//                  << " rightq[0] = " << LefttoRight[d][dotq][leftq][0] << std::endl;
//    
//        std::cout << "lw leftq " << leftq << " RowVector ="; 
//        for (int j=0; j < m.second.Ncols(); j++){
//            std::cout << j << " " << m.second.element(j); 
//        }
//        std::cout << std::endl;
//        //lsh test ======================================================== 

        for(int rightq : LefttoRight[d][dotq][leftq])
        {
          Matrix& rotM = SiteMatrices[d][dotq][std::make_pair(leftq,rightq)];
          double cgfactor;
          //cgfactor = cg(spin_vec[d+1][rightq], dotspin, spin_vec[d][leftq], Ms+dotMs, dotMs, Ms);
          cgfactor = cg(spin_vec[d][leftq], dotspin, spin_vec[d+1][rightq], Ms, dotMs, Ms+dotMs);

//          //lsh test
//          double cgt1 = cg(0.0, 1.0, 1.0, 0.0, -1.0, 1.0);  //  T = A + B
//          std::cout << " cg S = s1 + s2 " << cgt1 << std::endl;
//          double cgt2 = cg(1.0, 1.0, 0.0, -1.0, 1.0, 0.0);  //  A + B = T 
//          std::cout << " cg s1 + s2 = S " << cgt2 << std::endl; // this is correct 
//  
//          //lsh test
//          std::cout << " rotM" << std::endl ;
//          for (int i=0; i < rotM.Nrows(); i++){
//              std::cout << i << " row " ;
//              for (int j=0; j < rotM.Ncols(); j++){
//                  std::cout << j << " " << rotM.element(i, j);
//              }
//              std::cout << std::endl;
//          }
//          std::cout << "S left, dot, right Ms left, dot, right = " <<  spin_vec[d][leftq]<< " " << dotspin<< " " <<  spin_vec[d+1][rightq] << " " <<  Ms << " " <<  dotMs<< " " <<  Ms+dotMs << std::endl;
//          std::cout << " cgfactor = " << cgfactor << std::endl;

          SQ SQr = std::make_tuple(d, leftq, rightq, Ms+dotMs); 
          
          if(nw.find(SQr)==nw.end())
          {
            nw[SQr] = RowVector();
            nw[SQr].ReSize(rotM.Ncols());
            nw[SQr] = 0.0;
          }
          MatrixMultiply(m.second,'n',rotM,'n',nw[SQr], cgfactor);
    
//          //lsh test
//          std::cout << "nw leftq, rightq " << leftq << ", " << rightq << " RowVector= "; 
//          for (int j=0; j < nw[SQr].Ncols(); j++){
//              std::cout << j << " " << nw[SQr].element(j); 
//          }
//          std::cout << std::endl;
        }
     }

}

void NonAbelianmps::lw_dot_mat_new(const int d, const int j, partial_wv &lw, partial_wv &nw)
{
    //   RowVector(lw)   SiteMatrices                     RowVector(nw)
    //                        dotq
    //                         |                        =  
    // leftq --*-- leftq .   --*-- rightq = LefttoRight      --*--
    //         d              d+1                             d+1

    int dotq = (j+1)/2; 
    int dotspin = dotq%2;
    int dotMs = j/2-j%2;
    // dotq   : # of electron of dot
    // dotq    = 0,  1, 1, 2 (for j = 0, 1, 2, 3)  
    // dotspin: <S2> of dot
    // dotspin = 0,  1, 1, 0 (for j = 0, 1, 2, 3) 
    // dotMs  : 2*Sz of dot
    // dotMs   = 0, -1, 1, 0 (for j = 0, 1, 2, 3)
    // |0>, |beta>, |alpha>, |alpha beta>  (for j = 0, 1, 2, 3)

    for(auto &m : lw){
        int Ms    = std::get<3>(m.first);    // Ms for left 
        int leftq = std::get<2>(m.first);    // rightq in lw is leftq for nw

//        //lsh test ======================================================== 
//        std::cout << "d      = " << d << std::endl;
//        std::cout << "Ms     = " << Ms << std::endl;
//        std::cout << "leftq  = " << leftq << std::endl;
//        std::cout << "dotq   = " << dotq << std::endl;
//        std::cout << "rightq size = " << LefttoRight[d][dotq][leftq].size() << std::endl;
//        //std::cout << "rightq[0]   = " << LefttoRight[d][dotq][leftq][0] << std::endl;
//        std::cout << "RowVector   = ";
//        for (int j=0; j < m.second.Ncols(); j++){
//            std::cout << j << " " << m.second.element(j); 
//        }
//        std::cout << std::endl;
//        //lsh test ======================================================== 

        for(int rightq : LefttoRight[d][dotq][leftq])
        {
          Matrix& rotM = SiteMatrices[d][dotq][std::make_pair(leftq,rightq)];
          double cgfactor;
          //cgfactor = cg(spin_vec[d+1][rightq], dotspin, spin_vec[d][leftq], Ms+dotMs, dotMs, Ms);
          cgfactor = cg(spin_vec[d][leftq], dotspin, spin_vec[d+1][rightq], Ms, dotMs, Ms+dotMs);

//          //lsh test
//          double cgt1 = cg(0.0, 1.0, 1.0, 0.0, -1.0, 1.0);  //  T = A + B
//          std::cout << " cg S = s1 + s2 " << cgt1 << std::endl;
//          double cgt2 = cg(1.0, 1.0, 0.0, -1.0, 1.0, 0.0);  //  A + B = T 
//          std::cout << " cg s1 + s2 = S " << cgt2 << std::endl; // this is correct 
//  
//          //lsh test
//          std::cout << " rotM" << std::endl ;
//          for (int i=0; i < rotM.Nrows(); i++){
//              std::cout << i << " row " ;
//              for (int j=0; j < rotM.Ncols(); j++){
//                  std::cout << j << " " << rotM.element(i, j);
//              }
//              std::cout << std::endl;
//          }
//          std::cout << "S left, dot, right Ms left, dot, right = " <<  spin_vec[d][leftq]<< " " << dotspin<< " " <<  spin_vec[d+1][rightq] << " " <<  Ms << " " <<  dotMs<< " " <<  Ms+dotMs << std::endl;
//          std::cout << " cgfactor = " << cgfactor << std::endl;

          SQ SQr = std::make_tuple(d, leftq, rightq, Ms+dotMs); 

          if(nw.find(SQr)==nw.end())
          {
            nw[SQr] = RowVector();
            nw[SQr].ReSize(rotM.Ncols());
            nw[SQr] = 0.0;
          }
          MatrixMultiply(m.second,'n',rotM,'n',nw[SQr], cgfactor);
    
//          //lsh test
//          std::cout << "nw leftq, rightq " << leftq << ", " << rightq << " RowVector= "; 
//          for (int j=0; j < nw[SQr].Ncols(); j++){
//              std::cout << j << " " << nw[SQr].element(j); 
//          }
//          std::cout << std::endl;
        }
     }

}
