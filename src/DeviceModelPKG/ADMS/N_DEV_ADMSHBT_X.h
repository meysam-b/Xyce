
// -*-c++-*-
//-------------------------------------------------------------------------
//   Copyright 2002-2025 National Technology & Engineering Solutions of
//   Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
//   NTESS, the U.S. Government retains certain rights in this software.
//
//   This file is part of the Xyce(TM) Parallel Electrical Simulator.
//
//   Xyce(TM) is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   Xyce(TM) is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Xyce(TM).
//   If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//
// Purpose        :
//
// Special Notes  : Generated from verilog file fbhhbt-2.1_nonoise_limited_inductors_typed.va with ADMS
//                  interface for Xyce 7.9.0
//                  DO NOT EDIT THIS FILE DIRECTLY!  It may be overwritten!
//
// Creator        : admsXml-2.3.7
//
// Creation Date  : Thu, 23 May 2024 20:30:05
//
//-----------------------------------------------------------------------------
#ifndef Xyce_N_DEV_ADMSHBT_X_h
#define Xyce_N_DEV_ADMSHBT_X_h


#include <N_DEV_Configuration.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_BJT.h>


// Xyce_config.h contains a VERSION macro from autoconf, and some
// Verilog-A models like to define a variable of that name.  This can be
// a serious problem, and we don't need the VERSION macro.  Get rid of it.
// This must happen *after* all the includes of Xyce headers, each of which
// includes Xyce_config.h.  The implementation file must do this all again
// because it includes more Xyce headers *after* including this one.
#ifdef VERSION
#undef VERSION
#endif

namespace Xyce {
namespace Device {
namespace ADMSHBT_X {

class Model;
class Instance;
class InstanceSensitivity;

#ifdef Xyce_ADMS_SENSITIVITIES
//-----------------------------------------------------------------------------
// Class         : InstanceSensitivity
//
// Purpose       : This class is a functor for sensitivity
//
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
class InstanceSensitivity : public baseSensitivity
{
public:
InstanceSensitivity() :
baseSensitivity() {};

virtual ~InstanceSensitivity() {};

virtual void operator()(
const ParameterBase &entity,
const std::string &param,
std::vector<double> & dfdp,
std::vector<double> & dqdp,
std::vector<double> & dbdp,
std::vector<int> & Findices,
std::vector<int> & Qindices,
std::vector<int> & Bindices
) const ;
};

static InstanceSensitivity instSens;


//-----------------------------------------------------------------------------
// Class         : ModelSensitivity
//
// Purpose       : This class is a functor for sensitivity
//
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
class ModelSensitivity : public baseSensitivity
{
public:
ModelSensitivity() :
baseSensitivity() {};

virtual ~ModelSensitivity() {};

virtual void operator()(
const ParameterBase &entity,
const std::string &param,
std::vector<double> & dfdp,
std::vector<double> & dqdp,
std::vector<double> & dbdp,
std::vector<int> & Findices,
std::vector<int> & Qindices,
std::vector<int> & Bindices
) const ;
};

static ModelSensitivity modSens;
#endif // Xyce_ADMS_SENSITIVITIES

// general purpose free functions
// thermal voltage at kelvin temperature temp)
static inline double adms_vt(const double temp) {return(CONSTKoverQ*temp);};


#ifdef Xyce_ADMS_SENSITIVITIES
//-----------------------------------------------------------------------------
// "structs" to hold instance and model param/variable copies
//-----------------------------------------------------------------------------
class instanceSensStruct
{
public:
// instance parameters
// reals
double instancePar_Temp;
double d_instancePar_Temp_dX;
bool instancePar_given_Temp;
double instancePar_L;
double d_instancePar_L_dX;
bool instancePar_given_L;
double instancePar_W;
double d_instancePar_W_dX;
bool instancePar_given_W;
// non-reals(including hidden)
int instancePar_N;
bool instancePar_given_N;
};

class modelSensStruct
{
public:
// model parameters
// reals
double modelPar_Temp;
double d_modelPar_Temp_dX;
bool modelPar_given_Temp;
double modelPar_Rth;
double d_modelPar_Rth_dX;
bool modelPar_given_Rth;
double modelPar_Cth;
double d_modelPar_Cth_dX;
bool modelPar_given_Cth;
double modelPar_L;
double d_modelPar_L_dX;
bool modelPar_given_L;
double modelPar_W;
double d_modelPar_W_dX;
bool modelPar_given_W;
double modelPar_Jsf;
double d_modelPar_Jsf_dX;
bool modelPar_given_Jsf;
double modelPar_nf;
double d_modelPar_nf_dX;
bool modelPar_given_nf;
double modelPar_Vg;
double d_modelPar_Vg_dX;
bool modelPar_given_Vg;
double modelPar_Jse;
double d_modelPar_Jse_dX;
bool modelPar_given_Jse;
double modelPar_ne;
double d_modelPar_ne_dX;
bool modelPar_given_ne;
double modelPar_Rbxx;
double d_modelPar_Rbxx_dX;
bool modelPar_given_Rbxx;
double modelPar_Vgb;
double d_modelPar_Vgb_dX;
bool modelPar_given_Vgb;
double modelPar_Jsee;
double d_modelPar_Jsee_dX;
bool modelPar_given_Jsee;
double modelPar_nee;
double d_modelPar_nee_dX;
bool modelPar_given_nee;
double modelPar_Rbbxx;
double d_modelPar_Rbbxx_dX;
bool modelPar_given_Rbbxx;
double modelPar_Vgbb;
double d_modelPar_Vgbb_dX;
bool modelPar_given_Vgbb;
double modelPar_Jsr;
double d_modelPar_Jsr_dX;
bool modelPar_given_Jsr;
double modelPar_nr;
double d_modelPar_nr_dX;
bool modelPar_given_nr;
double modelPar_Vgr;
double d_modelPar_Vgr_dX;
bool modelPar_given_Vgr;
double modelPar_XCjc;
double d_modelPar_XCjc_dX;
bool modelPar_given_XCjc;
double modelPar_Jsc;
double d_modelPar_Jsc_dX;
bool modelPar_given_Jsc;
double modelPar_nc;
double d_modelPar_nc_dX;
bool modelPar_given_nc;
double modelPar_Rcxx;
double d_modelPar_Rcxx_dX;
bool modelPar_given_Rcxx;
double modelPar_Vgc;
double d_modelPar_Vgc_dX;
bool modelPar_given_Vgc;
double modelPar_Bf;
double d_modelPar_Bf_dX;
bool modelPar_given_Bf;
double modelPar_kBeta;
double d_modelPar_kBeta_dX;
bool modelPar_given_kBeta;
double modelPar_Br;
double d_modelPar_Br_dX;
bool modelPar_given_Br;
double modelPar_VAF;
double d_modelPar_VAF_dX;
bool modelPar_given_VAF;
double modelPar_VAR;
double d_modelPar_VAR_dX;
bool modelPar_given_VAR;
double modelPar_IKF;
double d_modelPar_IKF_dX;
bool modelPar_given_IKF;
double modelPar_IKR;
double d_modelPar_IKR_dX;
bool modelPar_given_IKR;
double modelPar_Mc;
double d_modelPar_Mc_dX;
bool modelPar_given_Mc;
double modelPar_BVceo;
double d_modelPar_BVceo_dX;
bool modelPar_given_BVceo;
double modelPar_kc;
double d_modelPar_kc_dX;
bool modelPar_given_kc;
double modelPar_BVebo;
double d_modelPar_BVebo_dX;
bool modelPar_given_BVebo;
double modelPar_Tr;
double d_modelPar_Tr_dX;
bool modelPar_given_Tr;
double modelPar_Trx;
double d_modelPar_Trx_dX;
bool modelPar_given_Trx;
double modelPar_Tf;
double d_modelPar_Tf_dX;
bool modelPar_given_Tf;
double modelPar_Tft;
double d_modelPar_Tft_dX;
bool modelPar_given_Tft;
double modelPar_Thcs;
double d_modelPar_Thcs_dX;
bool modelPar_given_Thcs;
double modelPar_Ahc;
double d_modelPar_Ahc_dX;
bool modelPar_given_Ahc;
double modelPar_Cje;
double d_modelPar_Cje_dX;
bool modelPar_given_Cje;
double modelPar_mje;
double d_modelPar_mje_dX;
bool modelPar_given_mje;
double modelPar_Vje;
double d_modelPar_Vje_dX;
bool modelPar_given_Vje;
double modelPar_Cjc;
double d_modelPar_Cjc_dX;
bool modelPar_given_Cjc;
double modelPar_mjc;
double d_modelPar_mjc_dX;
bool modelPar_given_mjc;
double modelPar_Vjc;
double d_modelPar_Vjc_dX;
bool modelPar_given_Vjc;
double modelPar_kjc;
double d_modelPar_kjc_dX;
bool modelPar_given_kjc;
double modelPar_Cmin;
double d_modelPar_Cmin_dX;
bool modelPar_given_Cmin;
double modelPar_J0;
double d_modelPar_J0_dX;
bool modelPar_given_J0;
double modelPar_XJ0;
double d_modelPar_XJ0_dX;
bool modelPar_given_XJ0;
double modelPar_Rci0;
double d_modelPar_Rci0_dX;
bool modelPar_given_Rci0;
double modelPar_Jk;
double d_modelPar_Jk_dX;
bool modelPar_given_Jk;
double modelPar_RJk;
double d_modelPar_RJk_dX;
bool modelPar_given_RJk;
double modelPar_Vces;
double d_modelPar_Vces_dX;
bool modelPar_given_Vces;
double modelPar_Rc;
double d_modelPar_Rc_dX;
bool modelPar_given_Rc;
double modelPar_Re;
double d_modelPar_Re_dX;
bool modelPar_given_Re;
double modelPar_Rb;
double d_modelPar_Rb_dX;
bool modelPar_given_Rb;
double modelPar_Rb2;
double d_modelPar_Rb2_dX;
bool modelPar_given_Rb2;
double modelPar_Lc;
double d_modelPar_Lc_dX;
bool modelPar_given_Lc;
double modelPar_Le;
double d_modelPar_Le_dX;
bool modelPar_given_Le;
double modelPar_Lb;
double d_modelPar_Lb_dX;
bool modelPar_given_Lb;
double modelPar_Cq;
double d_modelPar_Cq_dX;
bool modelPar_given_Cq;
double modelPar_Cpb;
double d_modelPar_Cpb_dX;
bool modelPar_given_Cpb;
double modelPar_Cpc;
double d_modelPar_Cpc_dX;
bool modelPar_given_Cpc;
double modelPar_Tnom;
double d_modelPar_Tnom_dX;
bool modelPar_given_Tnom;
// non-reals (including hidden)
int modelPar_Mode;
bool modelPar_given_Mode;
int modelPar_Noise;
bool modelPar_given_Noise;
int modelPar_Debug;
bool modelPar_given_Debug;
int modelPar_DebugPlus;
bool modelPar_given_DebugPlus;
int modelPar_N;
bool modelPar_given_N;
int modelPar_dtype;
};



//-----------------------------------------------------------------------------
// Free functions used by sensitivity
//
//-----------------------------------------------------------------------------
void evaluateModelEquations(
std::vector <double> & probeVars,
// probe constants
const int admsProbeID_V_t_ti,
const int admsProbeID_V_b_c,
const int admsProbeID_V_c_GND,
const int admsProbeID_V_b_GND,
const int admsProbeID_V_cx_bii,
const int admsProbeID_V_exx_bii,
const int admsProbeID_V_ex_bii,
const int admsProbeID_V_bii_bi,
const int admsProbeID_I_c_ci,
const int admsProbeID_I_e_ei,
const int admsProbeID_I_b_bi,
const int admsProbeID_V_ti_GND,
const int admsProbeID_V_ci_ei,
const int admsProbeID_V_exx_ei,
const int admsProbeID_V_cx_ci,
const int admsProbeID_V_ex_ei,
const int admsProbeID_V_bii_ei,
const int admsProbeID_V_bii_ci,
const int admsProbeID_V_bi_ci,
// node constants
const int admsNodeID_c,
const int admsNodeID_b,
const int admsNodeID_e,
const int admsNodeID_t,
const int admsNodeID_ei,
const int admsNodeID_bi,
const int admsNodeID_bii,
const int admsNodeID_ci,
const int admsNodeID_ti,
const int admsNodeID_ex,
const int admsNodeID_exx,
const int admsNodeID_cx,
const int admsBRA_ID_b_bi,
const int admsBRA_ID_e_ei,
const int admsBRA_ID_c_ci,
instanceSensStruct & instanceStruct,
modelSensStruct & modelStruct,
// basic variables
 double admsTemperature, double adms_vt_nom, double ADMSgmin_arg, std::vector <double> & d_staticContributions_dX, std::vector <double> & d_dynamicContributions_dX, const Instance & theInstance);

void evaluateInitialInstance(
instanceSensStruct & instanceStruct,
modelSensStruct & modelStruct,
 double admsTemperature,double adms_vt_nom, double ADMSgmin_arg, const Instance & theInstance);

void evaluateInitialModel(
modelSensStruct & modelStruct,
 double admsTemperature, double ADMSgmin_arg, const Instance & theInstance);

#endif // Xyce_ADMS_SENSITIVITIES


// Limited exponential --- NOT what verilog LRM says, but what qucs,
// ng-spice, and zspice do.

template <typename T>
T limexp(const T &x)
{
  if ((x) < 80.0)
  return (exp(x));
  else
  return (exp(80.0)*(x-79.0));
}


struct Traits: public DeviceTraits<Model, Instance, BJT::Traits>
{
  static const char *name() {return "FBH HBT_X v2.1";}
  static const char *deviceTypeName() {return "Q level 23";}

  static int numNodes() {return 4;}


  static bool modelRequired() {return true;}
  static bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

//-----------------------------------------------------------------------------
// Class         : Instance

//
// Purpose       : This class represents a single instance  of the
//                 device.  It mainly contains indices and pointers into
//                 the matrix equation (see the resistor instance class for
//                 more details).
//
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
#ifdef Xyce_ADMS_SENSITIVITIES
  friend class InstanceSensitivity;
  friend class ModelSensitivity;
#endif // Xyce_ADMS_SENSITIVITIES
  friend struct Traits;

  public:
    Instance(
      const Configuration &       configuration,
      const InstanceBlock &       instance_block,
      Model &                     model,
      const FactoryBlock &        factory_block);

    ~Instance();

private:
    Instance(const Instance &);
    Instance &operator=(const Instance &);

public:
    void registerLIDs( const LocalIdVector & intLIDVecRef,
                       const LocalIdVector & extLIDVecRef );
    void registerStoreLIDs( const LocalIdVector & stoLIDVecRef );
    void setupPointers();

    void loadNodeSymbols(Util::SymbolTable &symbol_table) const;

    const JacobianStamp & jacobianStamp() const;
    void registerJacLIDs( const JacobianStamp & jacLIDVec );

    void registerBranchDataLIDs(const std::vector<int> & branchLIDVecRef);

    bool processParams();
    bool updateTemperature ( const double & temp = -999.0 );
    bool updateIntermediateVars ();
    bool updatePrimaryState ();
    bool updateSecondaryState ();

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();


  private:

  public:
    // iterator reference to the HBT_X model which owns this instance.
    // Getters and setters
    Model &getModel()
    {
      return model_;
    }

  private:

    Model & model_;   //< Owning Model
    // Begin verilog Instance Variables
    //   Instance Parameters
    double Temp;
    int N;
    double L;
    double W;
    //  Variables of global_instance scope
    // end verilog Instance Variables=====
    // Nodal LID Variables
    int li_c;
    int li_b;
    int li_e;
    int li_t;
    int li_ei;
    int li_bi;
    int li_bii;
    int li_ci;
    int li_ti;
    int li_ex;
    int li_exx;
    int li_cx;
    // end Nodal LID Variables
    // Branch LID Variables
    int li_BRA_b_bi;
    int li_BRA_e_ei;
    int li_BRA_c_ci;
    // end Branch LID Variables
    // Lead (branch) LID Variables
    int li_branch_ic;
    int li_branch_ib;
    int li_branch_ie;
    int li_branch_it;
    // end Lead (branch) LID Variables
    // Jacobian  pointers
    double * f_bi_Equ_ti_Node_Ptr;
    double * f_ci_Equ_ti_Node_Ptr;
    double * f_bi_Equ_bi_Node_Ptr;
    double * f_bi_Equ_ci_Node_Ptr;
    double * f_ci_Equ_bi_Node_Ptr;
    double * f_ci_Equ_ci_Node_Ptr;
    double * f_bii_Equ_ti_Node_Ptr;
    double * f_bii_Equ_bii_Node_Ptr;
    double * f_bii_Equ_ci_Node_Ptr;
    double * f_ci_Equ_bii_Node_Ptr;
    double * f_bii_Equ_ei_Node_Ptr;
    double * f_ci_Equ_ei_Node_Ptr;
    double * f_ei_Equ_bii_Node_Ptr;
    double * f_ei_Equ_ci_Node_Ptr;
    double * f_ei_Equ_ti_Node_Ptr;
    double * f_ei_Equ_ei_Node_Ptr;
    double * f_ex_Equ_ti_Node_Ptr;
    double * f_ex_Equ_ex_Node_Ptr;
    double * f_ex_Equ_ei_Node_Ptr;
    double * f_ei_Equ_ex_Node_Ptr;
    double * f_exx_Equ_ti_Node_Ptr;
    double * f_exx_Equ_exx_Node_Ptr;
    double * f_exx_Equ_ei_Node_Ptr;
    double * f_ei_Equ_exx_Node_Ptr;
    double * f_cx_Equ_ti_Node_Ptr;
    double * f_cx_Equ_cx_Node_Ptr;
    double * f_cx_Equ_ci_Node_Ptr;
    double * f_ci_Equ_cx_Node_Ptr;
    double * f_bii_Equ_bi_Node_Ptr;
    double * f_bi_Equ_bii_Node_Ptr;
    double * f_ex_Equ_bii_Node_Ptr;
    double * f_bii_Equ_ex_Node_Ptr;
    double * f_exx_Equ_bii_Node_Ptr;
    double * f_bii_Equ_exx_Node_Ptr;
    double * f_cx_Equ_bii_Node_Ptr;
    double * f_bii_Equ_cx_Node_Ptr;
    double * f_b_Equ_b_Node_Ptr;
    double * f_c_Equ_c_Node_Ptr;
    double * f_b_Equ_c_Node_Ptr;
    double * f_c_Equ_b_Node_Ptr;
    double * f_ti_Equ_bi_Node_Ptr;
    double * f_ti_Equ_ci_Node_Ptr;
    double * f_ti_Equ_ei_Node_Ptr;
    double * f_ti_Equ_bii_Node_Ptr;
    double * f_ti_Equ_ti_Node_Ptr;
    double * f_t_Equ_t_Node_Ptr;
    double * f_t_Equ_ti_Node_Ptr;
    double * f_ti_Equ_t_Node_Ptr;
    double * f_c_Equ_BRA_c_ci_Var_Ptr;
    double * f_b_Equ_BRA_b_bi_Var_Ptr;
    double * f_e_Equ_BRA_e_ei_Var_Ptr;
    double * f_ei_Equ_BRA_e_ei_Var_Ptr;
    double * f_bi_Equ_BRA_b_bi_Var_Ptr;
    double * f_ci_Equ_BRA_c_ci_Var_Ptr;
    double * f_BRA_b_bi_Equ_b_Node_Ptr;
    double * f_BRA_b_bi_Equ_bi_Node_Ptr;
    double * f_BRA_b_bi_Equ_BRA_b_bi_Var_Ptr;
    double * f_BRA_e_ei_Equ_e_Node_Ptr;
    double * f_BRA_e_ei_Equ_ei_Node_Ptr;
    double * f_BRA_e_ei_Equ_BRA_e_ei_Var_Ptr;
    double * f_BRA_c_ci_Equ_c_Node_Ptr;
    double * f_BRA_c_ci_Equ_ci_Node_Ptr;
    double * f_BRA_c_ci_Equ_BRA_c_ci_Var_Ptr;
    double * q_bi_Equ_ti_Node_Ptr;
    double * q_ci_Equ_ti_Node_Ptr;
    double * q_bi_Equ_bi_Node_Ptr;
    double * q_bi_Equ_ci_Node_Ptr;
    double * q_ci_Equ_bi_Node_Ptr;
    double * q_ci_Equ_ci_Node_Ptr;
    double * q_bii_Equ_ti_Node_Ptr;
    double * q_bii_Equ_bii_Node_Ptr;
    double * q_bii_Equ_ci_Node_Ptr;
    double * q_ci_Equ_bii_Node_Ptr;
    double * q_bii_Equ_ei_Node_Ptr;
    double * q_ci_Equ_ei_Node_Ptr;
    double * q_ei_Equ_bii_Node_Ptr;
    double * q_ei_Equ_ci_Node_Ptr;
    double * q_ei_Equ_ti_Node_Ptr;
    double * q_ei_Equ_ei_Node_Ptr;
    double * q_ex_Equ_ti_Node_Ptr;
    double * q_ex_Equ_ex_Node_Ptr;
    double * q_ex_Equ_ei_Node_Ptr;
    double * q_ei_Equ_ex_Node_Ptr;
    double * q_exx_Equ_ti_Node_Ptr;
    double * q_exx_Equ_exx_Node_Ptr;
    double * q_exx_Equ_ei_Node_Ptr;
    double * q_ei_Equ_exx_Node_Ptr;
    double * q_cx_Equ_ti_Node_Ptr;
    double * q_cx_Equ_cx_Node_Ptr;
    double * q_cx_Equ_ci_Node_Ptr;
    double * q_ci_Equ_cx_Node_Ptr;
    double * q_bii_Equ_bi_Node_Ptr;
    double * q_bi_Equ_bii_Node_Ptr;
    double * q_ex_Equ_bii_Node_Ptr;
    double * q_bii_Equ_ex_Node_Ptr;
    double * q_exx_Equ_bii_Node_Ptr;
    double * q_bii_Equ_exx_Node_Ptr;
    double * q_cx_Equ_bii_Node_Ptr;
    double * q_bii_Equ_cx_Node_Ptr;
    double * q_b_Equ_b_Node_Ptr;
    double * q_c_Equ_c_Node_Ptr;
    double * q_b_Equ_c_Node_Ptr;
    double * q_c_Equ_b_Node_Ptr;
    double * q_ti_Equ_bi_Node_Ptr;
    double * q_ti_Equ_ci_Node_Ptr;
    double * q_ti_Equ_ei_Node_Ptr;
    double * q_ti_Equ_bii_Node_Ptr;
    double * q_ti_Equ_ti_Node_Ptr;
    double * q_t_Equ_t_Node_Ptr;
    double * q_t_Equ_ti_Node_Ptr;
    double * q_ti_Equ_t_Node_Ptr;
    double * q_c_Equ_BRA_c_ci_Var_Ptr;
    double * q_b_Equ_BRA_b_bi_Var_Ptr;
    double * q_e_Equ_BRA_e_ei_Var_Ptr;
    double * q_ei_Equ_BRA_e_ei_Var_Ptr;
    double * q_bi_Equ_BRA_b_bi_Var_Ptr;
    double * q_ci_Equ_BRA_c_ci_Var_Ptr;
    double * q_BRA_b_bi_Equ_b_Node_Ptr;
    double * q_BRA_b_bi_Equ_bi_Node_Ptr;
    double * q_BRA_b_bi_Equ_BRA_b_bi_Var_Ptr;
    double * q_BRA_e_ei_Equ_e_Node_Ptr;
    double * q_BRA_e_ei_Equ_ei_Node_Ptr;
    double * q_BRA_e_ei_Equ_BRA_e_ei_Var_Ptr;
    double * q_BRA_c_ci_Equ_c_Node_Ptr;
    double * q_BRA_c_ci_Equ_ci_Node_Ptr;
    double * q_BRA_c_ci_Equ_BRA_c_ci_Var_Ptr;
    // Jacobian offsets
    int A_bi_Equ_ti_NodeOffset;
    int A_ci_Equ_ti_NodeOffset;
    int A_bi_Equ_bi_NodeOffset;
    int A_bi_Equ_ci_NodeOffset;
    int A_ci_Equ_bi_NodeOffset;
    int A_ci_Equ_ci_NodeOffset;
    int A_bii_Equ_ti_NodeOffset;
    int A_bii_Equ_bii_NodeOffset;
    int A_bii_Equ_ci_NodeOffset;
    int A_ci_Equ_bii_NodeOffset;
    int A_bii_Equ_ei_NodeOffset;
    int A_ci_Equ_ei_NodeOffset;
    int A_ei_Equ_bii_NodeOffset;
    int A_ei_Equ_ci_NodeOffset;
    int A_ei_Equ_ti_NodeOffset;
    int A_ei_Equ_ei_NodeOffset;
    int A_ex_Equ_ti_NodeOffset;
    int A_ex_Equ_ex_NodeOffset;
    int A_ex_Equ_ei_NodeOffset;
    int A_ei_Equ_ex_NodeOffset;
    int A_exx_Equ_ti_NodeOffset;
    int A_exx_Equ_exx_NodeOffset;
    int A_exx_Equ_ei_NodeOffset;
    int A_ei_Equ_exx_NodeOffset;
    int A_cx_Equ_ti_NodeOffset;
    int A_cx_Equ_cx_NodeOffset;
    int A_cx_Equ_ci_NodeOffset;
    int A_ci_Equ_cx_NodeOffset;
    int A_bii_Equ_bi_NodeOffset;
    int A_bi_Equ_bii_NodeOffset;
    int A_ex_Equ_bii_NodeOffset;
    int A_bii_Equ_ex_NodeOffset;
    int A_exx_Equ_bii_NodeOffset;
    int A_bii_Equ_exx_NodeOffset;
    int A_cx_Equ_bii_NodeOffset;
    int A_bii_Equ_cx_NodeOffset;
    int A_b_Equ_b_NodeOffset;
    int A_c_Equ_c_NodeOffset;
    int A_b_Equ_c_NodeOffset;
    int A_c_Equ_b_NodeOffset;
    int A_ti_Equ_bi_NodeOffset;
    int A_ti_Equ_ci_NodeOffset;
    int A_ti_Equ_ei_NodeOffset;
    int A_ti_Equ_bii_NodeOffset;
    int A_ti_Equ_ti_NodeOffset;
    int A_t_Equ_t_NodeOffset;
    int A_t_Equ_ti_NodeOffset;
    int A_ti_Equ_t_NodeOffset;
    int A_c_Equ_BRA_c_ci_Var_Offset;
    int A_b_Equ_BRA_b_bi_Var_Offset;
    int A_e_Equ_BRA_e_ei_Var_Offset;
    int A_ei_Equ_BRA_e_ei_Var_Offset;
    int A_bi_Equ_BRA_b_bi_Var_Offset;
    int A_ci_Equ_BRA_c_ci_Var_Offset;
    int A_BRA_b_bi_Equ_b_Node_Offset;
    int A_BRA_b_bi_Equ_bi_Node_Offset;
    int A_BRA_b_bi_Equ_BRA_b_bi_Var_Offset;
    int A_BRA_e_ei_Equ_e_Node_Offset;
    int A_BRA_e_ei_Equ_ei_Node_Offset;
    int A_BRA_e_ei_Equ_BRA_e_ei_Var_Offset;
    int A_BRA_c_ci_Equ_c_Node_Offset;
    int A_BRA_c_ci_Equ_ci_Node_Offset;
    int A_BRA_c_ci_Equ_BRA_c_ci_Var_Offset;
    // end of Jacobian and pointers
   // node numbers
    static const int admsNodeID_c = 0;
    static const int admsNodeID_b = 1;
    static const int admsNodeID_e = 2;
    static const int admsNodeID_t = 3;
    static const int admsNodeID_ei = 0+4;
    static const int admsNodeID_bi = 1+4;
    static const int admsNodeID_bii = 2+4;
    static const int admsNodeID_ci = 3+4;
    static const int admsNodeID_ti = 4+4;
    static const int admsNodeID_ex = 5+4;
    static const int admsNodeID_exx = 6+4;
    static const int admsNodeID_cx = 7+4;
    static const int admsNodeID_GND = -1;
   // end node numbers
   // Additional IDs for branch equations
    static const int admsBRA_ID_b_bi = 12;
    static const int admsBRA_ID_e_ei = 13;
    static const int admsBRA_ID_c_ci = 14;
   // end branch numbers
   // Probe numbers
    static const int admsProbeID_V_t_ti = 0;
    static const int admsProbeID_V_b_c = 1;
    static const int admsProbeID_V_c_GND = 2;
    static const int admsProbeID_V_b_GND = 3;
    static const int admsProbeID_V_cx_bii = 4;
    static const int admsProbeID_V_exx_bii = 5;
    static const int admsProbeID_V_ex_bii = 6;
    static const int admsProbeID_V_bii_bi = 7;
    static const int admsProbeID_I_c_ci = 8;
    static const int admsProbeID_I_e_ei = 9;
    static const int admsProbeID_I_b_bi = 10;
    static const int admsProbeID_V_ti_GND = 11;
    static const int admsProbeID_V_ci_ei = 12;
    static const int admsProbeID_V_exx_ei = 13;
    static const int admsProbeID_V_cx_ci = 14;
    static const int admsProbeID_V_ex_ei = 15;
    static const int admsProbeID_V_bii_ei = 16;
    static const int admsProbeID_V_bii_ci = 17;
    static const int admsProbeID_V_bi_ci = 18;
   // end probe numbers
   // Store LIDs
    int li_store_admsProbeID_V_bi_ci;
    int li_store_admsProbeID_V_bii_ci;
    int li_store_admsProbeID_V_bii_ei;
   // end store LIDs
   // Store LIDs for output vars
   // end store LIDs for output vars
 // Arrays to hold probes
 std::vector < double > probeVars;
 std::vector < std::vector < double > > d_probeVars;
 // Arrays to hold contributions
 // dynamic contributions are differentiated w.r.t time
 std::vector < double > staticContributions;
 std::vector < std::vector < double > > d_staticContributions;
 std::vector < double > dynamicContributions;
 std::vector < std::vector < double > > d_dynamicContributions;

    // This array stores the differences between original and limited variables.
    std::vector<double> probeDiffs;
    // These store the Jdxp's for F and Q, respectively
    std::vector<double> Jdxp_static;
    std::vector<double> Jdxp_dynamic;

    // this is what we'll use when any model uses $temperature.  We'll
    // set it in updateTemperature, and initialize it to whatever
    // is in devOptions when the instance is constructed.
    double admsTemperature;

    // vt at $temperature;
    double adms_vt_nom;


    // This one is for the annoying bogus "XyceADMSInstTemp" parameter
    // that we need so we can set it from the device manager when there's no
    // "TEMP" parameter to use
    double admsInstTemp;


    static JacobianStamp jacStamp;
    static IdVector nodeMap;
    static PairMap pairToJacStampMap;

    // These instance-owned vectors are for storage of lead current data
    std::vector<double> leadCurrentF;
    std::vector<double> leadCurrentQ;


    };



namespace AnalogFunctions
{

      // Analog Function exp_soft
double exp_soft(double x);
// Derivative of Analog Function exp_soft
double d_exp_soft(double x  , double d_x  );
// Evaluator class for Analog Function exp_soft
class exp_softEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_x;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  exp_softEvaluator(double x);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  x);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double x  , double d_x);
private:
  returnType exp_softReturn_;
  returnType evaluator_(double x);
};


      // Analog Function Vt
double Vt(double U, double Ud);
// Derivative of Analog Function Vt
double d_Vt(double U , double Ud  , double d_U  , double d_Ud  );
// Evaluator class for Analog Function Vt
class VtEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_U;
     double deriv_WRT_Ud;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  VtEvaluator(double U, double Ud);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  U, double  Ud);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double U , double Ud  , double d_U, double d_Ud);
private:
  returnType VtReturn_;
  returnType evaluator_(double U, double Ud);
};


      // Analog Function diode
double diode(double U, double Is, double Ug, double N, double AREA, double TJ, double TNOM);
// Derivative of Analog Function diode
double d_diode(double U , double Is , double Ug , double N , double AREA , double TJ , double TNOM  , double d_U  , double d_Is  , double d_Ug  , double d_N  , double d_AREA  , double d_TJ  , double d_TNOM  );
// Evaluator class for Analog Function diode
class diodeEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_U;
     double deriv_WRT_Is;
     double deriv_WRT_Ug;
     double deriv_WRT_N;
     double deriv_WRT_AREA;
     double deriv_WRT_TJ;
     double deriv_WRT_TNOM;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  diodeEvaluator(double U, double Is, double Ug, double N, double AREA, double TJ, double TNOM);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  U, double  Is, double  Ug, double  N, double  AREA, double  TJ, double  TNOM);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double U , double Is , double Ug , double N , double AREA , double TJ , double TNOM  , double d_U, double d_Is, double d_Ug, double d_N, double d_AREA, double d_TJ, double d_TNOM);
private:
  returnType diodeReturn_;
  returnType evaluator_(double U, double Is, double Ug, double N, double AREA, double TJ, double TNOM);
};


      // Analog Function MM
double MM(double VBCI, double VCBO, double MC, double VCBLIN, double BF, double KC);
// Derivative of Analog Function MM
double d_MM(double VBCI , double VCBO , double MC , double VCBLIN , double BF , double KC  , double d_VBCI  , double d_VCBO  , double d_MC  , double d_VCBLIN  , double d_BF  , double d_KC  );
// Evaluator class for Analog Function MM
class MMEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_VBCI;
     double deriv_WRT_VCBO;
     double deriv_WRT_MC;
     double deriv_WRT_VCBLIN;
     double deriv_WRT_BF;
     double deriv_WRT_KC;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  MMEvaluator(double VBCI, double VCBO, double MC, double VCBLIN, double BF, double KC);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  VBCI, double  VCBO, double  MC, double  VCBLIN, double  BF, double  KC);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double VBCI , double VCBO , double MC , double VCBLIN , double BF , double KC  , double d_VBCI, double d_VCBO, double d_MC, double d_VCBLIN, double d_BF, double d_KC);
private:
  returnType MMReturn_;
  returnType evaluator_(double VBCI, double VCBO, double MC, double VCBLIN, double BF, double KC);
};


      // Analog Function charge
double charge(double U, double C0, double Ud, double m, double Area);
// Derivative of Analog Function charge
double d_charge(double U , double C0 , double Ud , double m , double Area  , double d_U  , double d_C0  , double d_Ud  , double d_m  , double d_Area  );
// Evaluator class for Analog Function charge
class chargeEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_U;
     double deriv_WRT_C0;
     double deriv_WRT_Ud;
     double deriv_WRT_m;
     double deriv_WRT_Area;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  chargeEvaluator(double U, double C0, double Ud, double m, double Area);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  U, double  C0, double  Ud, double  m, double  Area);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double U , double C0 , double Ud , double m , double Area  , double d_U, double d_C0, double d_Ud, double d_m, double d_Area);
private:
  returnType chargeReturn_;
  returnType evaluator_(double U, double C0, double Ud, double m, double Area);
};


      // Analog Function Vceff
double Vceff(double U, double VCES);
// Derivative of Analog Function Vceff
double d_Vceff(double U , double VCES  , double d_U  , double d_VCES  );
// Evaluator class for Analog Function Vceff
class VceffEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_U;
     double deriv_WRT_VCES;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  VceffEvaluator(double U, double VCES);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  U, double  VCES);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double U , double VCES  , double d_U, double d_VCES);
private:
  returnType VceffReturn_;
  returnType evaluator_(double U, double VCES);
};


      // Analog Function ICK
double ICK(double U, double RCI0, double VLIM, double InvVPT, double VCES);
// Derivative of Analog Function ICK
double d_ICK(double U , double RCI0 , double VLIM , double InvVPT , double VCES  , double d_U  , double d_RCI0  , double d_VLIM  , double d_InvVPT  , double d_VCES  );
// Evaluator class for Analog Function ICK
class ICKEvaluator
{
  struct returnType
  {
     double value;
     double deriv_WRT_U;
     double deriv_WRT_RCI0;
     double deriv_WRT_VLIM;
     double deriv_WRT_InvVPT;
     double deriv_WRT_VCES;
  };
public:
  // constructor takes all same arguments as regular templated function,
  // even though it doesn't USE the output args
  ICKEvaluator(double U, double RCI0, double VLIM, double InvVPT, double VCES);
  // function that returns the precomputed values.  This, too, takes
  // all the same arguments as the regular function, though it ONLY
  // uses the output arguments
  double getValues(double  U, double  RCI0, double  VLIM, double  InvVPT, double  VCES);
  // function that returns the total derivative of the function and its
  // output arguments with respect to some variable.  We pass in the
  // normal arguments(none of which are used) and the derivatives of those
  // arguments with respect to the desired variable.  We compute the
  // derivatives using the chain rule and our precomputed derivatives
  // with respect to input variables
double getDerivs(double U , double RCI0 , double VLIM , double InvVPT , double VCES  , double d_U, double d_RCI0, double d_VLIM, double d_InvVPT, double d_VCES);
private:
  returnType ICKReturn_;
  returnType evaluator_(double U, double RCI0, double VLIM, double InvVPT, double VCES);
};

}


//-----------------------------------------------------------------------------
// Class         : Model

// Purpose       :
// Special Notes :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
    typedef std::vector<Instance *> InstanceVector;

    friend class ParametricData<Model>;
    friend class Instance;
#ifdef Xyce_ADMS_SENSITIVITIES
    friend class InstanceSensitivity;
    friend class ModelSensitivity;
#endif // Xyce_ADMS_SENSITIVITIES
    friend struct Traits;

  public:
    Model(
      const Configuration &       configuration,
      const ModelBlock &          model_block,
      const FactoryBlock &        factory_block);

    ~Model();

private:
    Model(const Model &);
    Model &operator=(const Model &);

public:
    virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    virtual std::ostream &printOutInstances(std::ostream &os) const;
    bool processParams();
    bool processInstanceParams();

  private:

  public:
    void addInstance(Instance *instance)
    {
      instanceContainer.push_back(instance);
    }

    void setupBaseInstanceContainer()
    {
      std::vector<Instance*>::iterator iter = instanceContainer.begin();
      std::vector<Instance*>::iterator end   = instanceContainer.end();
      for ( ; iter!=end; ++iter)
      {
      Xyce::Device::DeviceModel::baseInstanceContainer.push_back( static_cast<Xyce::Device::DeviceInstance *>(*iter) );
    }
  }

  private:
    std::vector<Instance*> instanceContainer;

  private:

    // This one is for the annoying bogus "XyceADMSInstTemp" parameter
    // that we need so we can set it from the device manager when there's no
    // "TEMP" model parameter to use
    double admsModTemp;
// Begin verilog Model Variables
//   Model Parameters
    int Mode;
    int Noise;
    int Debug;
    int DebugPlus;
    double Temp;
    double Rth;
    double Cth;
    int N;
    double L;
    double W;
    double Jsf;
    double nf;
    double Vg;
    double Jse;
    double ne;
    double Rbxx;
    double Vgb;
    double Jsee;
    double nee;
    double Rbbxx;
    double Vgbb;
    double Jsr;
    double nr;
    double Vgr;
    double XCjc;
    double Jsc;
    double nc;
    double Rcxx;
    double Vgc;
    double Bf;
    double kBeta;
    double Br;
    double VAF;
    double VAR;
    double IKF;
    double IKR;
    double Mc;
    double BVceo;
    double kc;
    double BVebo;
    double Tr;
    double Trx;
    double Tf;
    double Tft;
    double Thcs;
    double Ahc;
    double Cje;
    double mje;
    double Vje;
    double Cjc;
    double mjc;
    double Vjc;
    double kjc;
    double Cmin;
    double J0;
    double XJ0;
    double Rci0;
    double Jk;
    double RJk;
    double Vces;
    double Rc;
    double Re;
    double Rb;
    double Rb2;
    double Lc;
    double Le;
    double Lb;
    double Cq;
    double Cpb;
    double Cpc;
    double Tnom;
    int dtype;
    //  Variables of global_model scope
    // end verilog model variables=====
};

void registerDevice(const DeviceCountMap& deviceMap = DeviceCountMap(),
                    const std::set<int>& levelSet = std::set<int>());

} // namespace ADMSHBT_X
} // namespace Device
} // namespace Xyce
#endif //Xyce_N_DEV_ADMSHBT_X_h
