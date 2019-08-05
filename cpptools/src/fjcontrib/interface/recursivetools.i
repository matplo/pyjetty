%module recursivetools
%{
  #include <fastjet/FunctionOfPseudoJet.hh>
  #include <fastjet/tools/Transformer.hh>
  #include <fastjet/ClusterSequence.hh>
  #include <fastjet/WrappedStructure.hh>
  #include <fastjet/tools/Transformer.hh>

  #include <RecursiveTools/Recluster.hh>
  #include <RecursiveTools/IteratedSoftDrop.hh>
  #include <RecursiveTools/RecursiveSymmetryCutBase.hh>
  #include <RecursiveTools/BottomUpSoftDrop.hh>
  #include <RecursiveTools/RecursiveSoftDrop.hh>
  #include <RecursiveTools/ModifiedMassDropTagger.hh>
  #include <RecursiveTools/SoftDrop.hh>
  #include <RecursiveTools/Util.hh>
%}

%include "std_string.i"
%include "std_vector.i"

// Process symbols in header
// %nodefaultctor Recluster;

namespace fastjet{
  namespace contrib{
    // %include "RecursiveTools/Recluster.hh"
    class Recluster : public fastjet::Transformer {
    public:
      Recluster(const fastjet::JetDefinition & subjet_def, bool single=true);
      Recluster(fastjet::JetAlgorithm subjet_alg, double subjet_radius, double subjet_extra,
                bool single=true);
      Recluster(fastjet::JetAlgorithm subjet_alg, double subjet_radius, bool single=true);
      Recluster(fastjet::JetAlgorithm subjet_alg, bool single=true);
      virtual ~Recluster();
      virtual fastjet::PseudoJet result(const fastjet::PseudoJet & jet) const;
      virtual std::string description() const;
      typedef CompositeJetStructure StructureType;
    };

    // %include "RecursiveTools/RecursiveSymmetryCutBase.hh"
    class RecursiveSymmetryCutBase : public fastjet::Transformer {
    public:
      enum SymmetryMeasure{scalar_z,   ///< \f$ \min(p_{ti}, p_{tj})/(p_{ti} + p_{tj}) \f$
                           vector_z,   ///< \f$ \min(p_{ti}, p_{tj})/p_{t(i+j)} \f$
                           y,          ///< \f$ \min(p_{ti}^2,p_{tj}^2) \Delta R_{ij}^2 / m_{ij}^2 \f$
                           theta_E,    ///< \f$ \min(E_i,E_j)/(E_i+E_j) \f$ with 3d angle (ee collisions)
                           cos_theta_E ///< \f$ \min(E_i,E_j)/(E_i+E_j) \f$ with
                                       ///  \f$ \sqrt{2[1-cos(theta)]}\f$ for angles (ee collisions)
      };

      enum RecursionChoice{larger_pt, ///< choose the subjet with larger \f$ p_t \f$
                           larger_mt, ///< choose the subjet with larger \f$ m_t \equiv (m^2+p_t^2)^{\frac12}] \f$
                           larger_m,  ///< choose the subjet with larger mass (deprecated)
                           larger_E   ///< choose the subjet with larger energy (meant for ee collisions)
      };

      RecursiveSymmetryCutBase(fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure = scalar_z,
                               double           mu_cut = std::numeric_limits<double>::infinity(),
                               fastjet::contrib::RecursiveSymmetryCutBase::RecursionChoice  recursion_choice = larger_pt,
                               const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0
                               );
      virtual ~RecursiveSymmetryCutBase();
      fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure symmetry_measure() const;
      double mu_cut() const;
      fastjet::contrib::RecursiveSymmetryCutBase::RecursionChoice recursion_choice() const;
      void set_input_jet_is_subtracted(bool is_subtracted);
      bool input_jet_is_subtracted() const;
      void set_subtractor(const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor_);
      const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor() const;
      void set_reclustering(bool do_reclustering=true, const Recluster *recluster=0);
      void set_grooming_mode(bool enable=true);
      void set_tagging_mode(bool enable=true);
      void set_verbose_structure(bool enable=true);
      bool has_verbose_structure() const;
      virtual fastjet::PseudoJet result(const fastjet::PseudoJet & j) const;
      virtual std::string description() const;
      double squared_geometric_distance(const fastjet::PseudoJet &j1,
                                        const fastjet::PseudoJet &j2) const;
    };

    class RecursiveSymmetryCutBase::StructureType : public fastjet::WrappedStructure {
    public:
      StructureType(const fastjet::PseudoJet & j);
      StructureType(const fastjet::PseudoJet & j, double delta_R_in, double symmetry_in, double mu_in=-1.0);
      double delta_R()  const;
      double thetag()   const;
      double symmetry() const;
      double zg()       const;
      double mu()       const;
      bool has_verbose() const;
      void set_verbose(bool value);
      bool has_substructure() const;
      int dropped_count(bool global=true);
      std::vector<double> dropped_delta_R(bool global=true) const;
      void set_dropped_delta_R(const std::vector<double> &v);
      std::vector<double> dropped_symmetry(bool global=true) const;
      void set_dropped_symmetry(const std::vector<double> &v);
      std::vector<double> dropped_mu(bool global=true) const;
      void set_dropped_mu(const std::vector<double> &v);
      double max_dropped_symmetry(bool global=true) const;
      std::vector<std::pair<double,double> > sorted_zg_and_thetag() const;
    };

    // %include "RecursiveTools/SoftDrop.hh"
    class SoftDrop : public fastjet::contrib::RecursiveSymmetryCutBase {
    public:
      SoftDrop(double beta,
               double symmetry_cut,
               double R0 = 1,
               const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      SoftDrop(double           beta,
               double           symmetry_cut,
               fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure,
               double           R0 = 1.0,
               double           mu_cut = std::numeric_limits<double>::infinity(),
               fastjet::contrib::RecursiveSoftDrop::RecursionChoice  recursion_choice = larger_pt,
               const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      virtual ~SoftDrop(){}
      double beta()         const;
      double symmetry_cut() const;
      double R0()           const;
    };

    struct SDinfo
    {
    public:
      SDinfo();
      ~SDinfo();
      double z;
      double dR;
      double mu;
    };

    fastjet::contrib::SDinfo get_SD_jet_info(const fastjet::PseudoJet &j);

    // %include "RecursiveTools/RecursiveSoftDrop.hh"
    class RecursiveSoftDrop : public fastjet::contrib::SoftDrop {
    public:
      RecursiveSoftDrop(double beta,
                        double symmetry_cut,
                        int n = -1,
                        double R0 = 1,
                        const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      RecursiveSoftDrop(double           beta,
                        double           symmetry_cut,
                        fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure,
                        int              n = -1,
                        double           R0 = 1.0,
                        double           mu_cut = std::numeric_limits<double>::infinity(),
                        fastjet::contrib::RecursiveSymmetryCutBase::RecursionChoice  recursion_choice = larger_pt,
                        const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      virtual ~RecursiveSoftDrop();
      int n() const;
      void set_defaults();
      void set_fixed_depth_mode(bool value=true);
      bool fixed_depth_mode() const;
      void set_dynamical_R0(bool value=true);
      bool use_dynamical_R0() const;
      void set_hardest_branch_only(bool value=true);
      bool use_hardest_branch_only() const;
      void set_min_deltaR_squared(double value=-1.0);
      double   min_deltaR_squared() const;
      virtual std::string description() const;
      virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const;
      fastjet::PseudoJet result_fixed_tags(const fastjet::PseudoJet &jet) const;
      fastjet::PseudoJet result_fixed_depth(const fastjet::PseudoJet &jet) const;
    };

    std::vector<fastjet::PseudoJet> recursive_soft_drop_prongs(const fastjet::PseudoJet & rsd_jet);

    // %include "RecursiveTools/IteratedSoftDrop.hh"
    class IteratedSoftDropInfo{
    public:
      IteratedSoftDropInfo();
      IteratedSoftDropInfo(std::vector<std::pair<double,double> > zg_thetag_in);
      const std::vector<std::pair<double,double> > &all_zg_thetag();
      const std::vector<std::pair<double,double> > & operator()();
      // const std::pair<double,double> & operator[](unsigned int i);
      double angularity(double alpha, double kappa=1.0) const;
      unsigned int multiplicity() const;
      unsigned int size() const;
    };

    %extend IteratedSoftDropInfo {
    const std::pair<double,double> __getitem__(unsigned int i) {
        return (*($self))[i];
    }
}

    class IteratedSoftDrop : public fastjet::FunctionOfPseudoJet<IteratedSoftDropInfo> {
    public:
      IteratedSoftDrop(double beta, double symmetry_cut, double angular_cut, double R0 = 1.0,
                       const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      IteratedSoftDrop(double  beta,
                       double  symmetry_cut,
                       fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure,
                       double  angular_cut,
                       double  R0 = 1.0,
                       double  mu_cut = std::numeric_limits<double>::infinity(),
                       fastjet::contrib::RecursiveSoftDrop::RecursionChoice  recursion_choice = fastjet::contrib::RecursiveSoftDrop::larger_pt,
                       const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0);
      virtual ~IteratedSoftDrop();
      void set_dynamical_R0(bool value=true);
      bool use_dynamical_R0() const;
      void set_subtractor(const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor_);
      const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor() const;
      IteratedSoftDropInfo result(const fastjet::PseudoJet& jet) const;
      void set_input_jet_is_subtracted(bool is_subtracted);
      bool input_jet_is_subtracted() const;
      void set_reclustering(bool do_reclustering=true, const Recluster *recluster=0);

      std::vector<std::pair<double,double> > all_zg_thetag(const fastjet::PseudoJet& jet) const;
      double angularity(const fastjet::PseudoJet& jet, double alpha, double kappa=1.0) const;
      double multiplicity(const fastjet::PseudoJet& jet) const;
      /// description of the class
      std::string description() const;
    };

    // %include "RecursiveTools/BottomUpSoftDrop.hh"
    class BottomUpSoftDropStructure : public fastjet::WrappedStructure{
    public:
      BottomUpSoftDropStructure(const fastjet::PseudoJet & result_jet);
      virtual std::string description();
      std::vector<fastjet::PseudoJet> rejected();
      std::vector<fastjet::PseudoJet> extra_jets() const;
      double beta() const;
      double symmetry_cut() const;
      double R0();
    };

    class BottomUpSoftDrop : public fastjet::Transformer {
    public:
      BottomUpSoftDrop(double beta, double symmetry_cut, double R0 = 1.0);
      BottomUpSoftDrop(const fastjet::JetAlgorithm jet_alg, double beta, double symmetry_cut,
           double R0 = 1.0);
      BottomUpSoftDrop(const fastjet::JetDefinition &jet_def, double beta, double symmetry_cut,
           double R0 = 1.0);
      virtual fastjet::PseudoJet result(const fastjet::PseudoJet &jet) const;
      virtual std::vector<fastjet::PseudoJet> global_grooming(const std::vector<fastjet::PseudoJet> & event) const;
      virtual std::string description() const;
      typedef BottomUpSoftDropStructure StructureType;
    };

    class BottomUpSoftDropRecombiner : public fastjet::JetDefinition::Recombiner {
    public:
      BottomUpSoftDropRecombiner(double beta, double symmetry_cut, double R0,
                                 const fastjet::JetDefinition::Recombiner *recombiner);
      virtual void recombine(const fastjet::PseudoJet &pa,
           const fastjet::PseudoJet &pb,
           fastjet::PseudoJet &pab) const;
      virtual std::string description();
      const std::vector<unsigned int> & rejected() const{ return _rejected;}
      void clear_rejected();
    };

    class BottomUpSoftDropPlugin : public fastjet::JetDefinition::Plugin {
    public:
      BottomUpSoftDropPlugin(const fastjet::JetDefinition &jet_def, double beta, double symmetry_cut,
                             double R0 = 1.0);
      virtual void run_clustering(fastjet::ClusterSequence &input_cs) const;
      virtual std::string description() const;
      virtual double R() const;
    };

    // %include "RecursiveTools/ModifiedMassDropTagger.hh"
    class ModifiedMassDropTagger : public RecursiveSymmetryCutBase {
    public:
      ModifiedMassDropTagger(double symmetry_cut,
                             const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0
                             );
      ModifiedMassDropTagger(double           symmetry_cut,
                             fastjet::contrib::RecursiveSymmetryCutBase::SymmetryMeasure  symmetry_measure,
                             double           mu_cut = std::numeric_limits<double>::infinity(),
                             fastjet::contrib::RecursiveSymmetryCutBase::RecursionChoice  recursion_choice = larger_pt,
                             const fastjet::FunctionOfPseudoJet<fastjet::PseudoJet> * subtractor = 0
                             );
      virtual ~ModifiedMassDropTagger();
      double symmetry_cut() const;
    };

  } // namespace contrib
} // namespace fastjet
