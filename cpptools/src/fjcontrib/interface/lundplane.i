%module lundplane
%{
  #include <fastjet/PseudoJet.hh>
  #include <fastjet/FunctionOfPseudoJet.hh>
  #include <LundPlane/LundGenerator.hh>
  #include <LundPlane/LundJSON.hh>
  #include <LundPlane/LundWithSecondary.hh>
  #include <LundPlane/SecondaryLund.hh>
%}
%include "std_string.i"
%include "std_vector.i"

// Process symbols in header
// %nodefaultctor Recluster;

namespace fastjet{
  namespace contrib{

// %include "LundPlane/LundGenerator.hh"
//----------------------------------------------------------------------

// class LundGenerator;

//----------------------------------------------------------------------
// %nodefaultctor LundDeclustering;
// %nodefaultdtor LundDeclustering;
class LundDeclustering {
public:
  const fastjet::PseudoJet & pair()  const;
  const fastjet::PseudoJet & harder() const;
  const fastjet::PseudoJet & softer() const;
  double m()         const;
  double Delta()     const;
  double z()         const;
  double kt()        const;
  double kappa()     const;
  double psi()       const;
  std::pair<double,double> const lund_coordinates() const;
  virtual ~LundDeclustering();
  LundDeclustering(const fastjet::PseudoJet& pair,
       const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2);
  // friend class fastjet::contrib::LundGenerator;
};

//----------------------------------------------------------------------
class LundGenerator : public fastjet::FunctionOfPseudoJet< std::vector<fastjet::contrib::LundDeclustering> > {
public:
  LundGenerator(fastjet::JetAlgorithm jet_alg = fastjet::Algorithm::cambridge_algorithm);
  LundGenerator(const fastjet::JetDefinition & jet_def);
  virtual ~LundGenerator();
  std::vector<fastjet::contrib::LundDeclustering> result(const fastjet::PseudoJet& jet) const;
  virtual std::string description() const;
};

// %include "LundPlane/LundJSON.hh"
//----------------------------------------------------------------------
// declaration of helper function
void lund_elements_to_json(std::ostream & ostr, const fastjet::contrib::LundDeclustering & d);

/// writes json to ostr for an individual declustering
std::ostream & lund_to_json(std::ostream & ostr, const fastjet::contrib::LundDeclustering & d);

/// writes json to ostr for a vector of declusterings
std::ostream & lund_to_json(std::ostream & ostr, const std::vector<fastjet::contrib::LundDeclustering> & d);

// helper function to write individual elements to json
void lund_elements_to_json(std::ostream & ostr, const fastjet::contrib::LundDeclustering & d);

// %include "LundPlane/LundWithSecondary.hh"
//----------------------------------------------------------------------
class LundWithSecondary {
public:
  /// LundWithSecondary constructor
  LundWithSecondary(fastjet::contrib::SecondaryLund * secondary_def = 0);
  /// LundWithSecondary constructor with jet alg
  LundWithSecondary(fastjet::JetAlgorithm jet_alg,
        fastjet::contrib::SecondaryLund * secondary_def = 0);
  /// LundWithSecondary constructor with jet def
  LundWithSecondary(const fastjet::JetDefinition & jet_def,
        fastjet::contrib::SecondaryLund * secondary_def = 0);
  /// destructor
  virtual ~LundWithSecondary();
  /// primary Lund declustering
  std::vector<fastjet::contrib::LundDeclustering> primary(const fastjet::PseudoJet& jet) const;
  /// secondary Lund declustering (slow)
  std::vector<fastjet::contrib::LundDeclustering> secondary(const fastjet::PseudoJet& jet) const;
  /// secondary Lund declustering with primary sequence as input
  std::vector<fastjet::contrib::LundDeclustering> secondary(
       const std::vector<fastjet::contrib::LundDeclustering> & declusts) const;
  /// return the index associated of the primary declustering that is to be
  /// used for the secondary plane.
  int secondary_index(const std::vector<fastjet::contrib::LundDeclustering> & declusts) const;
  /// description of the class
  std::string description() const;
};

// %include "LundPlane/SecondaryLund.hh"
//----------------------------------------------------------------------
class SecondaryLund {
public:
  /// SecondaryLund constructor
  SecondaryLund();
  /// destructor
  virtual ~SecondaryLund();
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<fastjet::contrib::LundDeclustering> & declusts);
  int operator()(const std::vector<fastjet::contrib::LundDeclustering> & declusts) const;
  /// description of the class
  virtual std::string description() const;
};

//----------------------------------------------------------------------
/// \class SecondaryLund_mMDT
/// Contains a definition for the leading emission using mMDTZ
class SecondaryLund_mMDT : public SecondaryLund {
public:
  /// SecondaryLund_mMDT constructor
  SecondaryLund_mMDT(double zcut = 0.025);
  /// destructor
  virtual ~SecondaryLund_mMDT();
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<fastjet::contrib::LundDeclustering> & declusts) const;
  /// description of the class
  virtual std::string description() const;
};

//----------------------------------------------------------------------
/// \class SecondaryLund_dotmMDT
/// Contains a definition for the leading emission using dotmMDT
class SecondaryLund_dotmMDT : public fastjet::contrib::SecondaryLund {
public:
  /// SecondaryLund_dotmMDT constructor
  SecondaryLund_dotmMDT(double zcut = 0.025);
  /// destructor
  virtual ~SecondaryLund_dotmMDT();
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<fastjet::contrib::LundDeclustering> & declusts) const;
  /// description of the class
  virtual std::string description() const;
private:
  /// zcut parameter
  double zcut_;
};

//----------------------------------------------------------------------
/// \class SecondaryLund_Mass
/// Contains a definition for the leading emission using mass
class SecondaryLund_Mass : public fastjet::contrib::SecondaryLund {
public:
  /// SecondaryLund_Mass constructor (default mass reference is W mass)
  SecondaryLund_Mass(double ref_mass = 80.4) : mref2_(ref_mass*ref_mass);
  /// destructor
  virtual ~SecondaryLund_Mass();
  /// returns the index of branch corresponding to the root of the secondary plane
  virtual int result(const std::vector<fastjet::contrib::LundDeclustering> & declusts);
  /// description of the class
  virtual std::string description() const;
};

  } // namespace contrib
} // namespace fastjet

namespace std
{
  %template(FastJetVLundDeclustering) std::vector<fastjet::contrib::LundDeclustering>;
}
