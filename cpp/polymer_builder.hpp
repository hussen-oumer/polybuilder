// polybuilder / polymer_builder.hpp  --  THE CONNECTOR ENGINE (do not edit)
// ==========================================================================
// C++ mirror of the Python engine. Same algorithm: open the terminal vinyl
// C=C, tag the two backbone carbons with dummy attachment points, and join
// tail(:2) of one unit to head(:1) of the next. Side chains (cationic centre,
// CH2 bridge, sulfonate, ...) are never touched.
//
// Build: needs RDKit C++ (>= 2022) and CMake. See CMakeLists.txt.
#pragma once

#include <string>
#include <vector>
#include <memory>

#include <GraphMol/RWMol.h>
#include <GraphMol/ROMol.h>

namespace polybuilder {

// internal atom-map tags for the two open chain ends and an initiator cap
constexpr int HEAD = 1;   // the =CH2 carbon
constexpr int TAIL = 2;   // the substituted backbone carbon
constexpr int CAP  = 9;   // the initiator attachment dummy

// ------------------------------------------------------------------ Monomer
class Monomer {
 public:
  Monomer(std::string name, std::string smiles);
  virtual ~Monomer() = default;

  // "vinyl", overridden by subclasses only to describe the chemistry
  virtual std::string kind() const { return "vinyl"; }

  // build one repeat unit with HEAD/TAIL dummy attachment points
  std::unique_ptr<RDKit::RWMol> repeatUnit() const;

  const std::string& name()   const { return name_; }
  const std::string& smiles() const { return smiles_; }

 protected:
  std::string name_;
  std::string smiles_;
  std::shared_ptr<RDKit::ROMol> mol_;
};

// subclasses: inheritance carries the chemical meaning, not new build logic
class ZwitterionicMonomer : public Monomer {
 public: using Monomer::Monomer;
  std::string kind() const override { return "zwitterionic"; }
};
class DicationicMonomer : public Monomer {
 public: using Monomer::Monomer;
  std::string kind() const override { return "dicationic"; }
};
class NeutralMonomer : public Monomer {
 public: using Monomer::Monomer;
  std::string kind() const override { return "neutral"; }
};

// ---------------------------------------------------------------- Initiator
class Initiator {
 public:
  // capSmiles must contain exactly one [*] dummy, e.g. "[*]OS(=O)(=O)[O-]"
  Initiator(std::string name, std::string capSmiles, double percent = 1.0);
  const RDKit::ROMol& cap() const { return *cap_; }
  const std::string& name() const { return name_; }
  double percent() const { return percent_; }
 private:
  std::string name_;
  double percent_;
  std::shared_ptr<RDKit::ROMol> cap_;
};

// ------------------------------------------------------------ PolymerBuilder
class PolymerBuilder {
 public:
  // core builders (return an owning RWMol with two open ends still tagged)
  static std::unique_ptr<RDKit::RWMol> homopolymer(const Monomer& m, int n);
  static std::unique_ptr<RDKit::RWMol> copolymer(
      const std::vector<const Monomer*>& sequence);
  static std::unique_ptr<RDKit::RWMol> alternating(
      const Monomer& a, const Monomer& b, int nPairs);
  static std::unique_ptr<RDKit::RWMol> block(
      const Monomer& a, int na, const Monomer& b, int nb);

  // end handling
  static std::unique_ptr<RDKit::RWMol> capHydrogen(
      std::unique_ptr<RDKit::RWMol> chain);
  static std::unique_ptr<RDKit::RWMol> addInitiator(
      std::unique_ptr<RDKit::RWMol> chain, const Initiator& init,
      const std::string& ends = "both");

  // lab helper: initiator wt% -> rough number-average DP
  static int estimateDP(double monomerMW, double initiatorMW,
                        double initiatorWtPercent,
                        double conversion = 0.99,
                        double radicalEfficiency = 0.5);

  static std::string toSmiles(const RDKit::ROMol& mol);

 private:
  static std::unique_ptr<RDKit::RWMol> connect(
      std::unique_ptr<RDKit::RWMol> chain, std::unique_ptr<RDKit::RWMol> unit);
  static std::unique_ptr<RDKit::RWMol> attachCap(
      std::unique_ptr<RDKit::RWMol> chain, int tag, const RDKit::ROMol& cap);
  static int firstNeighbor(const RDKit::RWMol& m, int idx);
};

}  // namespace polybuilder
