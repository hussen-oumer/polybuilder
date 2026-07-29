// polybuilder / polymer_builder.cpp  --  THE CONNECTOR ENGINE (do not edit)
// ==========================================================================
#include "polymer_builder.hpp"

#include <stdexcept>
#include <algorithm>

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>  // combineMols
#include <GraphMol/MolOps.h>

namespace polybuilder {

// terminal vinyl:  CH2=C<
static RDKit::ROMol* vinylQuery() {
  static RDKit::ROMol* q = RDKit::SmartsToMol("[CH2]=[CX3]");
  return q;
}

// ------------------------------------------------------------------ Monomer
Monomer::Monomer(std::string name, std::string smiles)
    : name_(std::move(name)), smiles_(std::move(smiles)) {
  RDKit::ROMol* m = RDKit::SmilesToMol(smiles_);
  if (!m) throw std::runtime_error("Invalid SMILES for " + name_);
  mol_.reset(m);
}

std::unique_ptr<RDKit::RWMol> Monomer::repeatUnit() const {
  std::vector<RDKit::MatchVectType> matches;
  if (!RDKit::SubstructMatch(*mol_, *vinylQuery(), matches) || matches.empty())
    throw std::runtime_error("No polymerisable C=C found in " + name_);
  int a = matches[0][0].second;   // =CH2 carbon  -> HEAD
  int b = matches[0][1].second;   // substituted C -> TAIL

  auto unit = std::make_unique<RDKit::RWMol>(*mol_);
  unit->getBondBetweenAtoms(a, b)->setBondType(RDKit::Bond::SINGLE);

  int da = unit->addAtom(new RDKit::Atom(0), true, true);
  unit->getAtomWithIdx(da)->setAtomMapNum(HEAD);
  int db = unit->addAtom(new RDKit::Atom(0), true, true);
  unit->getAtomWithIdx(db)->setAtomMapNum(TAIL);
  unit->addBond(a, da, RDKit::Bond::SINGLE);
  unit->addBond(b, db, RDKit::Bond::SINGLE);

  RDKit::MolOps::sanitizeMol(*unit);
  return unit;
}

// ---------------------------------------------------------------- Initiator
Initiator::Initiator(std::string name, std::string capSmiles, double percent)
    : name_(std::move(name)), percent_(percent) {
  RDKit::RWMol* c = RDKit::SmilesToMol(capSmiles);
  if (!c) throw std::runtime_error("Invalid initiator SMILES for " + name_);
  bool tagged = false;
  for (auto at : c->atoms())
    if (at->getAtomicNum() == 0) { at->setAtomMapNum(CAP); tagged = true; }
  if (!tagged)
    throw std::runtime_error("Initiator " + name_ + " needs one [*] dummy");
  cap_.reset(c);
}

// ------------------------------------------------------------ PolymerBuilder
int PolymerBuilder::firstNeighbor(const RDKit::RWMol& m, int idx) {
  const RDKit::Atom* a = m.getAtomWithIdx(idx);
  RDKit::ROMol::ADJ_ITER nbr, end;
  boost::tie(nbr, end) = m.getAtomNeighbors(a);
  return static_cast<int>(*nbr);
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::connect(
    std::unique_ptr<RDKit::RWMol> chain, std::unique_ptr<RDKit::RWMol> unit) {
  // chain first, unit second -> chain atoms keep the lower indices
  std::unique_ptr<RDKit::ROMol> combo(RDKit::combineMols(*chain, *unit));
  auto rw = std::make_unique<RDKit::RWMol>(*combo);

  std::vector<int> tails, heads;
  for (auto at : rw->atoms()) {
    if (at->getAtomMapNum() == TAIL) tails.push_back(at->getIdx());
    if (at->getAtomMapNum() == HEAD) heads.push_back(at->getIdx());
  }
  std::sort(tails.begin(), tails.end());
  std::sort(heads.begin(), heads.end());
  int tailDummy = tails.front();   // chain tail (lowest index)
  int headDummy = heads.back();    // unit head (highest index)

  int cTail = firstNeighbor(*rw, tailDummy);
  int cHead = firstNeighbor(*rw, headDummy);
  rw->addBond(cTail, cHead, RDKit::Bond::SINGLE);

  // remove the two consumed dummies, higher index first
  std::vector<int> rm{tailDummy, headDummy};
  std::sort(rm.rbegin(), rm.rend());
  for (int idx : rm) rw->removeAtom(idx);

  RDKit::MolOps::sanitizeMol(*rw);
  return rw;
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::homopolymer(
    const Monomer& m, int n) {
  auto chain = m.repeatUnit();
  for (int i = 1; i < n; ++i) chain = connect(std::move(chain), m.repeatUnit());
  return chain;
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::copolymer(
    const std::vector<const Monomer*>& seq) {
  if (seq.empty()) throw std::runtime_error("empty sequence");
  auto chain = seq[0]->repeatUnit();
  for (size_t i = 1; i < seq.size(); ++i)
    chain = connect(std::move(chain), seq[i]->repeatUnit());
  return chain;
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::alternating(
    const Monomer& a, const Monomer& b, int nPairs) {
  std::vector<const Monomer*> seq;
  for (int i = 0; i < nPairs; ++i) { seq.push_back(&a); seq.push_back(&b); }
  return copolymer(seq);
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::block(
    const Monomer& a, int na, const Monomer& b, int nb) {
  std::vector<const Monomer*> seq;
  for (int i = 0; i < na; ++i) seq.push_back(&a);
  for (int i = 0; i < nb; ++i) seq.push_back(&b);
  return copolymer(seq);
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::capHydrogen(
    std::unique_ptr<RDKit::RWMol> chain) {
  std::vector<int> dummies;
  for (auto at : chain->atoms())
    if (at->getAtomMapNum() == HEAD || at->getAtomMapNum() == TAIL)
      dummies.push_back(at->getIdx());
  std::sort(dummies.rbegin(), dummies.rend());
  for (int idx : dummies) chain->removeAtom(idx);
  RDKit::MolOps::sanitizeMol(*chain);
  return chain;
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::attachCap(
    std::unique_ptr<RDKit::RWMol> chain, int tag, const RDKit::ROMol& cap) {
  std::unique_ptr<RDKit::ROMol> combo(RDKit::combineMols(*chain, cap));
  auto rw = std::make_unique<RDKit::RWMol>(*combo);

  int endDummy = -1, capDummy = -1;
  for (auto at : rw->atoms()) {
    if (at->getAtomMapNum() == tag) endDummy = at->getIdx();
    if (at->getAtomMapNum() == CAP) capDummy = at->getIdx();
  }
  if (endDummy < 0) return rw;   // nothing to cap

  int cEnd = firstNeighbor(*rw, endDummy);
  int cCap = firstNeighbor(*rw, capDummy);
  rw->addBond(cEnd, cCap, RDKit::Bond::SINGLE);

  std::vector<int> rm{endDummy, capDummy};
  std::sort(rm.rbegin(), rm.rend());
  for (int idx : rm) rw->removeAtom(idx);

  RDKit::MolOps::sanitizeMol(*rw);
  return rw;
}

std::unique_ptr<RDKit::RWMol> PolymerBuilder::addInitiator(
    std::unique_ptr<RDKit::RWMol> chain, const Initiator& init,
    const std::string& ends) {
  if (ends == "both" || ends == "head")
    chain = attachCap(std::move(chain), HEAD, init.cap());
  if (ends == "both" || ends == "tail")
    chain = attachCap(std::move(chain), TAIL, init.cap());
  return capHydrogen(std::move(chain));
}

int PolymerBuilder::estimateDP(double monomerMW, double initiatorMW,
                               double initiatorWtPercent, double conversion,
                               double radicalEfficiency) {
  double molesRatio = (100.0 / initiatorWtPercent) * (initiatorMW / monomerMW);
  int dp = static_cast<int>(conversion * molesRatio / (2.0 * radicalEfficiency));
  return std::max(1, dp);
}

std::string PolymerBuilder::toSmiles(const RDKit::ROMol& mol) {
  return RDKit::MolToSmiles(mol);
}

}  // namespace polybuilder
