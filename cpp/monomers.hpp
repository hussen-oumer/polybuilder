// polybuilder / monomers.hpp  --  THE USER FILE (edit this one)
// =============================================================
// Add a new monomer or initiator by adding one line here. Nothing in this
// file knows how to build a polymer -- that stays in polymer_builder.*.
#pragma once

#include <map>
#include <memory>
#include <string>
#include "polymer_builder.hpp"

namespace polybuilder {

// Build the monomer library. Keep the polymerisable bond a terminal vinyl
// (C=C... or styrenic C=Cc1ccc...). The rest of the SMILES is the side chain.
inline std::map<std::string, std::shared_ptr<Monomer>> makeMonomers() {
  std::map<std::string, std::shared_ptr<Monomer>> lib;
  auto add = [&](const std::string& k, std::shared_ptr<Monomer> m) { lib[k] = m; };

  // zwitterionic sulfobetaines
  add("DMAPS", std::make_shared<ZwitterionicMonomer>(
      "DMAPS", "C=C(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"));
  add("M3295", std::make_shared<ZwitterionicMonomer>(
      "M3295", "C=C(C)C(=O)OCC[N+](C)(C)CCCCS(=O)(=O)[O-]"));
  add("A3367", std::make_shared<ZwitterionicMonomer>(
      "A3367", "C=CC(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"));
  add("A3361", std::make_shared<ZwitterionicMonomer>(
      "A3361", "C=CC(=O)NCCC[N+](C)(C)CCCS(=O)(=O)[O-]"));

  // dicationic DABCO (counter-ions added later in MD)
  add("DABCO", std::make_shared<DicationicMonomer>(
      "DABCO", "C=Cc1ccc(C[N+]23CC[N+](CCCC)(CC2)CC3)cc1"));

  // neutral comonomers
  add("MMA",   std::make_shared<NeutralMonomer>("MMA",   "C=C(C)C(=O)OC"));
  add("NIPAM", std::make_shared<NeutralMonomer>("NIPAM", "C=CC(=O)NC(C)C"));
  add("IBOA",  std::make_shared<NeutralMonomer>("IBOA",  "C=CC(=O)OC1CCC2(C)C(C)(C)C1C2"));
  return lib;
}

// Initiator library: cap fragment carries one [*], plus experimental wt%.
inline std::map<std::string, std::shared_ptr<Initiator>> makeInitiators() {
  std::map<std::string, std::shared_ptr<Initiator>> lib;
  lib["KPS"]  = std::make_shared<Initiator>("KPS",  "[*]OS(=O)(=O)[O-]", 1.0);
  lib["AIBA"] = std::make_shared<Initiator>("AIBA", "[*]C(N)=[NH2+]",    1.0);
  lib["TBHP"] = std::make_shared<Initiator>("TBHP", "[*]OC(C)(C)C",      1.0);
  return lib;
}

}  // namespace polybuilder
